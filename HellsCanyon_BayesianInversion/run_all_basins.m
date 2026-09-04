function summary = run_all_basins(varargin)
% RUN_ALL_BASINS  Run the Hells Canyon inversion over several basins.
%
%   summary = run_all_basins()
%   summary = run_all_basins('parallel', true, 'n_postburn', 5e4)
%
% Basins are INDEPENDENT: each gets its own chain, its own t_capture, its
% own K/ksn/n.  Cross-basin agreement in t_capture is then a result rather
% than something imposed -- "five tributaries independently converge on
% ~2.2 Ma" is a much stronger statement than fitting a shared parameter.
% A hierarchical model with shared timing is the natural follow-up if the
% independent runs justify it.
%
% Because basins are independent, this is embarrassingly parallel.  With
% Parallel Computing Toolbox, 'parallel' runs them under parfor; without
% it, MATLAB simply executes parfor serially, so the same code works.
%
% Options (name/value, all optional):
%   'parallel'   logical, default true
%   'basins'     struct array to use instead of hc_basin_list()
%   'output_dir' overrides the default
%   plus any field of hc_basin_defaults (e.g. 'n_postburn', 'prop_mode'),
%   applied globally to every basin.
%
% Each basin is wrapped in try/catch and saved as it completes, so one bad
% basin cannot cost you the others.
%
% See also: hc_basin_defaults, hc_basin_list, run_hc_inversion,
%           hc_multibasin_summary

%% ------------------------------------------------------------- options
p = inputParser;
p.KeepUnmatched = true;
addParameter(p, 'parallel',   true,  @(x) islogical(x) || isnumeric(x));
addParameter(p, 'basins',     [],    @(x) isempty(x) || isstruct(x));
parse(p, varargin{:});
use_par  = logical(p.Results.parallel);
basins   = p.Results.basins;

defaults = hc_basin_defaults();
% Any leftover name/value pairs override the global defaults.
extra = p.Unmatched;
fn = fieldnames(extra);
for i = 1:numel(fn)
    if ~isfield(defaults, fn{i})
        warning('run_all_basins:unknownOption', ...
                'Ignoring unrecognized option "%s".', fn{i});
        continue;
    end
    defaults.(fn{i}) = extra.(fn{i});
end

if isempty(basins), basins = hc_basin_list(); end

% Drop skipped / unconfigured basins up front.
keep = true(1, numel(basins));
for b = 1:numel(basins)
    if (isfield(basins(b),'skip') && ~isempty(basins(b).skip) && basins(b).skip) ...
            || isempty(basins(b).stream_data_file)
        keep(b) = false;
        fprintf('Skipping %s (no stream data configured).\n', basins(b).name);
    end
end
basins = basins(keep);
nb = numel(basins);
if nb == 0
    error('run_all_basins:noBasins', ...
          'No runnable basins. Edit hc_basin_list.m and set stream_data_file.');
end

if ~exist(defaults.output_dir, 'dir'), mkdir(defaults.output_dir); end

%% ------------------------------------------- build one cfg per basin
% Done BEFORE the loop so each worker gets a fully-formed, self-contained
% struct (parfor cannot see variables built inside the loop body).
cfgs = cell(1, nb);
for b = 1:nb
    c = defaults;
    bs = basins(b);

    c.name     = bs.name;
    c.fileTag  = pick(bs, 'fileTag', matlab.lang.makeValidName(bs.name));
    c.stream_data_file = bs.stream_data_file;
    c.cave_data = pick(bs, 'cave_data', []);
    c.rng_seed  = pick(bs, 'rng_seed',  []);

    % Per-basin ksn prior, if measured for this basin.
    kr = pick(bs, 'ksn_relict', []);
    if ~isempty(kr)
        c.cave_prior.ksn_ref_mean = kr;
        sd = pick(bs, 'ksn_relict_sd', []);
        if ~isempty(sd)
            % Inflated as in the single-basin case, to absorb the
            % theta_ref mismatch between ksnTable and the sampled m/n.
            c.cave_prior.ksn_ref_std = max(sd * 1.8, 25);
        end
        % Start ksn_ref at the measured value for this basin.
        c.params_init(3) = kr;
    end

    c.params_init = pick(bs, 'params_init', c.params_init);
    c.pilot_file  = pick(bs, 'pilot_file',  c.pilot_file);
    c.verbose     = ~use_par;   % parfor output interleaves into noise
    cfgs{b} = c;
end

fprintf('\n=== Multi-basin inversion: %d basin(s) ===\n', nb);
fprintf('  %d burn-in + %d post-burn-in, dt = %g yr, proposal = %s\n', ...
        defaults.n_burnin, defaults.n_postburn, defaults.dt_forward, ...
        defaults.prop_mode);
fprintf('  Output: %s\n', defaults.output_dir);
if use_par, fprintf('  Mode: parallel (serial if no Parallel Toolbox)\n');
else,       fprintf('  Mode: serial\n'); end

%% ------------------------------------------------------------- run
res  = cell(1, nb);
errs = cell(1, nb);
odir = defaults.output_dir;

t_all = tic;
if use_par
    parfor b = 1:nb
        [res{b}, errs{b}] = runOne(cfgs{b}, odir);
    end
else
    for b = 1:nb
        [res{b}, errs{b}] = runOne(cfgs{b}, odir);
    end
end
fprintf('\nAll basins finished in %.0f s.\n', toc(t_all));

%% ---------------------------------------------------------- report
ok = ~cellfun(@isempty, res);
for b = find(~ok)
    fprintf('  FAILED  %-20s %s\n', cfgs{b}.name, errs{b});
end

summary = hc_multibasin_summary(res(ok), odir);

save(fullfile(odir, 'multibasin_summary.mat'), 'summary');
fprintf('\nSummary saved to %s\n', fullfile(odir, 'multibasin_summary.mat'));

end

%% ========================================================================
function [r, errmsg] = runOne(c, odir)
% One basin, isolated so a failure cannot take down the others.
r = []; errmsg = '';
try
    fprintf('--> starting %s\n', c.name);
    r = run_hc_inversion(c);
    % Save immediately: a crash later must not cost a completed basin.
    saveBasin(fullfile(odir, ['basin_' c.fileTag '.mat']), r);
    fprintf('<-- finished %s (accept %.1f%%, RMS %.1f m)\n', ...
            c.name, 100*r.diag.accept_post, r.diag.rms_map);
catch ME
    errmsg = ME.message;
    fprintf('!!! %s FAILED: %s\n', c.name, ME.message);
end
end

function saveBasin(fname, r) %#ok<INUSD>
% Separate function so `save` has a clean scope under parfor.
save(fname, '-struct', 'r');
end

function v = pick(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
