function results = run_hc_inversion(cfg)
% RUN_HC_INVERSION  Bayesian MCMC inversion for ONE basin, as a function.
%
%   results = run_hc_inversion(cfg)
%
% Function form of main_hc_bayes_inversion.m, so the same inversion can be
% driven over many basins from run_all_basins.  Deliberately contains NO
% plotting, NO `clear`, and no global state, which is what makes it safe to
% call inside parfor.
%
% cfg comes from hc_basin_defaults, with per-basin overrides applied by
% run_all_basins.  Required fields: stream_data_file, prior_bounds,
% params_init, p_steps, cave_prior, n_burnin, n_postburn, dt_forward,
% stream_err, N_eff_stream.
%
% Returns a struct with the chain, the MAP model, and a .diag substruct of
% health checks (acceptance, failures, railed parameters, step ratios).
%
% !! SYNC NOTE !!
% This is a port of main_hc_bayes_inversion.m and duplicates its MCMC loop.
% Until the two are reconciled (see MULTIBASIN_README.md), any change to
% the model or sampler in main_hc_bayes_inversion.m must be mirrored here,
% or multi-basin runs will silently use stale science.
%
% See also: hc_basin_defaults, hc_basin_list, run_all_basins

verbose = ~isfield(cfg,'verbose') || cfg.verbose;
vp = @(varargin) verbosePrint(verbose, varargin{:});

%% ---------------------------------------------------------------- setup
if isempty(cfg.rng_seed)
    rng('shuffle');
else
    rng(cfg.rng_seed);   % reproducible, and gives each basin its own stream
end

n_params   = numel(cfg.params_init);
total_iter = cfg.n_burnin + cfg.n_postburn;

%% ------------------------------------------------------------ load data
if isempty(cfg.stream_data_file) || ~exist(cfg.stream_data_file, 'file')
    error('run_hc_inversion:noData', ...
          'Stream data file not found: %s', cfg.stream_data_file);
end
raw = load(cfg.stream_data_file);
if isfield(raw, 'stream_data'), sd = raw.stream_data; else, sd = raw; end

S    = sd.S;
Sz   = double(sd.Sz(:));
S_DA = double(sd.S_DA(:));

if ~isa(S, 'STREAMobj')
    error('run_hc_inversion:notStreamObj', ...
         ['S loaded as class "%s". Add TopoToolbox to the path BEFORE ' ...
          'loading, or the class cannot be instantiated.'], class(S));
end

% Same fail-fast validation as the single-basin script: a NaN here would
% propagate silently into the likelihood.
if numel(Sz) ~= numel(S.IXgrid) || numel(S_DA) ~= numel(S.IXgrid)
    error('run_hc_inversion:sizeMismatch', ...
          'Sz (%d) / S_DA (%d) do not match S.IXgrid (%d).', ...
          numel(Sz), numel(S_DA), numel(S.IXgrid));
end
if any(~isfinite(Sz)),   error('Sz has %d non-finite values.',   sum(~isfinite(Sz)));   end
if any(~isfinite(S_DA)), error('S_DA has %d non-finite values.', sum(~isfinite(S_DA))); end
if any(S_DA <= 0),       error('S_DA has %d non-positive values.', sum(S_DA <= 0));     end

Sz_norm  = Sz - min(Sz);
n_stream = numel(Sz_norm);

%% ---------------------------------------------------------- cave data
if isfield(cfg,'cave_data') && ~isempty(cfg.cave_data)
    cave_ages       = cfg.cave_data(:,1);
    cave_heights    = cfg.cave_data(:,2);
    cave_height_err = cfg.cave_data(:,4);
    has_caves       = true;
else
    cave_ages = []; cave_heights = []; cave_height_err = [];
    has_caves = false;
end

% Likelihood weighting, kept IDENTICAL whether or not a basin has caves.
% hc_loglikelihood applies the N_eff down-weighting only when cave data is
% present; with no caves it would return the full-weight stream term, so a
% cave-less basin would carry ~50x more stream weight than a cave-bearing
% one and the two would not be comparable.  Folding the same factor into
% sigma instead reproduces the weighting exactly -- and inflating errors is
% how Gallen balances datasets in the first place.
if has_caves
    stream_sigma = cfg.stream_err * ones(n_stream, 1);
else
    inflate      = sqrt(n_stream / cfg.N_eff_stream);
    stream_sigma = cfg.stream_err * inflate * ones(n_stream, 1);
end

%% ------------------------------------------------------------- proposal
prop_L     = [];
prop_scale = 1.0;
if strcmpi(cfg.prop_mode, 'covariance')
    if isempty(cfg.pilot_file) || ~exist(cfg.pilot_file, 'file')
        error('run_hc_inversion:noPilot', ...
              'prop_mode=''covariance'' needs a pilot chain; got "%s".', ...
              cfg.pilot_file);
    end
    pilot = load(cfg.pilot_file, 'params_post');
    Sig   = cov(pilot.params_post);
    if size(Sig,1) == n_params
        Sigma = Sig;
    elseif size(Sig,1) == n_params - 1
        Sigma = zeros(n_params);
        Sigma(1:end-1,1:end-1) = Sig;
        Sigma(end,end)         = cfg.p_steps(end)^2;
    else
        error('Pilot has %d parameters; expected %d or %d.', ...
              size(Sig,1), n_params, n_params-1);
    end
    Sigma = Sigma + diag(max(diag(Sigma), realmin) * 1e-10);
    [prop_L, pd_fail] = chol((2.38^2/n_params) * Sigma, 'lower');
    if pd_fail ~= 0
        error('Pilot covariance not positive definite (minor %d).', pd_fail);
    end
end

%% --------------------------------------------------------------- init
params     = zeros(total_iter, n_params);
logL_chain = zeros(total_iter, 1);
accepted   = zeros(total_iter, 1);
params(1,:) = cfg.params_init;

p0     = cfg.params_init;
K_init = hc_derive_K(p0);
m_init = p0(5) * p0(4);
U_pat  = hc_uplift_pattern(S, p0(7));

Z_mod = hc_river_forward_model(S, S_DA, [p0(1) p0(2)], p0(6), ...
            K_init, m_init, p0(4), cfg.dt_forward, U_pat);

if has_caves
    cave_pred = cave_forward_model(cave_ages, [p0(1) p0(2)], p0(6));
else
    cave_pred = [];
end

logL_init = hc_loglikelihood(Sz_norm, Z_mod, stream_sigma, ...
                             cave_heights, cave_pred, cave_height_err);
if ~isfinite(logL_init)
    error('run_hc_inversion:badStart', ...
          'Initial log-likelihood is %s for basin %s.', ...
          num2str(logL_init), cfg.name);
end

logL_chain(1) = logL_init;
accepted(1)   = 1;
logP_map      = logL_init + logprior_hc(p0, cfg.prior_bounds, cfg.cave_prior);
params_map    = p0;
Z_mod_map     = Z_mod;
cave_pred_map = cave_pred;

vp('[%s] %d nodes, %d caves, initial logL = %.2f\n', ...
   cfg.name, n_stream, numel(cave_ages), logL_init);

%% ---------------------------------------------------------------- MCMC
t_start         = tic;
n_accept        = 0;
n_ffail         = 0;
n_nonfinite     = 0;
accept_window   = 0;
first_err_shown = false;
first_err_msg   = '';

for i = 2:total_iter

    if verbose && mod(i, 10000) == 0
        el = toc(t_start);
        vp('  [%s] %d/%d (%.1f/s, ~%.0f s left, accept=%.1f%%)\n', ...
           cfg.name, i, total_iter, i/el, (total_iter-i)/(i/el), ...
           100*n_accept/(i-1));
    end

    if cfg.adapt_steps && i <= cfg.n_burnin && i > 2 && ...
            mod(i-1, cfg.tune_interval) == 0
        win_acc    = accept_window / cfg.tune_interval;
        scale      = min(max(exp(win_acc - cfg.target_accept), 0.5), 2.0);
        prop_scale = prop_scale * scale;
        accept_window = 0;
    end

    current = params(i-1, :);

    if isempty(prop_L)
        candidate = current + prop_scale * (cfg.p_steps .* randn(1, n_params));
    else
        candidate = current + prop_scale * (prop_L * randn(n_params, 1))';
    end

    lp_cand = logprior_hc(candidate, cfg.prior_bounds, cfg.cave_prior);
    if lp_cand == -Inf
        params(i,:) = current;  logL_chain(i) = logL_chain(i-1);
        accepted(i) = 0;        continue;
    end
    lp_current = logprior_hc(current, cfg.prior_bounds, cfg.cave_prior);

    K_cand   = hc_derive_K(candidate);
    m_cand   = candidate(5) * candidate(4);
    U_pat_c  = hc_uplift_pattern(S, candidate(7));

    try
        Z_cand = hc_river_forward_model(S, S_DA, ...
                    [candidate(1) candidate(2)], candidate(6), ...
                    K_cand, m_cand, candidate(4), cfg.dt_forward, U_pat_c);
        if has_caves
            cave_cand = cave_forward_model(cave_ages, ...
                            [candidate(1) candidate(2)], candidate(6));
        else
            cave_cand = [];
        end
        logL_cand = hc_loglikelihood(Sz_norm, Z_cand, stream_sigma, ...
                                     cave_heights, cave_cand, cave_height_err);
    catch ME
        if ~first_err_shown
            first_err_shown = true;
            first_err_msg   = ME.message;
            vp('  [%s] first forward-model failure at iter %d: %s\n', ...
               cfg.name, i, ME.message);
        end
        params(i,:) = current;  logL_chain(i) = logL_chain(i-1);
        accepted(i) = 0;        n_ffail = n_ffail + 1;  continue;
    end

    % MATLAB's min() omits NaN, so min(NaN,0) is 0 and a NaN likelihood
    % would be accepted unconditionally, poisoning the rest of the chain.
    if ~isfinite(logL_cand)
        params(i,:) = current;  logL_chain(i) = logL_chain(i-1);
        accepted(i) = 0;        n_nonfinite = n_nonfinite + 1;  continue;
    end

    % Proposal is symmetric in both modes, so these cancel exactly.
    log_alpha = (lp_cand + logL_cand) - (lp_current + logL_chain(i-1));
    if isnan(log_alpha), log_alpha = -Inf; end

    if log(rand) < min(log_alpha, 0)
        params(i,:)   = candidate;
        logL_chain(i) = logL_cand;
        accepted(i)   = 1;
        n_accept      = n_accept + 1;
        accept_window = accept_window + 1;
        if (lp_cand + logL_cand) > logP_map
            logP_map      = lp_cand + logL_cand;
            params_map    = candidate;
            Z_mod_map     = Z_cand;
            cave_pred_map = cave_cand;
        end
    else
        params(i,:) = current;  logL_chain(i) = logL_chain(i-1);
        accepted(i) = 0;
    end
end
elapsed = toc(t_start);

%% ------------------------------------------------------------- package
params_post = params(cfg.n_burnin+1:end, :);
logL_post   = logL_chain(cfg.n_burnin+1:end);

results = struct();
results.name          = cfg.name;
results.cfg           = cfg;
results.params        = params;
results.params_post   = params_post;
results.params_map    = params_map;
results.logL_chain    = logL_chain;
results.logL_post     = logL_post;
results.accepted      = accepted;
results.Z_mod_map     = Z_mod_map;
results.cave_pred_map = cave_pred_map;
results.logP_map      = logP_map;
results.Sz_norm       = Sz_norm;

d = struct();
d.elapsed_s     = elapsed;
d.iter_per_s    = total_iter / max(elapsed, eps);
d.accept_all    = n_accept / total_iter;
d.accept_post   = mean(accepted(cfg.n_burnin+1:end));
d.n_ffail       = n_ffail;
d.n_nonfinite   = n_nonfinite;
d.first_err     = first_err_msg;
d.prop_scale    = prop_scale;
d.n_stream      = n_stream;
d.has_caves     = has_caves;
d.rms_map       = sqrt(mean((Sz_norm - Z_mod_map).^2));

% Step / posterior-sigma ratio (diagonal mode only): far below ~0.2 means
% the chain is crawling while acceptance looks deceptively healthy.
d.post_sd = std(params_post, 0, 1);
if isempty(prop_L)
    d.step_ratio = (prop_scale * cfg.p_steps) ./ max(d.post_sd, realmin);
else
    d.step_ratio = nan(1, n_params);
end

% Prior-bound railing: a posterior pressed against a wall is set by the
% prior, not the data, and drags its correlates with it.
d.railed = false(1, n_params);
d.ci95   = zeros(n_params, 2);
for j = 1:n_params
    lo = cfg.prior_bounds(j,1);  hi = cfg.prior_bounds(j,2);
    q  = prctile(params_post(:,j), [2.5, 97.5]);
    d.ci95(j,:) = q;
    d.railed(j) = (q(1)-lo) < 0.02*(hi-lo) || (hi-q(2)) < 0.02*(hi-lo);
end

d.K_post = hc_derive_K(params_post);
results.diag = d;

vp('[%s] done in %.0f s (%.1f it/s), accept %.1f%%, RMS %.1f m\n', ...
   cfg.name, elapsed, d.iter_per_s, 100*d.accept_post, d.rms_map);

end

%% ------------------------------------------------------------------------
function verbosePrint(flag, varargin)
if flag, fprintf(varargin{:}); end
end
