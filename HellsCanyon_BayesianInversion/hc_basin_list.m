function basins = hc_basin_list()
% HC_BASIN_LIST  Per-basin definitions for a multi-basin inversion.
%
% Each entry overrides only what genuinely differs between basins; the
% shared settings (burn-in, timestep, bounds, proposal tuning, cave priors)
% come from hc_basin_defaults and are deliberately NOT repeated here.
%
% TEMPLATE -- the paths below are placeholders.  Fill in one basin, run it,
% confirm it reproduces your single-basin result, and only then add the
% rest.  See MULTIBASIN_README.md.
%
% Fields (all optional except name / stream_data_file):
%   name            Human-readable basin name
%   stream_data_file  hc_stream_data-style .mat from prepare_hc_stream_data
%   fileTag         Prefix for output files (defaults to a cleaned name)
%   cave_data       [N x 4] age, height, age_err, height_err -- or [] for
%                   basins with no caves (the common case for tributaries)
%   ksn_relict      Mean ksn of above-knickpoint segments, for this basin
%   ksn_relict_sd   Its standard deviation
%   ksn_adjusted    Mean ksn of below-knickpoint segments
%   params_init     Per-basin starting values (else the global default)
%   rng_seed        Integer for reproducibility (else shuffle)
%   skip            true to leave a basin out without deleting its entry
%
% NOTE: ksnTable.xlsx carries a stream_ID column, so per-basin ksn values
% can be pulled from it directly rather than typed in by hand -- see
% ksn_erosion_analysis.m.
%
% See also: hc_basin_defaults, run_all_basins

root = 'C:\Users\mmorriss\Desktop\Side_projects\Hells_Canyon_Inversion';

basins = struct([]);
k = 0;

%% --- Basin 1: the trunk basin already validated -------------------------
k = k + 1;
basins(k).name             = 'Granite_Creek';
basins(k).stream_data_file = fullfile(root, 'HellsCanyon_BayesianInversion', ...
                                      'hc_stream_data.mat');
basins(k).fileTag          = 'granite';
basins(k).cave_data = [
    1.5e6,   242,  0.73e6,  20;
    2.71e6,  343,  1.4e6,   20;
    5.5e6,   375,  3.8e6,   20;
];
basins(k).ksn_relict    = 113.3;
basins(k).ksn_relict_sd = 33.8;
basins(k).ksn_adjusted  = 227.4;
basins(k).rng_seed      = 101;
basins(k).skip          = false;

%% --- Basins 2-5: tributaries (PLACEHOLDERS -- fill in) ------------------
% Tributaries will typically have cave_data = [], so their incision rates
% are carried entirely by the shared cave priors in hc_basin_defaults.
trib_names = {'Sheep_Creek', 'Basin_03', 'Basin_04', 'Basin_05'};
trib_files = {'', '', '', ''};   % <-- point these at real stream data

for t = 1:numel(trib_names)
    k = k + 1;
    basins(k).name             = trib_names{t};
    basins(k).stream_data_file = trib_files{t};
    basins(k).fileTag          = lower(trib_names{t});
    basins(k).cave_data        = [];      % no caves in this basin
    basins(k).ksn_relict       = [];      % fill from ksnTable per stream_ID
    basins(k).ksn_relict_sd    = [];
    basins(k).ksn_adjusted     = [];
    basins(k).rng_seed         = 100 + k;
    % Skipped until a real stream_data_file is supplied, so the driver can
    % be run end-to-end on basin 1 alone without editing this file.
    basins(k).skip             = isempty(trib_files{t});
end

end
