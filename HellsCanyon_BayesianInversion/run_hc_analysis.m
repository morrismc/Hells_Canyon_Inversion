% RUN_HC_ANALYSIS  Quick-start script for Hells Canyon inversion analyses.
%
% This script demonstrates how to run the three main analyses:
%   (1) Ksn vs erosion rate analysis (estimate K and n)
%   (2) Linear inversion with cave constraints (n=1, multiple K values)
%   (3) Full Bayesian MCMC inversion (free n, cave priors)
%
% UPDATE THE PATHS BELOW before running.
%
% Author: Adapted for Hells Canyon from Gallen (2018, 2021) framework

clear; close all; clc;

%% ========================================================================
%  PATHS - UPDATE THESE
%  ========================================================================

% Path to ksnTable.xlsx
ksn_file = '';  % <-- UPDATE: e.g., 'C:\Users\...\ksnTable.xlsx'

% Path to DEM (.mat or .tif)
dem_file = '';  % <-- UPDATE: e.g., 'C:\Users\...\Basin_80_Data.mat'

% Output directory
output_dir = fullfile(pwd, 'results');
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

% Add code directories to path
addpath(genpath(fullfile(fileparts(mfilename('fullpath')))));
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', ...
    'Block_Uplift_Linear_Inversion_Models')));

%% ========================================================================
%  CAVE CONSTRAINTS (from Morriss et al. 2025, PNAS)
%  ========================================================================

% Cave burial ages and heights above modern Snake River
% UPDATE THESE with your actual cave data:
cave_data = struct();
cave_data.ages      = [5.5e6, 3.5e6, 1.5e6];  % years before present
cave_data.heights   = [100,   120,   250];      % meters above river
cave_data.age_err   = [0.5e6, 0.5e6, 0.3e6];   % 1-sigma age uncertainty
cave_data.height_err= [20,    20,    20];        % 1-sigma height uncertainty

% Derived incision rates from caves
cave_data.U_pre     = 0.01;    % mm/yr (background, from upper caves)
cave_data.U_pre_err = 0.005;   % mm/yr
cave_data.U_post    = 0.125;   % mm/yr (post-capture, from lower caves; range 0.09-0.16)
cave_data.U_post_err= 0.035;   % mm/yr
cave_data.t_capture = 2.1e6;   % years (drainage capture timing)
cave_data.t_capture_err = 1.0e6; % years

fprintf('============================================\n');
fprintf('  HELLS CANYON INVERSION ANALYSIS\n');
fprintf('============================================\n\n');

%% ========================================================================
%  ANALYSIS 1: Ksn vs Erosion Rate (Estimate K and n)
%  ========================================================================

fprintf('\n--- ANALYSIS 1: Ksn vs Erosion Rate ---\n\n');

if ~isempty(ksn_file) && exist(ksn_file, 'file')
    ksn_results = ksn_erosion_analysis(ksn_file, ...
        'E_relict',      cave_data.U_pre, ...
        'E_relict_err',  cave_data.U_pre_err, ...
        'E_adjusted',    cave_data.U_post, ...
        'E_adjusted_err',cave_data.U_post_err, ...
        'n_samples',     100000, ...
        'output_dir',    output_dir);

    fprintf('\n  Best K (Bayesian): %.2e (n=%.1f)\n', ...
        ksn_results.K_median, ksn_results.n_median);
    fprintf('  Best K (n=1, adjusted): %.2e\n', ksn_results.K_n1_adjusted);
    fprintf('  Best K (n=1, relict):   %.2e\n', ksn_results.K_n1_relict);

    % Store K values for downstream analyses
    K_low  = ksn_results.K_ci95(1);
    K_mid  = ksn_results.K_n1_adjusted;  % Use n=1 adjusted for linear inversion
    K_high = ksn_results.K_ci95(2);
else
    fprintf('  SKIPPED: No ksn_file specified. Using default K values.\n');
    % Default K values from previous analysis
    K_low  = 7.0e-8;
    K_mid  = 5.73e-7;
    K_high = 7.0e-7;
    ksn_results = [];
end

fprintf('\n  K range for inversions: [%.1e, %.1e, %.1e]\n', K_low, K_mid, K_high);

%% ========================================================================
%  ANALYSIS 2: Linear Inversion with Cave Overlay (n=1)
%  ========================================================================

fprintf('\n--- ANALYSIS 2: Linear Inversion (n=1, Gallen method) ---\n\n');

if ~isempty(dem_file) && exist(dem_file, 'file')
    % Load DEM
    [~, ~, ext] = fileparts(dem_file);
    if strcmpi(ext, '.mat')
        tmp = load(dem_file);
        % Try common variable names
        if isfield(tmp, 'DEMoc')
            DEM = tmp.DEMoc;
        elseif isfield(tmp, 'tDEM')
            DEM = tmp.tDEM;
        elseif isfield(tmp, 'DEM')
            DEM = tmp.DEM;
        else
            fn = fieldnames(tmp);
            for fi = 1:length(fn)
                if isa(tmp.(fn{fi}), 'GRIDobj')
                    DEM = tmp.(fn{fi});
                    break;
                end
            end
        end
    else
        DEM = GRIDobj(dem_file);
    end

    % Run linear inversions at multiple K values
    K_test = [K_low, K_mid, K_high];

    lin_results = hc_linear_inversion_with_caves(DEM, K_test, ...
        'tau_inc',       1.5e6, ...
        'crita',         1e6, ...
        'mn',            0.45, ...
        'flowOption',    'fill', ...
        'cave_ages',     cave_data.ages, ...
        'cave_heights',  cave_data.heights, ...
        'cave_err',      cave_data.height_err, ...
        't_capture',     cave_data.t_capture, ...
        't_capture_err', cave_data.t_capture_err, ...
        'output_dir',    output_dir);

    % Report which K values give capture signal consistent with caves
    fprintf('\n  Linear inversion results:\n');
    for k = 1:length(K_test)
        if isfield(lin_results(k), 'capture_in_window')
            fprintf('    K = %.1e: peak at %.1f Ma (%s cave window)\n', ...
                lin_results(k).K, lin_results(k).t_peak/1e6, ...
                ternary(lin_results(k).capture_in_window, 'IN', 'OUTSIDE'));
        end
    end
else
    fprintf('  SKIPPED: No dem_file specified.\n');
    fprintf('  To run this analysis, set dem_file to your Basin DEM path.\n');
    lin_results = [];
end

%% ========================================================================
%  ANALYSIS 3: Bayesian MCMC Inversion (free n)
%  ========================================================================

fprintf('\n--- ANALYSIS 3: Bayesian MCMC Inversion ---\n\n');

if ~isempty(dem_file) && exist(dem_file, 'file')
    fprintf('  To run the full MCMC, edit and run main_hc_bayes_inversion.m\n');
    fprintf('  Key settings to update:\n');
    fprintf('    stream_data_file = ''path/to/hc_stream_data.mat''\n');
    fprintf('    cave_data array with your actual cave ages/heights\n');
    fprintf('    n_burnin / n_postburn for production runs\n');
    fprintf('\n  First prepare stream data:\n');
    fprintf('    sd = prepare_hc_stream_data(''%s'');\n', dem_file);

    % Optionally prepare stream data now
    % Uncomment to run:
    % sd = prepare_hc_stream_data(dem_file, 'output_file', ...
    %     fullfile(output_dir, 'hc_stream_data.mat'));
else
    fprintf('  SKIPPED: No dem_file specified.\n');
end

%% ========================================================================
%  SUMMARY
%  ========================================================================

fprintf('\n============================================\n');
fprintf('  ANALYSIS SUMMARY\n');
fprintf('============================================\n\n');

fprintf('Cave constraints (Morriss et al. 2025):\n');
fprintf('  Capture timing: %.1f +/- %.1f Ma\n', ...
    cave_data.t_capture/1e6, cave_data.t_capture_err/1e6);
fprintf('  Pre-capture rate:  %.3f +/- %.3f mm/yr\n', ...
    cave_data.U_pre, cave_data.U_pre_err);
fprintf('  Post-capture rate: %.3f +/- %.3f mm/yr\n', ...
    cave_data.U_post, cave_data.U_post_err);

if ~isempty(ksn_results)
    fprintf('\nErodibility from ksn analysis:\n');
    fprintf('  K (n=1, adjusted):  %.2e\n', ksn_results.K_n1_adjusted);
    fprintf('  K (Bayesian, n=%.1f): %.2e\n', ksn_results.n_median, ksn_results.K_median);
    fprintf('  Evidence for n!=1: implied n = %.1f (%.1f to %.1f)\n', ...
        ksn_results.n_uniform_K, ksn_results.n_uniform_K_range);
end

fprintf('\nOutput directory: %s\n', output_dir);
fprintf('\nRecommended next steps:\n');
fprintf('  1. Update cave_data with your actual PNAS cave ages/heights\n');
fprintf('  2. Review ksn_erosion_analysis results for K and n\n');
fprintf('  3. Run linear inversion to check timing consistency\n');
fprintf('  4. Run full MCMC for simultaneous K, n, timing estimation\n');
fprintf('  5. Compare Bayesian n posterior with n=1 assumption\n');

function s = ternary(cond, a, b)
if cond; s = a; else; s = b; end
end
