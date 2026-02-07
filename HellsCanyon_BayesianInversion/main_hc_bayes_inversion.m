% MAIN_HC_BAYES_INVERSION  Bayesian MCMC inversion of river profiles and
% cave burial ages for Hells Canyon drainage capture timing.
%
% Adapted from Gallen & Fernandez-Blanco (2021) bayes_profiler.
% Estimates stream power parameters (K, n, m/n) and capture timing
% simultaneously using Metropolis-Hastings MCMC.
%
% Two-phase tectonic model:
%   Phase 1: Pre-capture slow incision at rate U_pre
%   Phase 2: Post-capture rapid incision at rate U_post
%   Transition at time t_capture (years before present)
%
% Data constraints:
%   (1) River profile elevations from tributary DEMs
%   (2) Cave burial ages and heights above modern river
%
% Cave-derived priors from Morriss et al. (2025) PNAS:
%   t_capture ~ 2.1 +/- 1.0 Ma
%   U_pre     ~ 0.01 mm/yr (background rate from upper caves)
%   U_post    ~ 0.09-0.16 mm/yr (from lower caves)
%
% Parameters estimated [6 total]:
%   (1) U_pre      - Pre-capture incision rate (m/yr)
%   (2) U_post     - Post-capture incision rate (m/yr)
%   (3) log10(K)   - Log-erodibility
%   (4) n          - Slope exponent
%   (5) m/n        - Concavity ratio
%   (6) t_capture  - Capture timing (years)
%
% Required files on MATLAB path:
%   hc_river_forward_model.m
%   cave_forward_model.m
%   hc_loglikelihood.m
%   logprior_hc.m
%   logproposal_hc.m
%   prepare_hc_stream_data.m (to create input data)
%
% References:
%   Gallen & Fernandez-Blanco (2021), JGR-Earth Surface
%   Morriss et al. (2025), PNAS
%   Braun & Willett (2013), Geomorphology

clear; close all; clc;

%% ========================================================================
%  SECTION 1: USER CONFIGURATION
%  ========================================================================

% --- File paths ---
% Stream data: either pre-extracted .mat or DEM file
% If using prepare_hc_stream_data.m, point to the output .mat file
stream_data_file = '';  % <-- UPDATE: path to hc_stream_data.mat

% Output tag for saving results
fileTag = 'HC_capture';

% Output directory
output_dir = pwd;  % <-- UPDATE if desired

% --- Cave Data ---
% Cave burial ages (years before present) and heights above modern river (m)
% From Morriss et al. (2025) PNAS
% Format: [age_yr, height_m, age_uncertainty_yr, height_uncertainty_m]
%
% UPDATE THESE WITH YOUR ACTUAL CAVE DATA:
% Example values based on published constraints:
cave_data = [
    5.5e6,  100,  0.5e6,  20;   % Oldest cave (high, slow incision phase)
    3.5e6,  120,  0.5e6,  20;   % Middle cave
    1.5e6,  250,  0.3e6,  20;   % Youngest cave (rapid incision phase)
];
% Column 1: burial age (years)
% Column 2: height above modern river (m)
% Column 3: 1-sigma age uncertainty (years)
% Column 4: 1-sigma height uncertainty (m)

cave_ages       = cave_data(:,1);
cave_heights    = cave_data(:,2);
cave_age_err    = cave_data(:,3);
cave_height_err = cave_data(:,4);

% --- Whether to use cave data as PRIORS (informative) or just LIKELIHOOD ---
use_informative_priors = true;  % true = Gaussian priors from caves on timing/rates

% --- MCMC Settings ---
n_burnin   = 1e4;    % Burn-in iterations (increase for production: 3e5)
n_postburn = 1e5;    % Post-burn-in iterations (increase for production: 3e6)
dt_forward = 25000;  % Forward model time step (years)

% --- Run time for forward model ---
% Should be long enough to capture the full pre-capture steady state
run_time = 10e6;  % 10 Ma (generous to allow pre-capture equilibrium)

%% ========================================================================
%  SECTION 2: PARAMETER SETUP
%  ========================================================================

% Parameter ordering:
% [1] U_pre     - Pre-capture incision rate (m/yr)
% [2] U_post    - Post-capture incision rate (m/yr)
% [3] log10(K)  - Log-erodibility
% [4] n         - Slope exponent
% [5] m/n       - Concavity ratio (theta)
% [6] t_capture - Capture timing (years before present)

% --- Prior bounds (uniform) ---
prior_bounds = [
    1e-6,   5e-4;     % U_pre:     ~0.001 to 0.5 mm/yr
    1e-5,   5e-3;     % U_post:    ~0.01 to 5 mm/yr
    -9,     -4;       % log10(K):  wide range
    0.5,    10;       % n:         0.5 to 10 (allow high n like Gallen found)
    0.3,    0.7;      % m/n:       typical concavity range
    0.5e6,  5e6;      % t_capture: 0.5 to 5 Ma
];

% --- Informative priors from cave constraints ---
cave_prior = struct();
cave_prior.use_informative = use_informative_priors;
cave_prior.t_capture_mean  = 2.1e6;    % 2.1 Ma from PNAS paper
cave_prior.t_capture_std   = 1.0e6;    % +/- 1.0 Ma
cave_prior.U_pre_mean      = 1e-5;     % 0.01 mm/yr
cave_prior.U_pre_std       = 5e-6;     % +/- 0.005 mm/yr
cave_prior.U_post_mean     = 1.25e-4;  % 0.125 mm/yr (midpoint of 0.09-0.16)
cave_prior.U_post_std      = 3.5e-5;   % +/- 0.035 mm/yr

% --- Starting values (near expected MAP for faster convergence) ---
params_init = [
    1e-5,     ... % U_pre = 0.01 mm/yr
    1.25e-4,  ... % U_post = 0.125 mm/yr
    -6.24,    ... % log10(K) ~ 5.7e-7
    1.0,      ... % n = 1 (start at linear)
    0.45,     ... % m/n = 0.45
    2.1e6     ... % t_capture = 2.1 Ma
];

% --- MCMC step sizes (tune for ~25-50% acceptance) ---
p_steps = [
    2e-6,    ... % U_pre
    1e-5,    ... % U_post
    0.05,    ... % log10(K)
    0.05,    ... % n
    0.01,    ... % m/n
    5e4      ... % t_capture (50 kyr steps)
];

n_params = length(params_init);

%% ========================================================================
%  SECTION 3: LOAD STREAM DATA
%  ========================================================================

if ~isempty(stream_data_file) && exist(stream_data_file, 'file')
    fprintf('Loading stream data from: %s\n', stream_data_file);
    tmp = load(stream_data_file);
    if isfield(tmp, 'stream_data')
        sd = tmp.stream_data;
    else
        sd = tmp;
    end

    S    = sd.S;
    Sz   = sd.Sz;
    S_DA = sd.S_DA;
else
    fprintf('WARNING: No stream data file specified.\n');
    fprintf('Please run prepare_hc_stream_data.m first to create hc_stream_data.mat\n');
    fprintf('Then set stream_data_file at the top of this script.\n');
    fprintf('\nExample:\n');
    fprintf('  sd = prepare_hc_stream_data(''path/to/Basin_80_Data.mat'');\n');
    fprintf('  % Then set stream_data_file = ''hc_stream_data.mat'';\n');
    return;
end

% Normalize elevations relative to outlet
Sz_norm = Sz - min(Sz);

% Stream data error (meters)
stream_err_base = 5;  % 5 m base uncertainty
n_stream = length(Sz_norm);
n_cave   = length(cave_ages);

% Weight adjustment for balanced likelihood
Ws = sqrt(n_cave / n_stream);
stream_err = stream_err_base / Ws;  % inflate stream error to balance

fprintf('Data loaded: %d stream nodes, %d cave observations\n', n_stream, n_cave);

%% ========================================================================
%  SECTION 4: INITIALIZE MCMC
%  ========================================================================

total_iter = n_burnin + n_postburn;

% Storage arrays
params      = zeros(total_iter, n_params);
logL_chain  = zeros(total_iter, 1);
accepted    = zeros(total_iter, 1);

% Initialize
params(1,:) = params_init;

% Run initial forward model
K_init = 10^params_init(3);
m_init = params_init(5) * params_init(4);  % m = (m/n) * n

fprintf('Running initial forward model...\n');
Z_mod = hc_river_forward_model(S, S_DA, params_init(1), params_init(2), ...
    params_init(6), K_init, m_init, params_init(4), run_time, dt_forward);

% Initial cave predictions
cave_pred = cave_forward_model(cave_ages, params_init(6), ...
    params_init(1), params_init(2));

% Initial likelihood
[logL_init, ~, ~] = hc_loglikelihood(Sz_norm, Z_mod, ...
    stream_err * ones(n_stream, 1), cave_heights, cave_pred, cave_height_err);

logL_chain(1) = logL_init;
accepted(1)   = 1;

% Track MAP (maximum a posteriori)
logP_map = logL_init + logprior_hc(params_init, prior_bounds, cave_prior);
params_map = params_init;
Z_mod_map  = Z_mod;
cave_pred_map = cave_pred;

fprintf('Initial log-likelihood: %.2f\n', logL_init);
fprintf('\nStarting MCMC: %d burn-in + %d post-burn-in iterations\n', ...
    n_burnin, n_postburn);

%% ========================================================================
%  SECTION 5: RUN MCMC
%  ========================================================================

tic;
n_accept = 0;

for i = 2:total_iter
    % Progress report
    if mod(i, 1000) == 0
        elapsed = toc;
        rate = i / elapsed;
        remaining = (total_iter - i) / rate;
        fprintf('  Iter %d/%d (%.0f/s, ~%.0f s remaining, accept=%.1f%%)\n', ...
            i, total_iter, rate, remaining, 100*n_accept/(i-1));
    end

    % Current parameters
    current = params(i-1, :);

    % Propose new parameters (Gaussian random walk)
    candidate = current + p_steps .* randn(1, n_params);

    % Check prior
    lp_cand = logprior_hc(candidate, prior_bounds, cave_prior);

    if lp_cand == -Inf
        % Rejected by prior - skip forward model (expensive)
        params(i,:) = current;
        logL_chain(i) = logL_chain(i-1);
        accepted(i) = 0;
        continue;
    end

    % Compute prior for current
    lp_current = logprior_hc(current, prior_bounds, cave_prior);

    % Run forward model with candidate parameters
    K_cand = 10^candidate(3);
    m_cand = candidate(5) * candidate(4);  % m = (m/n) * n

    try
        Z_cand = hc_river_forward_model(S, S_DA, candidate(1), candidate(2), ...
            candidate(6), K_cand, m_cand, candidate(4), run_time, dt_forward);

        cave_cand = cave_forward_model(cave_ages, candidate(6), ...
            candidate(1), candidate(2));

        % Likelihood
        [logL_cand, ~, ~] = hc_loglikelihood(Sz_norm, Z_cand, ...
            stream_err * ones(n_stream, 1), cave_heights, cave_cand, cave_height_err);
    catch
        % Forward model failed - reject
        params(i,:) = current;
        logL_chain(i) = logL_chain(i-1);
        accepted(i) = 0;
        continue;
    end

    % Proposal probabilities (symmetric, so they cancel)
    lq_fwd = logproposal_hc(current, candidate, p_steps);
    lq_rev = logproposal_hc(candidate, current, p_steps);

    % Log acceptance ratio
    log_alpha = (lp_cand + logL_cand + lq_rev) - ...
                (lp_current + logL_chain(i-1) + lq_fwd);

    % Metropolis-Hastings accept/reject
    if log(rand) < min(log_alpha, 0)
        % Accept
        params(i,:) = candidate;
        logL_chain(i) = logL_cand;
        accepted(i) = 1;
        n_accept = n_accept + 1;

        % Update MAP
        logP_cand = lp_cand + logL_cand;
        if logP_cand > logP_map
            logP_map = logP_cand;
            params_map = candidate;
            Z_mod_map = Z_cand;
            cave_pred_map = cave_cand;
        end
    else
        % Reject
        params(i,:) = current;
        logL_chain(i) = logL_chain(i-1);
        accepted(i) = 0;
    end
end

elapsed_total = toc;
fprintf('\nMCMC complete: %.0f seconds (%.1f iter/s)\n', elapsed_total, ...
    total_iter / elapsed_total);
fprintf('Overall acceptance rate: %.1f%%\n', 100 * n_accept / total_iter);

%% ========================================================================
%  SECTION 6: EXTRACT POST-BURN-IN RESULTS
%  ========================================================================

params_post = params(n_burnin+1:end, :);
logL_post   = logL_chain(n_burnin+1:end);

% Compute statistics
param_names = {'U_{pre} (m/yr)', 'U_{post} (m/yr)', 'log_{10}(K)', ...
               'n', 'm/n', 't_{capture} (yr)'};
param_units = {'mm/yr', 'mm/yr', '', '', '', 'Ma'};
param_scale = [1e3, 1e3, 1, 1, 1, 1e-6];  % for display

fprintf('\n========== POSTERIOR SUMMARY ==========\n');
fprintf('%-20s %12s %12s %20s %20s\n', 'Parameter', 'MAP', 'Median', '68% CI', '95% CI');
fprintf('%s\n', repmat('-', 1, 84));

for j = 1:n_params
    map_val  = params_map(j) * param_scale(j);
    med_val  = median(params_post(:,j)) * param_scale(j);
    ci68     = prctile(params_post(:,j), [16, 84]) * param_scale(j);
    ci95     = prctile(params_post(:,j), [2.5, 97.5]) * param_scale(j);

    fprintf('%-20s %12.4g %12.4g [%8.4g, %8.4g] [%8.4g, %8.4g]\n', ...
        param_names{j}, map_val, med_val, ci68(1), ci68(2), ci95(1), ci95(2));
end

% K in real space
K_post = 10.^params_post(:,3);
fprintf('\nK (real space): median = %.2e, 68%% CI = [%.2e, %.2e]\n', ...
    median(K_post), prctile(K_post, 16), prctile(K_post, 84));

%% ========================================================================
%  SECTION 7: SAVE RESULTS
%  ========================================================================

save(fullfile(output_dir, ['params_' fileTag '.mat']), 'params', 'params_post', ...
    'params_map', 'logL_chain', 'logL_post', 'accepted', 'prior_bounds', ...
    'cave_prior', 'p_steps', 'n_burnin', 'n_postburn', 'param_names');

save(fullfile(output_dir, ['mMAP_' fileTag '.mat']), 'params_map', ...
    'Z_mod_map', 'cave_pred_map', 'logP_map');

fprintf('\nResults saved to: %s\n', output_dir);

%% ========================================================================
%  SECTION 8: BASIC DIAGNOSTIC PLOTS
%  ========================================================================

plot_hc_results(params, logL_chain, n_burnin, params_map, Z_mod_map, ...
    Sz_norm, S, cave_ages, cave_heights, cave_height_err, cave_pred_map, ...
    prior_bounds, cave_prior, param_names, param_scale, output_dir, fileTag);

fprintf('\nDone. Run plot_hc_results.m for additional visualizations.\n');
