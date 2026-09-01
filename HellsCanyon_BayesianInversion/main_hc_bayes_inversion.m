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
% The two-phase model captures the primary capture-driven incision
% signal.  Residual structure in the trunk profile (additional knickpoints
% higher in the catchment) likely reflects pre-capture changes in uplift
% or incision rate that are not the focus of this analysis but are
% highlighted in the trunk residual plot.
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
%   (3) ksn_ref    - Relict-channel steepness index
%   (4) n          - Slope exponent
%   (5) m/n        - Concavity ratio
%   (6) t_capture  - Capture timing (years)
%
% PARAMETERIZATION NOTE -- K is DERIVED, not sampled:
%
%       K = U_pre / ksn_ref^n
%
% Earlier versions sampled log10(K) directly.  That is badly conditioned
% here: at steady state the relict channel has ksn = (U_pre/K)^(1/n), so
% the K that fits the data moves by ~6 orders of magnitude as n varies.
% Sampling log10(K) and n as independent coordinates therefore proposes
% ACROSS a very narrow curved ridge, which produced integrated
% autocorrelation times of ~5400 (ESS ~9 from 50,000 samples), and the
% [-9,-4] box on log10(K) excluded the solution the ksn data implies
% (log10 K ~ -12.5 at n ~ 3.7) altogether.
%
% Sampling ksn_ref instead removes both problems: it is directly
% observable (measured 113.3 +/- 33.8 from 38 relict channel segments),
% it is only weakly correlated with n, and K can no longer be boxed out
% of its physically required range because it is computed rather than
% bounded.
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

addpath(genpath("C:\Users\mmorriss\Documents\matlab\topotoolbox-master"))
addpath(genpath("C:\Users\mmorriss\Documents\matlab\Topographic-Analysis-Kit-master"))
addpath(genpath("C:\Users\mmorriss\Documents\matlab\TT_shortcourse_6august2016"))
addpath(genpath("C:\Users\mmorriss\Documents\matlab\knickdisttest_morriss"))
%% ========================================================================
%  SECTION 1: USER CONFIGURATION
%  ========================================================================

% --- File paths ---
stream_data_file = 'C:\Users\mmorriss\Desktop\Side_projects\Hells_Canyon_Inversion\HellsCanyon_BayesianInversion\hc_stream_data.mat';

% Output tag for saving results
fileTag = 'HC_capture';

% Output directory
output_dir = pwd;

% --- Cave Data ---
% Cave burial ages (years before present) and heights above modern river (m)
% From Morriss et al. (2025) PNAS
cave_data = [
    1.5e6,   242,  0.73e6,  20;   % Youngest cave (rapid incision phase)
    2.71e6,  343,  1.4e6,   20;   % Middle cave
    5.5e6,   375,  3.8e6,   20;   % Oldest cave (high, slow incision phase)
];

cave_ages       = cave_data(:,1);
cave_heights    = cave_data(:,2);
cave_age_err    = cave_data(:,3);  % NOT used in the likelihood: ages are
                                   % treated as exact (as Gallen treats
                                   % terrace ages); age uncertainty enters
                                   % via the informative t_capture prior.
cave_height_err = cave_data(:,4);

% --- Whether to use cave data as PRIORS (informative) ---
use_informative_priors = true;

% --- MCMC Settings ---
% Production settings.  Sized from the measured IACT (~450 for the worst
% parameter, m/n and ksn_ref) with headroom: see the truncation caveat in
% ess_check.  n_burnin raised to 2e4 because n needed several thousand
% iterations to travel from its 3.7 start to ~6 and settle.
n_burnin   = 2e4;    % Burn-in iterations
n_postburn = 2e5;    % Post-burn-in iterations

% dt = 100,000 validated against a 6,250 yr reference AT n ~ 6 (the
% shock-forming regime), giving 0.41 m RMS / 3.90 m max deviation against
% a stream_err of 200 m -- a ~500x margin -- for a 2.7x cheaper forward
% model.  Do NOT raise this further without re-running dt_sensitivity.
dt_forward = 100000; % Forward model time step (years)

% --- Adaptive proposal tuning (burn-in ONLY) ---
% NOTE: Gallen & Fernandez-Blanco (2021) hand-tune FIXED step sizes by
% trial and error ("a tedious process") until the acceptance rate falls in
% the ~20-60% range.  Here that manual tuning is automated: step sizes are
% rescaled during burn-in to drive acceptance toward the target, then
% FROZEN once burn-in ends so the post-burn-in chain is a proper
% stationary Markov chain (adapting during sampling would break detailed
% balance).  The post-burn-in sampler is therefore identical in form to
% Gallen's fixed-step Metropolis-Hastings.
adapt_steps   = true;   % enable burn-in step-size tuning
tune_interval = 200;    % iterations between step-size updates
target_accept = 0.30;   % target acceptance rate (optimal ~0.23-0.44)

%% ========================================================================
%  SECTION 2: PARAMETER SETUP
%  ========================================================================

% Parameter ordering [6 total]:
% [1] U_pre      - Pre-capture incision rate (m/yr)
% [2] U_post     - Post-capture incision rate (m/yr)
% [3] ksn_ref    - Relict channel steepness  (K = U_pre / ksn_ref^n)
% [4] n          - Slope exponent
% [5] m/n        - Concavity ratio (theta)
% [6] t_capture  - Capture timing (years BP)

% --- Prior bounds (uniform hard walls) ---
prior_bounds = [
    1e-6,   5e-4;     % U_pre:     ~0.001 to 0.5 mm/yr
    1e-5,   5e-3;     % U_post:    ~0.01 to 5 mm/yr
    20,     500;      % ksn_ref:   relict steepness (measured ~113)
    0.5,    12;       % n:         raised from 8; posterior sits at ~6 and
                      %            must not be able to press the ceiling
                      %            (Gallen's bayes_profiler allowed n to 20)
    0.3,    0.8;      % m/n:       typical concavity range
    0.5e6,  5e6;      % t_capture: 0.5 to 5 Ma
];

% --- Informative priors from cave constraints ---
cave_prior = struct();
cave_prior.use_informative = use_informative_priors;

% Capture timing from cave burial dating
cave_prior.t_capture_mean  = 2.1e6;    % 2.1 Ma from PNAS paper
cave_prior.t_capture_std   = 1.0e6;    % +/- 1.0 Ma

% Background incision rate from pre-capture caves
cave_prior.U_pre_mean      = 1e-5;     % 0.01 mm/yr
cave_prior.U_pre_std       = 5e-6;     % +/- 0.005 mm/yr

% Post-capture incision rate from lower caves
cave_prior.U_post_mean     = 1.25e-4;  % 0.125 mm/yr (midpoint of 0.09-0.16)
cave_prior.U_post_std      = 3.5e-5;   % +/- 0.035 mm/yr

% Relict channel steepness, measured from 38 above-knickpoint segments
% in ksnTable.xlsx (mean 113.3, std 33.8; see ksn_erosion_analysis.m).
%
% The std is deliberately inflated to 60 rather than the measured 33.8.
% ksn depends on the reference concavity used to compute it, and the
% ksnTable was built with TAK's theta_ref (0.45 by default) whereas the
% model's own normalization is m/n, which the chain samples freely.  The
% widened prior absorbs that normalization mismatch and keeps ksn_ref
% acting mainly as a well-scaled sampling coordinate rather than a hard
% constraint.
%
% ACTION: confirm the theta_ref actually used to build ksnTable.xlsx.  If
% it differs materially from the posterior m/n, widen this further or
% disable it entirely by setting ksn_ref_std = 0.
cave_prior.ksn_ref_mean    = 113.3;
cave_prior.ksn_ref_std     = 60.0;

% NO informative prior on n.  Previous runs used n ~ N(1, 0.5), which was
% fighting the data: the measured ksn ratio (227.4/113.3 = 2.01) together
% with the cave rate ratio (U_post/U_pre = 13.2) implies
%       n = ln(13.2)/ln(2.01) = 3.7  [3.3 - 4.2]
% i.e. 5.4 sigma from that prior's mean, while the old n bound of 3
% excluded it outright.  Leaving n uninformative lets the river profile
% determine it independently, so the ksn-derived 3.7 remains a genuine
% external check rather than an assumption folded into the prior.
% (Set n_std > 0 below to re-enable a Gaussian prior on n.)
cave_prior.n_mean          = 3.7;
cave_prior.n_std           = 0;     % 0 = disabled (uniform within bounds)

% --- Starting values ---
% Started near the approximate solution, as Gallen does in bayes_profiler.
% ksn_ref and n come from the independent ksn analysis; the rates and
% timing from the previous run's posterior medians.
params_init = [
    1.17e-5,  ... % U_pre   = 0.0117 mm/yr (posterior median)
    1.55e-4,  ... % U_post  = 0.155 mm/yr  (posterior median)
    113.3,    ... % ksn_ref = measured relict steepness
    3.7,      ... % n       = implied by ksn ratio + cave rate ratio
    0.6,      ... % m/n     = previous posterior median
    2.19e6    ... % t_capture = 2.19 Ma (posterior median)
];

% --- MCMC step sizes ---
% Only the RELATIVE proportions matter: the burn-in tuner rescales the
% whole vector by one shared factor to hit the target acceptance rate.
p_steps = [
    2e-6,    ... % U_pre   (tightly constrained by cave prior)
    1.5e-5,  ... % U_post
    5.0,     ... % ksn_ref (~0.08 of its prior std)
    0.10,    ... % n       (wider bound now: [0.5, 8])
    0.006,   ... % m/n
    7.5e4    ... % t_capture (75 kyr steps)
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
    Sz   = sd.Sz(:);
    S_DA = sd.S_DA(:);

    % Validate that node counts are consistent
    n_nodes = numel(S.IXgrid);
    if length(Sz) ~= n_nodes
        warning('Sz (%d) and S.IXgrid (%d) have different lengths. Re-extracting from S.', ...
            length(Sz), n_nodes);
        if isfield(sd, 'DEM') && isa(sd.DEM, 'GRIDobj')
            Sz   = double(sd.DEM.Z(S.IXgrid));
        else
            Sz = Sz(1:min(length(Sz), n_nodes));
            if length(Sz) < n_nodes
                error('Cannot reconcile Sz length with S: %d vs %d nodes', length(Sz), n_nodes);
            end
        end
    end
    if length(S_DA) ~= n_nodes
        warning('S_DA (%d) and S.IXgrid (%d) have different lengths.', ...
            length(S_DA), n_nodes);
        S_DA = S_DA(1:min(length(S_DA), n_nodes));
        if length(S_DA) < n_nodes
            error('Cannot reconcile S_DA length with S: %d vs %d nodes', length(S_DA), n_nodes);
        end
    end
else
    fprintf('WARNING: No stream data file specified.\n');
    fprintf('Please run prepare_hc_stream_data.m first to create hc_stream_data.mat\n');
    return;
end

% --- Validate and condition the input data -----------------------------
% TopoToolbox GRIDobj data is normally SINGLE precision.  Left as single,
% every downstream product (residuals, log-likelihood) is computed in
% single precision, and a summed chi-square over thousands of nodes loses
% meaningful accuracy.  Cast to double once, here.
Sz   = double(Sz(:));
S_DA = double(S_DA(:));

% A NaN or Inf anywhere in the input silently propagates into the modeled
% profile and then into the likelihood.  Catch it now rather than after
% hours of sampling.
if any(~isfinite(Sz))
    error(['Sz contains %d non-finite value(s). Fix the DEM / stream ' ...
           'extraction before running the MCMC.'], sum(~isfinite(Sz)));
end
if any(~isfinite(S_DA))
    error(['S_DA contains %d non-finite value(s). Fix the drainage-area ' ...
           'grid before running the MCMC.'], sum(~isfinite(S_DA)));
end
% Drainage area enters as (1./S_DA).^(m/n); a zero or negative area gives
% Inf/complex values in the steady-state profile.
if any(S_DA <= 0)
    error(['S_DA contains %d non-positive value(s). Drainage area must ' ...
           'be > 0 at every stream node.'], sum(S_DA <= 0));
end

% Normalize elevations relative to outlet
Sz_norm = Sz - min(Sz);

% Stream data error (meters)
stream_err = 200;  % reflects model inadequacy
n_stream = length(Sz_norm);
n_cave   = length(cave_ages);

% Hoisted out of the MCMC loop: this vector is constant, so rebuilding it
% on every iteration is pure overhead.
stream_sigma = stream_err * ones(n_stream, 1);

fprintf('Data loaded: %d stream nodes, %d cave observations\n', n_stream, n_cave);

%% ========================================================================
%  SECTION 4: INITIALIZE MCMC
%  ========================================================================

total_iter = n_burnin + n_postburn;

% Make sure new MATLAB sessions use different random numbers (as in
% Gallen's master script: "rng shuffle").  Without this, every fresh
% MATLAB session would reproduce the identical chain.
rng('shuffle');

% Storage arrays
params      = zeros(total_iter, n_params);
logL_chain  = zeros(total_iter, 1);
accepted    = zeros(total_iter, 1);

% Initialize
params(1,:) = params_init;

% Run initial forward model.
% K is DERIVED from the sampled steepness: K = U_pre / ksn_ref^n.
K_init = params_init(1) / params_init(3)^params_init(4);
m_init = params_init(5) * params_init(4);  % m = (m/n) * n
fprintf('Derived K at start: %.3e  (ksn_ref=%.1f, n=%.2f)\n', ...
        K_init, params_init(3), params_init(4));

% Build rate/transition vectors for the generalized forward model
U_rates_init = [params_init(1), params_init(2)];
t_trans_init = params_init(6);

fprintf('Running initial forward model...\n');
Z_mod = hc_river_forward_model(S, S_DA, U_rates_init, t_trans_init, ...
    K_init, m_init, params_init(4), dt_forward);

% Initial cave predictions
cave_pred = cave_forward_model(cave_ages, U_rates_init, t_trans_init);

% Initial likelihood
[logL_init, ~, ~] = hc_loglikelihood(Sz_norm, Z_mod, ...
    stream_sigma, cave_heights, cave_pred, cave_height_err);

% If the starting model is already non-finite the whole chain is
% meaningless, so stop now instead of after hours of sampling.
if ~isfinite(logL_init)
    error(['Initial log-likelihood is %s. The forward model returned a ' ...
           'non-finite profile at the starting parameters -- check ' ...
           'params_init, prior_bounds and the stream data.'], ...
           num2str(logL_init));
end

logL_chain(1) = logL_init;
accepted(1)   = 1;

% Track MAP (maximum a posteriori)
logP_map = logL_init + logprior_hc(params_init, prior_bounds, cave_prior);
params_map = params_init;
Z_mod_map  = Z_mod;
cave_pred_map = cave_pred;

fprintf('Initial log-likelihood: %.2f\n', logL_init);

% --- Diagnostics ---
fprintf('\n--- Initial Model vs Observations ---\n');
fprintf('  Observed elevation range: [%.1f, %.1f] m  (relief = %.1f m)\n', ...
    min(Sz_norm), max(Sz_norm), max(Sz_norm) - min(Sz_norm));
fprintf('  Model elevation range:    [%.1f, %.1f] m  (relief = %.1f m)\n', ...
    min(Z_mod), max(Z_mod), max(Z_mod) - min(Z_mod));
fprintf('  RMS residual: %.1f m  (ratio to relief: %.1f%%)\n', ...
    sqrt(mean((Sz_norm - Z_mod).^2)), ...
    100 * sqrt(mean((Sz_norm - Z_mod).^2)) / max(max(Sz_norm), 1));
fprintf('  Cave predictions vs observed:\n');
for ci = 1:length(cave_ages)
    fprintf('    Cave %d: obs=%.0f m, mod=%.0f m, resid=%.0f m\n', ...
        ci, cave_heights(ci), cave_pred(ci), cave_heights(ci) - cave_pred(ci));
end
fprintf('  stream_err = %.0f m, n_stream = %d, n_cave = %d\n', ...
    stream_err, n_stream, n_cave);
fprintf('--- End Initial Diagnostics ---\n');

fprintf('\nStarting MCMC: %d burn-in + %d post-burn-in iterations\n', ...
    n_burnin, n_postburn);

%% ========================================================================
%  SECTION 5: RUN MCMC
%  ========================================================================

tic;
n_accept     = 0;
n_ffail      = 0;   % forward-model failures (caught and rejected)
n_nonfinite  = 0;   % candidates with non-finite likelihood (rejected)
accept_window = 0;  % accepts within the current tuning window
first_err_shown = false;  % print the first forward-model error only

for i = 2:total_iter
    % Progress report
    if mod(i, 1000) == 0
        elapsed = toc;
        rate = i / elapsed;
        remaining = (total_iter - i) / rate;
        fprintf('  Iter %d/%d (%.0f/s, ~%.0f s remaining, accept=%.1f%%)\n', ...
            i, total_iter, rate, remaining, 100*n_accept/(i-1));
    end

    % --- Adaptive step-size tuning during burn-in ---
    % Automates the manual step tuning Gallen does by trial and error.
    % Runs at the top of the iteration so it always fires on schedule,
    % regardless of prior-rejection "continue" statements below.
    if adapt_steps && i <= n_burnin && i > 2 && mod(i-1, tune_interval) == 0
        win_acc = accept_window / tune_interval;
        scale   = exp(win_acc - target_accept);   % >1 too high, <1 too low
        scale   = min(max(scale, 0.5), 2.0);      % clamp per-update change
        p_steps = p_steps .* scale;
        accept_window = 0;
        if mod(i-1, tune_interval*5) == 0
            fprintf('    [tune] iter %d: window accept=%.1f%%, step scale=%.3f\n', ...
                i-1, 100*win_acc, scale);
        end
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

    % Run forward model with candidate parameters.
    % K is derived from ksn_ref and n, never sampled directly.
    K_cand = candidate(1) / candidate(3)^candidate(4);
    m_cand = candidate(5) * candidate(4);  % m = (m/n) * n

    U_rates_cand = [candidate(1), candidate(2)];
    t_trans_cand = candidate(6);

    try
        Z_cand = hc_river_forward_model(S, S_DA, U_rates_cand, ...
            t_trans_cand, K_cand, m_cand, candidate(4), dt_forward);

        cave_cand = cave_forward_model(cave_ages, U_rates_cand, t_trans_cand);

        % Likelihood
        [logL_cand, ~, ~] = hc_loglikelihood(Sz_norm, Z_cand, ...
            stream_sigma, cave_heights, cave_cand, cave_height_err);
    catch ME
        % Forward model failed - reject.  Surface the FIRST error so a
        % systematic fault (e.g. a typo or a bad data field) cannot
        % silently reject every candidate for days on end.
        if ~first_err_shown
            first_err_shown = true;
            fprintf(['\n  NOTE: first forward-model failure at iter %d:\n' ...
                     '        %s\n        (further failures counted silently)\n\n'], ...
                     i, ME.message);
        end
        params(i,:) = current;
        logL_chain(i) = logL_chain(i-1);
        accepted(i) = 0;
        n_ffail = n_ffail + 1;
        continue;
    end

    % Guard against a non-finite likelihood BEFORE it reaches the
    % acceptance test.  MATLAB's min() omits NaN, so min(NaN, 0) returns
    % 0 and "log(rand) < 0" would accept unconditionally; the NaN would
    % then enter logL_chain and every later ratio would also be NaN,
    % turning the remainder of the run into an unconstrained random walk.
    if ~isfinite(logL_cand)
        params(i,:) = current;
        logL_chain(i) = logL_chain(i-1);
        accepted(i) = 0;
        n_nonfinite = n_nonfinite + 1;
        continue;
    end

    % Proposal probabilities (symmetric, so they cancel)
    lq_fwd = logproposal_hc(current, candidate, p_steps);
    lq_rev = logproposal_hc(candidate, current, p_steps);

    % Log acceptance ratio
    log_alpha = (lp_cand + logL_cand + lq_rev) - ...
                (lp_current + logL_chain(i-1) + lq_fwd);

    % Belt-and-braces: an undefined ratio must never be treated as a
    % favourable one (see the min()/NaN note above).
    if isnan(log_alpha)
        log_alpha = -Inf;
    end

    % Metropolis-Hastings accept/reject
    if log(rand) < min(log_alpha, 0)
        % Accept
        params(i,:) = candidate;
        logL_chain(i) = logL_cand;
        accepted(i) = 1;
        n_accept = n_accept + 1;
        accept_window = accept_window + 1;

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

% Acceptance split (burn-in is inflated/tuned; post-burn-in is the chain
% that actually matters for the posterior).
acc_post = mean(accepted(n_burnin+1:end));
fprintf('Post-burn-in acceptance rate: %.1f%%  (target ~%.0f%%)\n', ...
    100 * acc_post, 100 * target_accept);

if adapt_steps
    fprintf('Tuned step sizes (frozen after burn-in): %s\n', ...
        mat2str(p_steps, 3));
end

if n_ffail > 0
    fprintf('WARNING: forward model failed on %d/%d iterations (%.2f%%).\n', ...
        n_ffail, total_iter, 100 * n_ffail / total_iter);
    fprintf('         If this fraction is large, check parameter bounds / dt.\n');
end

if n_nonfinite > 0
    fprintf('WARNING: %d/%d candidates (%.2f%%) had a non-finite likelihood\n', ...
        n_nonfinite, total_iter, 100 * n_nonfinite / total_iter);
    fprintf('         and were rejected. A large fraction means the forward\n');
    fprintf('         model is overflowing somewhere in the sampled range.\n');
end

%% ========================================================================
%  SECTION 6: EXTRACT POST-BURN-IN RESULTS
%  ========================================================================

params_post = params(n_burnin+1:end, :);
logL_post   = logL_chain(n_burnin+1:end);

% Compute statistics.  Units in param_names must match param_scale:
% rates are displayed in mm/yr (scale 1e3) and t_capture in Ma (1e-6).
param_names = {'U_{pre} (mm/yr)', 'U_{post} (mm/yr)', 'k_{sn} relict', ...
               'n', 'm/n', 't_{capture} (Ma)'};
param_scale = [1e3, 1e3, 1, 1, 1, 1e-6];

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

% K is a DERIVED quantity: K = U_pre / ksn_ref^n, evaluated per sample so
% its posterior correctly propagates the U_pre / ksn_ref / n covariance.
K_post   = params_post(:,1) ./ params_post(:,3).^params_post(:,4);
K_map    = params_map(1) / params_map(3)^params_map(4);
logK_post = log10(K_post);

fprintf('\nDerived K = U_pre / ksn_ref^n :\n');
fprintf('  MAP    : %.3e   (log10 K = %.2f)\n', K_map, log10(K_map));
fprintf('  median : %.3e   (log10 K = %.2f)\n', ...
    median(K_post), median(logK_post));
fprintf('  68%% CI : [%.3e, %.3e]\n', ...
    prctile(K_post, 16), prctile(K_post, 84));
fprintf('  95%% CI : [%.3e, %.3e]\n', ...
    prctile(K_post, 2.5), prctile(K_post, 97.5));

% --- Independent cross-check against the ksn analysis -------------------
% Under uniform K and steady state in both reaches,
%       ksn_adjusted/ksn_relict = (U_post/U_pre)^(1/n)
% so the measured steepness ratio implies
%       n = ln(U_post/U_pre) / ln(ksn_adjusted/ksn_relict).
%
% This is recomputed from THIS run's posterior rates, per sample, rather
% than hardcoded: the implied n depends on U_post/U_pre, which the profile
% is free to move.  (Hardcoding it against the old cave-prior rates made
% the two estimates look far more discrepant than they are.)
ksn_adj_obs = 227.4;    % mean of 38 below-knickpoint segments
ksn_rel_obs = 113.3;    % mean of 38 above-knickpoint segments
ksn_rel_sd  = 33.8;

n_implied = log(params_post(:,2) ./ params_post(:,1)) ./ ...
            log(ksn_adj_obs / ksn_rel_obs);

fprintf('\nCross-check vs ksn_erosion_analysis:\n');
fprintf('  ksn-implied n (from posterior rates): %.2f [%.2f - %.2f]\n', ...
    median(n_implied), prctile(n_implied, 2.5), prctile(n_implied, 97.5));
fprintf('  profile posterior n                 : %.2f [%.2f - %.2f]\n', ...
    median(params_post(:,4)), prctile(params_post(:,4), 2.5), ...
    prctile(params_post(:,4), 97.5));
fprintf('  measured ksn relict                 : %.1f +/- %.1f\n', ...
    ksn_rel_obs, ksn_rel_sd);
fprintf('  posterior ksn_ref                   : %.1f [%.1f - %.1f]\n', ...
    median(params_post(:,3)), prctile(params_post(:,3), 2.5), ...
    prctile(params_post(:,3), 97.5));
fprintf(['  NOTE: ksn above was measured at TAK theta_ref = 0.45, whereas\n' ...
         '        the model normalizes at m/n = %.2f.  The two segments sit\n' ...
         '        at different drainage areas, so the RATIO is not invariant\n' ...
         '        to that choice -- treat this as an approximate check.\n'], ...
         median(params_post(:,5)));

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
    prior_bounds, cave_prior, param_names, param_scale, output_dir, fileTag, S_DA);

fprintf('\nDone. Run plot_hc_results.m for additional visualizations.\n');
