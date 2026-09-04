function cfg = hc_basin_defaults()
% HC_BASIN_DEFAULTS  Settings shared by EVERY basin in a multi-basin run.
%
% This is the single place to set the things that should be identical
% across basins -- burn-in length, timestep, proposal tuning, prior bounds,
% likelihood weighting.  Per-basin overrides live in hc_basin_list.m and
% are applied on top of these.
%
% Keeping these global matters scientifically as much as practically: if
% each basin were tuned separately, a spread in recovered t_capture could
% just reflect a spread in settings.  Held fixed, cross-basin agreement in
% t_capture becomes evidence for a regional capture event.
%
% Usage:
%   cfg = hc_basin_defaults();
%   cfg.n_postburn = 5e4;          % override globally
%
% See also: hc_basin_list, run_hc_inversion, run_all_basins

cfg = struct();

%% --- MCMC settings ------------------------------------------------------
cfg.n_burnin   = 2e4;
cfg.n_postburn = 2e5;

% dt = 100,000 validated against a 6,250 yr reference at n ~ 6 (0.36-0.41 m
% RMS against a stream_err of 200 m).  Re-run dt_sensitivity if the
% posterior for a new basin lands somewhere very different, since knickpoint
% celerity scales with K and the required step scales with n.
cfg.dt_forward = 100000;

%% --- Proposal ----------------------------------------------------------
cfg.adapt_steps   = true;
cfg.tune_interval = 200;
cfg.target_accept = 0.30;

% 'diagonal' | 'covariance'.  For a first pass over new basins use
% 'diagonal'; switch to 'covariance' per basin once each has a pilot chain.
cfg.prop_mode  = 'diagonal';
cfg.pilot_file = '';

%% --- Likelihood weighting ----------------------------------------------
cfg.stream_err   = 200;   % metres; reflects model inadequacy, not DEM error
cfg.N_eff_stream = 100;   % effective independent stream constraints

%% --- Parameter layout [7] ----------------------------------------------
% [1] U_pre  [2] U_post  [3] ksn_ref  [4] n  [5] m/n  [6] t_capture  [7] U_grad
% K is DERIVED as U_pre / ksn_ref^n -- never sampled.  See hc_derive_K.
cfg.prior_bounds = [
    1e-6,   5e-4;     % U_pre     (m/yr, at the outlet)
    1e-5,   5e-3;     % U_post    (m/yr, at the outlet)
    20,     500;      % ksn_ref   relict channel steepness
    0.5,    12;       % n         slope exponent
    0.3,    0.8;      % m/n       concavity
    0.5e6,  5e6;      % t_capture (yr BP)
   -0.99,   2.0;      % U_grad    0 = spatially uniform uplift
];

% Starting values.  Per-basin overrides are strongly recommended once you
% have a sense of each basin -- starting near the solution is what Gallen
% does, and it shortens burn-in considerably.
cfg.params_init = [4.9e-6, 1.53e-4, 142.6, 5.55, 0.595, 2.40e6, 0.0];

% Only the RELATIVE proportions matter; the burn-in tuner rescales the
% whole vector by one shared scalar.  These are ~0.4x the posterior sigma
% measured on the trunk basin.
cfg.p_steps = [1.1e-6, 5.4e-6, 17.0, 0.40, 0.012, 1.2e5, 0.10];

%% --- Cave-derived priors -----------------------------------------------
% Applied to EVERY basin, including tributaries with no caves of their own.
% That is defensible here: these streams are graded to the Snake, so their
% base-level history IS the Snake's incision history, which is what the
% caves record.  Without them, U_pre and U_post are unidentifiable in a
% cave-less basin -- the profile alone cannot separate U from K, since
% relict steepness is (U_pre/K)^(1/n) and, with ksn_ref sampled, U_pre and
% K trade off exactly.
cfg.cave_prior = struct( ...
    'use_informative', true, ...
    't_capture_mean',  2.1e6, 't_capture_std', 1.0e6, ...
    'U_pre_mean',      1e-5,  'U_pre_std',     5e-6,  ...
    'U_post_mean',     1.25e-4, 'U_post_std',  3.5e-5, ...
    'ksn_ref_mean',    113.3, 'ksn_ref_std',   60.0,  ...
    'n_mean',          3.7,   'n_std',         0);
% n_std = 0 disables the prior on n, leaving it uniform within bounds so
% the profile determines it and the ksn-derived estimate stays an
% independent check rather than an assumption.

%% --- Output ------------------------------------------------------------
cfg.output_dir = fullfile(pwd, 'results_multibasin');
cfg.rng_seed   = [];      % [] = shuffle; set an integer for reproducibility
cfg.verbose    = true;    % per-iteration progress (noisy under parfor)

end
