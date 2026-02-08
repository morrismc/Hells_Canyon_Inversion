function results = ksn_erosion_analysis(ksn_file, varargin)
% KSN_EROSION_ANALYSIS  Estimate erodibility K and slope exponent n from
% ksn observations and cave-derived erosion rates.
%
% Uses the steady-state stream power relationship: E = K * ksn^n
% Combined with cave-derived incision rates for above/below knickpoint
% segments to constrain both K and n via Bayesian Monte Carlo.
%
% Inputs:
%   ksn_file  - path to ksnTable.xlsx (columns: Basin, Segment, ksn, ...)
%               Segment=1 is below knickpoint (adjusted), Segment=2 is above
%
%   Optional name-value pairs:
%     'E_relict'      - Pre-capture erosion rate [mm/yr] (default: 0.01)
%     'E_relict_err'  - 1-sigma uncertainty on E_relict (default: 0.005)
%     'E_adjusted'    - Post-capture erosion rate [mm/yr] (default: 0.125)
%     'E_adjusted_err'- 1-sigma uncertainty on E_adjusted (default: 0.035)
%     'mn'            - m/n reference concavity (default: 0.45)
%     'n_samples'     - Number of Monte Carlo samples (default: 100000)
%     'n_prior'       - Prior bounds on n [low, high] (default: [0.5, 10])
%     'logK_prior'    - Prior bounds on log10(K) [low, high] (default: [-10, -4])
%     'output_dir'    - Directory for output files (default: pwd)
%
% Outputs:
%   results - structure with fields:
%     .K_map, .n_map         - MAP estimates
%     .K_median, .n_median   - Median estimates
%     .K_ci68, .n_ci68       - 68% credible intervals
%     .K_ci95, .n_ci95       - 95% credible intervals
%     .ksn_relict, .ksn_adjusted  - Mean ksn for each segment
%     .n_uniform_K            - n implied by uniform-K assumption
%     .samples                - [n_accepted x 2] matrix of [logK, n] samples
%
% Reference:
%   Gallen (2018), EPSL: E = K * ksn^n at steady state
%   Morriss et al. (2025), PNAS: Cave-derived incision rates
%
% Author: Adapted for Hells Canyon from Gallen (2018) framework

%% Parse inputs
p = inputParser;
addRequired(p, 'ksn_file', @ischar);
addParameter(p, 'E_relict', 0.01, @isscalar);         % mm/yr
addParameter(p, 'E_relict_err', 0.005, @isscalar);     % mm/yr
addParameter(p, 'E_adjusted', 0.125, @isscalar);       % mm/yr
addParameter(p, 'E_adjusted_err', 0.035, @isscalar);   % mm/yr
addParameter(p, 'mn', 0.45, @isscalar);
addParameter(p, 'n_samples', 100000, @isscalar);
addParameter(p, 'n_prior', [0.5, 10], @(x) numel(x)==2);
addParameter(p, 'logK_prior', [-10, -4], @(x) numel(x)==2);
addParameter(p, 'output_dir', pwd, @ischar);

parse(p, ksn_file, varargin{:});
opts = p.Results;

%% Load ksn data
fprintf('Loading ksn data from: %s\n', ksn_file);
T = readtable(ksn_file);

% Extract ksn values by segment
% Segment 1 = below knickpoint (adjusted/responding to capture)
% Segment 2 = above knickpoint (relict/pre-capture)
idx_adjusted = T.Segment == 1;
idx_relict   = T.Segment == 2;

ksn_adj = T.ksn(idx_adjusted);
ksn_rel = T.ksn(idx_relict);

% Remove NaN/zero values
ksn_adj = ksn_adj(ksn_adj > 0 & ~isnan(ksn_adj));
ksn_rel = ksn_rel(ksn_rel > 0 & ~isnan(ksn_rel));

fprintf('  Adjusted segments (below KP): n=%d, mean ksn=%.1f, std=%.1f\n', ...
    length(ksn_adj), mean(ksn_adj), std(ksn_adj));
fprintf('  Relict segments (above KP):   n=%d, mean ksn=%.1f, std=%.1f\n', ...
    length(ksn_rel), mean(ksn_rel), std(ksn_rel));

%% Quick diagnostic: implied n assuming uniform K
ksn_ratio = mean(ksn_adj) / mean(ksn_rel);
E_ratio_low  = (opts.E_adjusted - opts.E_adjusted_err) / ...
               (opts.E_relict + opts.E_relict_err);
E_ratio_mid  = opts.E_adjusted / opts.E_relict;
E_ratio_high = (opts.E_adjusted + opts.E_adjusted_err) / ...
               max(opts.E_relict - opts.E_relict_err, 1e-4);

n_implied_low  = log(E_ratio_low)  / log(ksn_ratio);
n_implied_mid  = log(E_ratio_mid)  / log(ksn_ratio);
n_implied_high = log(E_ratio_high) / log(ksn_ratio);

fprintf('\n--- Diagnostic: Implied n (uniform K assumption) ---\n');
fprintf('  ksn ratio (adjusted/relict): %.2f\n', ksn_ratio);
fprintf('  E ratio range: %.1f to %.1f\n', E_ratio_low, E_ratio_high);
fprintf('  Implied n range: %.1f to %.1f (central: %.1f)\n', ...
    n_implied_low, n_implied_high, n_implied_mid);
if n_implied_mid > 3
    fprintf('  WARNING: n > 3 under uniform K - consider whether K varies spatially\n');
end

%% Bayesian Monte Carlo estimation of K and n
% Model: E = K * ksn^n
% log(E) = log(K) + n * log(ksn)
%
% We treat each ksn observation as providing a constraint:
%   For adjusted segments: log(E_adj) = log(K) + n * log(ksn_adj_i)
%   For relict segments:   log(E_rel) = log(K) + n * log(ksn_rel_i)
%
% with E_adj and E_rel drawn from their cave-constrained distributions.

fprintf('\nRunning Bayesian Monte Carlo (%d samples)...\n', opts.n_samples);

% Convert erosion rates to m/yr for K units consistency
E_adj_m = opts.E_adjusted * 1e-3;    % m/yr
E_rel_m = opts.E_relict * 1e-3;      % m/yr
E_adj_err_m = opts.E_adjusted_err * 1e-3;
E_rel_err_m = opts.E_relict_err * 1e-3;

% Collect all observations
n_adj = length(ksn_adj);
n_rel = length(ksn_rel);

% Preallocate accepted samples
max_accept = opts.n_samples;
accepted_logK = zeros(max_accept, 1);
accepted_n    = zeros(max_accept, 1);
accepted_logL = zeros(max_accept, 1);
n_accepted = 0;

% Fractional uncertainty on ksn observations (15% relative)
ksn_frac_err = 0.15;

for i = 1:opts.n_samples
    % Draw candidate parameters from prior
    logK_cand = opts.logK_prior(1) + ...
                (opts.logK_prior(2) - opts.logK_prior(1)) * rand;
    n_cand    = opts.n_prior(1) + ...
                (opts.n_prior(2) - opts.n_prior(1)) * rand;

    K_cand = 10^logK_cand;

    % Draw erosion rates from their distributions (truncated normal, E>0)
    E_adj_draw = -1;
    while E_adj_draw <= 0
        E_adj_draw = E_adj_m + E_adj_err_m * randn;
    end
    E_rel_draw = -1;
    while E_rel_draw <= 0
        E_rel_draw = E_rel_m + E_rel_err_m * randn;
    end

    % Compute log-likelihood
    logL = 0;

    % Adjusted segments: E_adj = K * ksn_adj^n
    for j = 1:n_adj
        % Draw ksn with uncertainty
        ksn_draw = ksn_adj(j) * (1 + ksn_frac_err * randn);
        if ksn_draw <= 0; ksn_draw = ksn_adj(j); end

        E_pred = K_cand * ksn_draw^n_cand;
        if E_pred <= 0; continue; end

        % Log-space likelihood (errors are multiplicative)
        log_resid = (log(E_adj_draw) - log(E_pred))^2;
        sigma_log = sqrt((E_adj_err_m/E_adj_draw)^2 + (ksn_frac_err * n_cand)^2);
        logL = logL - 0.5 * log_resid / sigma_log^2;
    end

    % Relict segments: E_rel = K * ksn_rel^n
    for j = 1:n_rel
        ksn_draw = ksn_rel(j) * (1 + ksn_frac_err * randn);
        if ksn_draw <= 0; ksn_draw = ksn_rel(j); end

        E_pred = K_cand * ksn_draw^n_cand;
        if E_pred <= 0; continue; end

        log_resid = (log(E_rel_draw) - log(E_pred))^2;
        sigma_log = sqrt((E_rel_err_m/E_rel_draw)^2 + (ksn_frac_err * n_cand)^2);
        logL = logL - 0.5 * log_resid / sigma_log^2;
    end

    % Accept all samples (importance sampling - we'll weight later)
    n_accepted = n_accepted + 1;
    accepted_logK(n_accepted) = logK_cand;
    accepted_n(n_accepted)    = n_cand;
    accepted_logL(n_accepted) = logL;
end

% Trim to actual accepted count
accepted_logK = accepted_logK(1:n_accepted);
accepted_n    = accepted_n(1:n_accepted);
accepted_logL = accepted_logL(1:n_accepted);

% Importance weights (normalized likelihood)
log_weights = accepted_logL - max(accepted_logL);  % for numerical stability
weights = exp(log_weights);
weights = weights / sum(weights);

%% Compute weighted statistics
% MAP estimate
[~, map_idx] = max(accepted_logL);
logK_map = accepted_logK(map_idx);
n_map    = accepted_n(map_idx);

% Weighted percentiles using resampling
n_resamp = 100000;
resamp_idx = randsample(n_accepted, n_resamp, true, weights);
logK_resamp = accepted_logK(resamp_idx);
n_resamp_vals = accepted_n(resamp_idx);

logK_median = median(logK_resamp);
n_median    = median(n_resamp_vals);

logK_ci68 = prctile(logK_resamp, [16, 84]);
logK_ci95 = prctile(logK_resamp, [2.5, 97.5]);
n_ci68    = prctile(n_resamp_vals, [16, 84]);
n_ci95    = prctile(n_resamp_vals, [2.5, 97.5]);

%% Also compute K assuming n=1 (for comparison with linear inversion)
% K = E / ksn when n=1
K_n1_adj = E_adj_m ./ ksn_adj;
K_n1_rel = E_rel_m ./ ksn_rel;

K_n1_adj_mean = mean(K_n1_adj);
K_n1_rel_mean = mean(K_n1_rel);

fprintf('\n========== RESULTS ==========\n');
fprintf('\n--- K and n estimates (Bayesian, uniform K assumed) ---\n');
fprintf('  MAP:    log10(K) = %.2f,  n = %.2f\n', logK_map, n_map);
fprintf('  Median: log10(K) = %.2f,  n = %.2f\n', logK_median, n_median);
fprintf('  K (MAP):    %.2e m^(1-2m)/yr\n', 10^logK_map);
fprintf('  K (median): %.2e m^(1-2m)/yr\n', 10^logK_median);
fprintf('  n 68%% CI: [%.2f, %.2f]\n', n_ci68);
fprintf('  n 95%% CI: [%.2f, %.2f]\n', n_ci95);
fprintf('  K 68%% CI: [%.2e, %.2e]\n', 10.^logK_ci68);
fprintf('  K 95%% CI: [%.2e, %.2e]\n', 10.^logK_ci95);

fprintf('\n--- K assuming n=1 (for linear inversion comparison) ---\n');
fprintf('  K from adjusted landscape: %.2e m^0.1/yr  (mean of %d obs)\n', ...
    K_n1_adj_mean, n_adj);
fprintf('  K from relict landscape:   %.2e m^0.1/yr  (mean of %d obs)\n', ...
    K_n1_rel_mean, n_rel);
fprintf('  K ratio (adjusted/relict): %.1f\n', K_n1_adj_mean / K_n1_rel_mean);
if K_n1_adj_mean / K_n1_rel_mean > 2
    fprintf('  NOTE: K_adj >> K_rel suggests either K varies or n != 1\n');
end

fprintf('\n--- Implied n (uniform K diagnostic) ---\n');
fprintf('  n = %.1f (%.1f to %.1f)\n', n_implied_mid, n_implied_low, n_implied_high);

%% Package results
results = struct();
results.K_map     = 10^logK_map;
results.n_map     = n_map;
results.logK_map  = logK_map;
results.K_median  = 10^logK_median;
results.n_median  = n_median;
results.logK_median = logK_median;
results.K_ci68    = 10.^logK_ci68;
results.K_ci95    = 10.^logK_ci95;
results.n_ci68    = n_ci68;
results.n_ci95    = n_ci95;
results.ksn_adjusted_mean = mean(ksn_adj);
results.ksn_adjusted_std  = std(ksn_adj);
results.ksn_relict_mean   = mean(ksn_rel);
results.ksn_relict_std    = std(ksn_rel);
results.n_uniform_K       = n_implied_mid;
results.n_uniform_K_range = [n_implied_low, n_implied_high];
results.K_n1_adjusted     = K_n1_adj_mean;
results.K_n1_relict       = K_n1_rel_mean;
results.samples           = [logK_resamp, n_resamp_vals];
results.E_adjusted        = opts.E_adjusted;
results.E_relict          = opts.E_relict;

%% Plot results
fig = figure('Position', [100, 100, 1400, 900]);

% Panel 1: ksn distributions
subplot(2,3,1)
histogram(ksn_rel, 'Normalization', 'pdf', 'FaceColor', [0.3 0.6 0.9], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.7); hold on
histogram(ksn_adj, 'Normalization', 'pdf', 'FaceColor', [0.9 0.3 0.3], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.7);
xlabel('k_{sn} (m^{0.9})'); ylabel('PDF');
legend('Relict (above KP)', 'Adjusted (below KP)');
title('Channel Steepness by Segment');

% Panel 2: log(E) vs log(ksn) with n-slope
subplot(2,3,2)
ksn_plot = logspace(log10(min([ksn_rel; ksn_adj])*0.5), ...
                    log10(max([ksn_rel; ksn_adj])*2), 100);
% Plot data points
plot(mean(ksn_rel), opts.E_relict, 'bs', 'MarkerSize', 12, ...
    'MarkerFaceColor', [0.3 0.6 0.9], 'LineWidth', 1.5); hold on
plot(mean(ksn_adj), opts.E_adjusted, 'rs', 'MarkerSize', 12, ...
    'MarkerFaceColor', [0.9 0.3 0.3], 'LineWidth', 1.5);
% Error bars
errorbar(mean(ksn_rel), opts.E_relict, opts.E_relict_err, opts.E_relict_err, ...
    std(ksn_rel), std(ksn_rel), 'b', 'LineWidth', 1.5);
errorbar(mean(ksn_adj), opts.E_adjusted, opts.E_adjusted_err, opts.E_adjusted_err, ...
    std(ksn_adj), std(ksn_adj), 'r', 'LineWidth', 1.5);
% Best-fit lines
K_fit = 10^logK_median;
E_n1   = K_n1_adj_mean * ksn_plot;
E_nfit = K_fit * ksn_plot.^n_median;
plot(ksn_plot, E_n1*1e3, 'k--', 'LineWidth', 1);
plot(ksn_plot, E_nfit*1e3, 'r-', 'LineWidth', 2);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('k_{sn} (m^{0.9})'); ylabel('E (mm/yr)');
legend('Relict', 'Adjusted', '', '', ...
    sprintf('n=1 (K=%.1e)', K_n1_adj_mean), ...
    sprintf('n=%.1f (K=%.1e)', n_median, K_fit), ...
    'Location', 'northwest');
title('Erosion Rate vs k_{sn}');

% Panel 3: Joint posterior of logK and n
subplot(2,3,3)
% 2D histogram
[N, xedges, yedges] = histcounts2(logK_resamp, n_resamp_vals, 50);
imagesc(xedges(1:end-1), yedges(1:end-1), N'); hold on
set(gca, 'YDir', 'normal');
colormap(gca, flipud(bone));
plot(logK_map, n_map, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
xlabel('log_{10}(K)'); ylabel('n');
title('Joint Posterior: K and n');
colorbar;

% Panel 4: Marginal posterior for n
subplot(2,3,4)
histogram(n_resamp_vals, 50, 'Normalization', 'pdf', ...
    'FaceColor', [0.4 0.7 0.4], 'EdgeColor', 'none'); hold on
xline(1, 'k--', 'n=1', 'LineWidth', 2, 'FontSize', 12);
xline(n_median, 'r-', sprintf('n=%.1f', n_median), 'LineWidth', 2, 'FontSize', 12);
xline(n_ci68(1), 'r--', 'LineWidth', 1);
xline(n_ci68(2), 'r--', 'LineWidth', 1);
xlabel('n (slope exponent)'); ylabel('PDF');
title('Marginal Posterior: n');

% Panel 5: Marginal posterior for K
subplot(2,3,5)
histogram(logK_resamp, 50, 'Normalization', 'pdf', ...
    'FaceColor', [0.7 0.4 0.7], 'EdgeColor', 'none'); hold on
xline(logK_median, 'r-', sprintf('K=%.1e', 10^logK_median), 'LineWidth', 2, 'FontSize', 12);
xline(logK_ci68(1), 'r--', 'LineWidth', 1);
xline(logK_ci68(2), 'r--', 'LineWidth', 1);
xline(log10(K_n1_adj_mean), 'b--', 'K (n=1, adj)', 'LineWidth', 1.5);
xline(log10(K_n1_rel_mean), 'c--', 'K (n=1, rel)', 'LineWidth', 1.5);
xlabel('log_{10}(K)'); ylabel('PDF');
title('Marginal Posterior: K');

% Panel 6: n sensitivity analysis
subplot(2,3,6)
n_test = linspace(0.5, 6, 100);
K_test = zeros(size(n_test));
for ii = 1:length(n_test)
    % For each n, compute best-fit K using both data points
    K_from_adj = E_adj_m / mean(ksn_adj)^n_test(ii);
    K_from_rel = E_rel_m / mean(ksn_rel)^n_test(ii);
    K_test(ii) = sqrt(K_from_adj * K_from_rel);  % geometric mean
end
plot(n_test, log10(K_test), 'k-', 'LineWidth', 2); hold on
% Mark the n=1 and best-fit points
K_at_n1 = sqrt((E_adj_m/mean(ksn_adj)) * (E_rel_m/mean(ksn_rel)));
plot(1, log10(K_at_n1), 'bs', 'MarkerSize', 12, 'MarkerFaceColor', [0.3 0.6 0.9]);
plot(n_median, logK_median, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
xlabel('n'); ylabel('log_{10}(K)');
title('K-n Trade-off');
legend('Geometric mean K', 'n=1', sprintf('Best fit (n=%.1f)', n_median), ...
    'Location', 'northeast');
grid on

sgtitle('Hells Canyon: k_{sn} vs Erosion Rate Analysis', 'FontSize', 14);

% Save figure
if ~exist(opts.output_dir, 'dir')
    mkdir(opts.output_dir);
end
saveas(fig, fullfile(opts.output_dir, 'ksn_erosion_analysis.png'));
fprintf('\nFigure saved to: %s\n', fullfile(opts.output_dir, 'ksn_erosion_analysis.png'));

end
