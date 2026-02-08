function plot_hc_results(params, logL_chain, n_burnin, params_map, ...
    Z_mod_map, Sz_obs, S, cave_ages, cave_heights, cave_height_err, ...
    cave_pred_map, prior_bounds, cave_prior, param_names, param_scale, ...
    output_dir, fileTag, S_DA)
% PLOT_HC_RESULTS  Diagnostic and summary plots for Hells Canyon MCMC inversion.
%
% Creates three figures:
%   Fig 1: MCMC diagnostics (chains, acceptance, likelihood)
%   Fig 2: Parameter posteriors with priors
%   Fig 3: Model fit (river profile + cave data)
%
% Inputs:
%   params       - [total_iter x n_params] full MCMC chain
%   logL_chain   - [total_iter x 1] log-likelihood chain
%   n_burnin     - Number of burn-in iterations
%   params_map   - [1 x n_params] MAP parameter vector
%   Z_mod_map    - Modeled stream elevations at MAP
%   Sz_obs       - Observed stream elevations (normalized)
%   S            - TopoToolbox STREAMobj
%   cave_ages    - Cave burial ages (years)
%   cave_heights - Cave heights above river (m)
%   cave_height_err - Cave height uncertainties (m)
%   cave_pred_map - MAP cave height predictions (m)
%   prior_bounds - [n_params x 2] prior bounds
%   cave_prior   - struct with informative prior settings
%   param_names  - cell array of parameter names
%   param_scale  - [1 x n_params] display scaling factors
%   output_dir   - Directory for saving figures
%   fileTag      - Tag for output filenames
%   S_DA         - (optional) Drainage area at stream nodes (m^2) for chi/tau

n_params = size(params, 2);
total_iter = size(params, 1);
params_post = params(n_burnin+1:end, :);

%% Figure 1: MCMC Diagnostics
fig1 = figure('Position', [50, 50, 1400, 800]);

% Trace plots
for j = 1:n_params
    subplot(n_params + 1, 2, 2*j - 1)
    plot(params(:,j) * param_scale(j), 'Color', [0.5 0.5 0.5 0.3]);
    hold on
    xline(n_burnin, 'r--', 'LineWidth', 1.5);
    yline(params_map(j) * param_scale(j), 'b-', 'LineWidth', 1);
    ylabel(param_names{j}, 'FontSize', 8);
    if j == 1; title('Parameter Traces'); end
    if j < n_params; set(gca, 'XTickLabel', []); end
end
xlabel('Iteration');

% Log-likelihood trace
subplot(n_params + 1, 2, 2*(n_params+1) - 1)
plot(logL_chain, 'Color', [0.5 0.5 0.5 0.3]); hold on
xline(n_burnin, 'r--', 'LineWidth', 1.5);
xlabel('Iteration'); ylabel('Log-likelihood');

% Posterior histograms
for j = 1:n_params
    subplot(n_params + 1, 2, 2*j)
    histogram(params_post(:,j) * param_scale(j), 50, ...
        'Normalization', 'pdf', 'FaceColor', [0.4 0.6 0.8], 'EdgeColor', 'none');
    hold on

    % MAP value
    xline(params_map(j) * param_scale(j), 'r-', 'LineWidth', 2);

    % Prior bounds
    xline(prior_bounds(j,1) * param_scale(j), 'k--');
    xline(prior_bounds(j,2) * param_scale(j), 'k--');

    % Informative prior (if applicable)
    if cave_prior.use_informative
        x_range = linspace(prior_bounds(j,1), prior_bounds(j,2), 200);
        switch j
            case 1  % U_pre
                if cave_prior.U_pre_std > 0
                    prior_pdf = normpdf(x_range, cave_prior.U_pre_mean, cave_prior.U_pre_std);
                    prior_pdf = prior_pdf / max(prior_pdf) * max(ylim) * 0.5;
                    plot(x_range * param_scale(j), prior_pdf, 'g-', 'LineWidth', 1.5);
                end
            case 2  % U_post
                if cave_prior.U_post_std > 0
                    prior_pdf = normpdf(x_range, cave_prior.U_post_mean, cave_prior.U_post_std);
                    prior_pdf = prior_pdf / max(prior_pdf) * max(ylim) * 0.5;
                    plot(x_range * param_scale(j), prior_pdf, 'g-', 'LineWidth', 1.5);
                end
            case 6  % t_capture
                if cave_prior.t_capture_std > 0
                    prior_pdf = normpdf(x_range, cave_prior.t_capture_mean, cave_prior.t_capture_std);
                    prior_pdf = prior_pdf / max(prior_pdf) * max(ylim) * 0.5;
                    plot(x_range * param_scale(j), prior_pdf, 'g-', 'LineWidth', 1.5);
                end
        end
    end

    xlabel(param_names{j});
    if j == 1; title('Posterior Distributions'); end
end

% Acceptance rate over time (detect changes in any parameter as proxy)
subplot(n_params + 1, 2, 2*(n_params+1))
win = 1000;
any_changed = any(diff(params, 1, 1) ~= 0, 2);
accept_rate = movmean(any_changed, win);
plot(accept_rate * 100, 'k-');
xlabel('Iteration'); ylabel('Accept Rate (%)');
yline(25, 'r--'); yline(50, 'r--');
title('Running Acceptance Rate');

sgtitle(sprintf('MCMC Diagnostics: %s', fileTag), 'FontSize', 14);
saveas(fig1, fullfile(output_dir, ['diagnostics_' fileTag '.png']));

%% Figure 2: Key posterior pairs and correlations
fig2 = figure('Position', [100, 100, 1200, 800]);

% K vs n (the key trade-off)
subplot(2,3,1)
scatter(params_post(:,3), params_post(:,4), 2, logL_chain(n_burnin+1:end), ...
    'filled', 'MarkerFaceAlpha', 0.1);
hold on
plot(params_map(3), params_map(4), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
xlabel('log_{10}(K)'); ylabel('n');
title('K vs n Trade-off');
colorbar; colormap(parula);

% t_capture vs U_post
subplot(2,3,2)
scatter(params_post(:,6) * 1e-6, params_post(:,2) * 1e3, 2, ...
    logL_chain(n_burnin+1:end), 'filled', 'MarkerFaceAlpha', 0.1);
hold on
plot(params_map(6)*1e-6, params_map(2)*1e3, 'rp', 'MarkerSize', 15, ...
    'MarkerFaceColor', 'r');
xlabel('t_{capture} (Ma)'); ylabel('U_{post} (mm/yr)');
title('Capture Time vs Post-Capture Rate');

% t_capture histogram with cave prior
subplot(2,3,3)
histogram(params_post(:,6) * 1e-6, 50, 'Normalization', 'pdf', ...
    'FaceColor', [0.8 0.4 0.4], 'EdgeColor', 'none'); hold on
if cave_prior.use_informative
    x = linspace(0.5, 5, 200);
    y = normpdf(x, cave_prior.t_capture_mean*1e-6, cave_prior.t_capture_std*1e-6);
    y = y / max(y) * max(ylim) * 0.7;
    plot(x, y, 'g-', 'LineWidth', 2);
    legend('Posterior', 'Cave Prior');
end
xlabel('t_{capture} (Ma)'); ylabel('PDF');
title('Capture Timing Posterior');

% n histogram
subplot(2,3,4)
histogram(params_post(:,4), 50, 'Normalization', 'pdf', ...
    'FaceColor', [0.4 0.7 0.4], 'EdgeColor', 'none'); hold on
xline(1, 'k--', 'n=1', 'LineWidth', 2, 'FontSize', 12);
xline(median(params_post(:,4)), 'r-', 'LineWidth', 2);
xlabel('n (slope exponent)'); ylabel('PDF');
title('Slope Exponent Posterior');

% K histogram
subplot(2,3,5)
histogram(params_post(:,3), 50, 'Normalization', 'pdf', ...
    'FaceColor', [0.7 0.4 0.7], 'EdgeColor', 'none');
xlabel('log_{10}(K)'); ylabel('PDF');
title('Erodibility Posterior');

% U_pre vs U_post
subplot(2,3,6)
scatter(params_post(:,1)*1e3, params_post(:,2)*1e3, 2, ...
    'filled', 'MarkerFaceAlpha', 0.1);
hold on
plot(params_map(1)*1e3, params_map(2)*1e3, 'rp', 'MarkerSize', 15, ...
    'MarkerFaceColor', 'r');
plot([0 0.5], [0 0.5], 'k--');
xlabel('U_{pre} (mm/yr)'); ylabel('U_{post} (mm/yr)');
title('Pre vs Post Capture Rates');

sgtitle(sprintf('Parameter Posteriors: %s', fileTag), 'FontSize', 14);
saveas(fig2, fullfile(output_dir, ['posteriors_' fileTag '.png']));

%% Figure 3: Model Fit
fig3 = figure('Position', [150, 150, 1200, 600]);

% River profile fit
subplot(1,2,1)
if nargin >= 18 && ~isempty(S_DA)
    Stau_map = compute_tau_from_chi(S, S_DA, params_map);
else
    Stau_map = [];
end
if ~isempty(Stau_map)
    plot(Stau_map / 1e6, Sz_obs, '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 3); hold on
    plot(Stau_map / 1e6, Z_mod_map, 'k.', 'MarkerSize', 3);
    xlabel('\tau (Ma)'); ylabel('Elevation (m)');
else
    plot(Sz_obs, '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 3); hold on
    plot(Z_mod_map, 'k.', 'MarkerSize', 3);
    xlabel('Node index'); ylabel('Elevation (m)');
end
legend('Observed', 'MAP Model', 'Location', 'northwest');
title('River Profile Fit');

% Cave fit
subplot(1,2,2)
errorbar(cave_ages / 1e6, cave_heights, cave_height_err, 'bo', ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.6 0.9], 'LineWidth', 1.5);
hold on
plot(cave_ages / 1e6, cave_pred_map, 'rs', 'MarkerSize', 10, ...
    'MarkerFaceColor', [0.9 0.3 0.3], 'LineWidth', 1.5);

% Plot model incision history
t_plot = linspace(0, max(cave_ages)*1.1, 200);
h_plot = cave_forward_model(t_plot, params_map(6), params_map(1), params_map(2));
plot(t_plot / 1e6, h_plot, 'r-', 'LineWidth', 2);

% Mark capture time
xline(params_map(6) / 1e6, 'k--', sprintf('t_{cap}=%.1f Ma', params_map(6)/1e6), ...
    'LineWidth', 1.5, 'FontSize', 11);

xlabel('Age (Ma)'); ylabel('Height above river (m)');
legend('Cave data', 'MAP predictions', 'Incision model', 'Location', 'northwest');
title('Cave Data Fit');
set(gca, 'XDir', 'reverse');

sgtitle(sprintf('Model Fit: %s', fileTag), 'FontSize', 14);
saveas(fig3, fullfile(output_dir, ['model_fit_' fileTag '.png']));

fprintf('Figures saved to: %s\n', output_dir);

end

%% Helper function
function Stau = compute_tau_from_chi(S, S_DA, params_map)
% Compute tau = chi / K using MAP parameters for plotting.
% chi = integral of (1/A)^(m/n) dx from outlet upstream.
try
    K_map  = 10^params_map(3);
    mn_map = params_map(5);  % m/n

    % Compute chi via upstream integration
    Schi = zeros(size(S.distance));
    Six  = S.ix;
    Sixc = S.ixc;
    Sx   = S.distance;
    Sa   = (1 ./ S_DA).^mn_map;

    for lp = numel(Six):-1:1
        Schi(Six(lp)) = Schi(Sixc(lp)) + ...
            (Sa(Sixc(lp)) + (Sa(Six(lp)) - Sa(Sixc(lp))) / 2) * ...
            abs(Sx(Sixc(lp)) - Sx(Six(lp)));
    end

    Stau = Schi / K_map;
catch
    Stau = [];
end
end
