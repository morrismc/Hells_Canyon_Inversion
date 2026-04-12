function plot_hc_results(params, logL_chain, n_burnin, params_map, ...
    Z_mod_map, Sz_obs, S, cave_ages, cave_heights, cave_height_err, ...
    cave_pred_map, prior_bounds, cave_prior, param_names, param_scale, ...
    output_dir, fileTag, S_DA)
% PLOT_HC_RESULTS  Diagnostic and summary plots for Hells Canyon MCMC
% inversion (3-phase model with 8 parameters).
%
% Creates three figures:
%   Fig 1: MCMC diagnostics (chains, acceptance, likelihood)
%   Fig 2: Parameter posteriors with priors
%   Fig 3: Model fit (river profile + cave data)
%
% Supports both 6-parameter (2-phase) and 8-parameter (3-phase) models.

n_params = size(params, 2);
total_iter = size(params, 1);
params_post = params(n_burnin+1:end, :);

% Detect model type from parameter count
is_3phase = (n_params >= 8);

% Parameter indices for key quantities (adapts to 6- or 8-param layout)
if is_3phase
    idx_K  = 4;  idx_n  = 5;  idx_mn = 6;
    idx_t1 = 7;  idx_t2 = 8;
    idx_Upre = 1; idx_Upost = 3;
else
    idx_K  = 3;  idx_n  = 4;  idx_mn = 5;
    idx_t1 = 6;  idx_t2 = [];
    idx_Upre = 1; idx_Upost = 2;
end

%% Figure 1: MCMC Diagnostics
fig1 = figure('Position', [50, 50, 1400, max(800, 100*n_params)]);

% Trace plots
for j = 1:n_params
    subplot(n_params + 1, 2, 2*j - 1)
    plot(params(:,j) * param_scale(j), 'Color', [0.5 0.5 0.5 0.3]);
    hold on
    xline(n_burnin, 'r--', 'LineWidth', 1.5);
    yline(params_map(j) * param_scale(j), 'b-', 'LineWidth', 1);
    ylabel(param_names{j}, 'FontSize', 7);
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

    % Overlay informative priors (green)
    if cave_prior.use_informative
        x_range = linspace(prior_bounds(j,1), prior_bounds(j,2), 200);
        prior_pdf = [];
        if j == idx_Upre && isfield(cave_prior, 'U_pre_mean') && cave_prior.U_pre_std > 0
            prior_pdf = normpdf(x_range, cave_prior.U_pre_mean, cave_prior.U_pre_std);
        elseif j == idx_Upost && isfield(cave_prior, 'U_post_mean') && cave_prior.U_post_std > 0
            prior_pdf = normpdf(x_range, cave_prior.U_post_mean, cave_prior.U_post_std);
        elseif j == idx_n && isfield(cave_prior, 'n_mean') && isfield(cave_prior, 'n_std') && cave_prior.n_std > 0
            prior_pdf = normpdf(x_range, cave_prior.n_mean, cave_prior.n_std);
        elseif j == idx_t1
            if isfield(cave_prior, 't1_mean') && isfield(cave_prior, 't1_std') && cave_prior.t1_std > 0
                prior_pdf = normpdf(x_range, cave_prior.t1_mean, cave_prior.t1_std);
            elseif isfield(cave_prior, 't_capture_mean') && isfield(cave_prior, 't_capture_std') && cave_prior.t_capture_std > 0
                prior_pdf = normpdf(x_range, cave_prior.t_capture_mean, cave_prior.t_capture_std);
            end
        end
        if ~isempty(prior_pdf) && max(prior_pdf) > 0
            prior_pdf = prior_pdf / max(prior_pdf) * max(ylim) * 0.5;
            plot(x_range * param_scale(j), prior_pdf, 'g-', 'LineWidth', 1.5);
        end
    end

    xlabel(param_names{j}, 'FontSize', 7);
    if j == 1; title('Posterior Distributions'); end
end

% Acceptance rate
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

% K vs n
subplot(2,3,1)
scatter(params_post(:,idx_K), params_post(:,idx_n), 2, ...
    logL_chain(n_burnin+1:end), 'filled', 'MarkerFaceAlpha', 0.1);
hold on
plot(params_map(idx_K), params_map(idx_n), 'rp', 'MarkerSize', 15, ...
    'MarkerFaceColor', 'r');
xlabel('log_{10}(K)'); ylabel('n');
title('K vs n Trade-off');
colorbar; colormap(parula);

% t1 vs U_post
subplot(2,3,2)
scatter(params_post(:,idx_t1) * 1e-6, params_post(:,idx_Upost) * 1e3, 2, ...
    logL_chain(n_burnin+1:end), 'filled', 'MarkerFaceAlpha', 0.1);
hold on
plot(params_map(idx_t1)*1e-6, params_map(idx_Upost)*1e3, 'rp', ...
    'MarkerSize', 15, 'MarkerFaceColor', 'r');
xlabel('t_1 (Ma)'); ylabel('U_{post} (mm/yr)');
title('t_1 vs Most-Recent Rate');

% t1 histogram
subplot(2,3,3)
histogram(params_post(:,idx_t1) * 1e-6, 50, 'Normalization', 'pdf', ...
    'FaceColor', [0.8 0.4 0.4], 'EdgeColor', 'none'); hold on
if is_3phase && ~isempty(idx_t2)
    histogram(params_post(:,idx_t2) * 1e-6, 50, 'Normalization', 'pdf', ...
        'FaceColor', [0.4 0.4 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
if cave_prior.use_informative
    t1_mean_Ma = []; t1_std_Ma = [];
    if isfield(cave_prior, 't1_mean') && isfield(cave_prior, 't1_std')
        t1_mean_Ma = cave_prior.t1_mean * 1e-6;
        t1_std_Ma  = cave_prior.t1_std * 1e-6;
    elseif isfield(cave_prior, 't_capture_mean')
        t1_mean_Ma = cave_prior.t_capture_mean * 1e-6;
        t1_std_Ma  = cave_prior.t_capture_std * 1e-6;
    end
    if ~isempty(t1_mean_Ma)
        x = linspace(0, 8, 200);
        y = normpdf(x, t1_mean_Ma, t1_std_Ma);
        y = y / max(y) * max(ylim) * 0.7;
        plot(x, y, 'g-', 'LineWidth', 2);
    end
end
xlabel('Time (Ma)'); ylabel('PDF');
if is_3phase
    legend('t_1 (older)', 't_2 (younger)', 'Cave Prior');
else
    legend('t_{capture}', 'Cave Prior');
end
title('Transition Times');

% n histogram
subplot(2,3,4)
histogram(params_post(:,idx_n), 50, 'Normalization', 'pdf', ...
    'FaceColor', [0.4 0.7 0.4], 'EdgeColor', 'none'); hold on
xline(1, 'k--', 'n=1', 'LineWidth', 2, 'FontSize', 12);
xline(median(params_post(:,idx_n)), 'r-', 'LineWidth', 2);
if isfield(cave_prior, 'n_mean') && isfield(cave_prior, 'n_std') && cave_prior.n_std > 0
    x = linspace(prior_bounds(idx_n,1), prior_bounds(idx_n,2), 200);
    y = normpdf(x, cave_prior.n_mean, cave_prior.n_std);
    y = y / max(y) * max(ylim) * 0.7;
    plot(x, y, 'g-', 'LineWidth', 1.5);
end
xlabel('n (slope exponent)'); ylabel('PDF');
title('Slope Exponent Posterior');

% K histogram
subplot(2,3,5)
histogram(params_post(:,idx_K), 50, 'Normalization', 'pdf', ...
    'FaceColor', [0.7 0.4 0.7], 'EdgeColor', 'none');
xlabel('log_{10}(K)'); ylabel('PDF');
title('Erodibility Posterior');

% Uplift rates comparison
subplot(2,3,6)
if is_3phase
    boxplot([params_post(:,1)*1e3, params_post(:,2)*1e3, params_post(:,3)*1e3], ...
        'Labels', {'U_{pre}', 'U_{mid}', 'U_{post}'});
    hold on
    plot([1 2 3], params_map([1 2 3])*1e3, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
    ylabel('Incision rate (mm/yr)');
    title('Uplift Rate Posteriors');
else
    scatter(params_post(:,1)*1e3, params_post(:,2)*1e3, 2, ...
        'filled', 'MarkerFaceAlpha', 0.1);
    hold on
    plot(params_map(1)*1e3, params_map(2)*1e3, 'rp', 'MarkerSize', 15, ...
        'MarkerFaceColor', 'r');
    plot([0 0.5], [0 0.5], 'k--');
    xlabel('U_{pre} (mm/yr)'); ylabel('U_{post} (mm/yr)');
    title('Pre vs Post Capture Rates');
end

sgtitle(sprintf('Parameter Posteriors: %s', fileTag), 'FontSize', 14);
saveas(fig2, fullfile(output_dir, ['posteriors_' fileTag '.png']));

%% Figure 3: Model Fit
fig3 = figure('Position', [150, 150, 1600, 700]);

% Identify trunk stream nodes
try
    S_trunk_obj = trunk(klargestconncomps(S, 1));
    [~, trunk_node_idx] = ismember(S_trunk_obj.IXgrid, S.IXgrid);
    trunk_node_idx = trunk_node_idx(trunk_node_idx > 0);
catch
    trunk_node_idx = [];
end

% River profile fit in chi space (avoids K-dependent tau distortion)
subplot(1,3,1)
Schi = [];
if nargin >= 18 && ~isempty(S_DA)
    Schi = compute_chi(S, S_DA, params_map(idx_mn));
end
if ~isempty(Schi)
    plot(Schi, Sz_obs, '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', 3); hold on
    plot(Schi, Z_mod_map, '.', 'Color', [0.2 0.2 0.2], 'MarkerSize', 3);
    xlabel('\chi (m)'); ylabel('Elevation (m)');
else
    plot(Sz_obs, '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', 3); hold on
    plot(Z_mod_map, '.', 'Color', [0.2 0.2 0.2], 'MarkerSize', 3);
    xlabel('Node index'); ylabel('Elevation (m)');
end
legend('Observed (all nodes)', 'MAP model (all nodes)', ...
    'Location', 'northwest');
title('River Profile Fit (\chi space)');
grid on

% Trunk-stream profile in along-river distance
subplot(1,3,2)
if ~isempty(trunk_node_idx)
    x_trunk = S.distance(trunk_node_idx);
    [~, idx_sort] = sort(x_trunk, 'descend');
    x_plot = (max(x_trunk) - x_trunk(idx_sort)) / 1e3;

    plot(x_plot, Sz_obs(trunk_node_idx(idx_sort)), '-', ...
        'Color', [0.3 0.5 0.8], 'LineWidth', 1.5); hold on
    plot(x_plot, Z_mod_map(trunk_node_idx(idx_sort)), '-', ...
        'Color', [0.85 0.2 0.2], 'LineWidth', 1.5);
    xlabel('Distance from outlet (km)'); ylabel('Elevation (m)');
    legend('Observed trunk', 'MAP model trunk', 'Location', 'northwest');
    title('Trunk Profile');
    grid on
else
    text(0.5, 0.5, 'Trunk extraction unavailable', ...
        'HorizontalAlignment', 'center', 'Units', 'normalized');
    title('Trunk Profile');
end

% Cave fit with multi-phase incision history
subplot(1,3,3)
errorbar(cave_ages / 1e6, cave_heights, cave_height_err, 'bo', ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.6 0.9], 'LineWidth', 1.5);
hold on
plot(cave_ages / 1e6, cave_pred_map, 'rs', 'MarkerSize', 10, ...
    'MarkerFaceColor', [0.9 0.3 0.3], 'LineWidth', 1.5);

% Plot model incision history
t_plot = linspace(0, max(cave_ages)*1.1, 200);
if is_3phase
    U_rates_map = [params_map(1), params_map(2), params_map(3)];
    t_trans_map = [params_map(7), params_map(8)];
else
    U_rates_map = [params_map(1), params_map(2)];
    t_trans_map = params_map(6);
end
h_plot = cave_forward_model(t_plot, U_rates_map, t_trans_map);
plot(t_plot / 1e6, h_plot, 'r-', 'LineWidth', 2);

% Mark transition times
if is_3phase
    xline(params_map(7) / 1e6, 'k--', sprintf('t_1=%.1f Ma', params_map(7)/1e6), ...
        'LineWidth', 1.5, 'FontSize', 10);
    xline(params_map(8) / 1e6, 'b--', sprintf('t_2=%.1f Ma', params_map(8)/1e6), ...
        'LineWidth', 1.5, 'FontSize', 10);
else
    xline(params_map(6) / 1e6, 'k--', sprintf('t_{cap}=%.1f Ma', params_map(6)/1e6), ...
        'LineWidth', 1.5, 'FontSize', 11);
end

xlabel('Age (Ma)'); ylabel('Height above river (m)');
legend('Cave data', 'MAP predictions', 'Incision model', 'Location', 'northwest');
title('Cave Data Fit');
set(gca, 'XDir', 'reverse');

sgtitle(sprintf('Model Fit: %s', fileTag), 'FontSize', 14);
saveas(fig3, fullfile(output_dir, ['model_fit_' fileTag '.png']));

fprintf('Figures saved to: %s\n', output_dir);

end

%% Helper functions
function Schi = compute_chi(S, S_DA, mn)
% Compute chi coordinate for the stream network.
% chi = integral from outlet to x of (1/A)^(m/n) dx
try
    Schi = zeros(size(S.distance));
    Six  = S.ix;
    Sixc = S.ixc;
    Sx   = S.distance;
    Sa   = (1 ./ S_DA).^mn;

    for lp = numel(Six):-1:1
        Schi(Six(lp)) = Schi(Sixc(lp)) + ...
            (Sa(Sixc(lp)) + (Sa(Six(lp)) - Sa(Sixc(lp))) / 2) * ...
            abs(Sx(Sixc(lp)) - Sx(Six(lp)));
    end
catch
    Schi = [];
end
end
