function results = hc_linear_inversion_with_caves(DEM, K_range, varargin)
% HC_LINEAR_INVERSION_WITH_CAVES  Run Gallen's linear inversion code at
% multiple K values to bracket the capture timing, with cave age overlay.
%
% This function wraps Sean Gallen's linear_inversion_block_uplift_erodibility
% to run the inversion across a range of K values (from the ksn analysis)
% and overlay the cave-constrained capture timing for comparison.
%
% This is the n=1 (linear) approach. For n!=1, use main_hc_bayes_inversion.m.
%
% Inputs:
%   DEM      - TopoToolbox GRIDobj clipped to watershed
%   K_range  - [K_low, K_mid, K_high] erodibility values to test
%
%   Optional name-value pairs:
%     'tau_inc'       - Response time increment (years, default: 1.5e6)
%     'crita'         - Critical drainage area (m^2, default: 1e6)
%     'mn'            - m/n ratio (default: 0.45)
%     'flowOption'    - 'fill', 'carve', or '' (default: 'fill')
%     'cave_ages'     - Cave burial ages (years, default: [5.5e6, 3.5e6, 1.5e6])
%     'cave_heights'  - Heights above river (m, default: [100, 120, 250])
%     'cave_err'      - Height uncertainties (m, default: [20, 20, 20])
%     't_capture'     - Cave-constrained capture time (years, default: 2.1e6)
%     't_capture_err' - Uncertainty on capture time (years, default: 1.0e6)
%     'output_dir'    - Output directory (default: pwd)
%
% Outputs:
%   results  - struct array (one per K value) with fields:
%     .K, .A, .Umod, .S, .Stau, .tau_steps
%     .capture_in_window  - logical, does the inversion show capture signal
%                           within the cave-constrained time window?
%
% Requires:
%   linear_inversion_block_uplift_erodibility.m (Gallen's code) on path
%   TopoToolbox 2.4+

%% Parse inputs
p = inputParser;
addRequired(p, 'DEM', @(x) isa(x, 'GRIDobj'));
addRequired(p, 'K_range', @(x) numel(x) >= 1);
addParameter(p, 'tau_inc', 1.5e6, @isscalar);
addParameter(p, 'crita', 1e6, @isscalar);
addParameter(p, 'mn', 0.45, @isscalar);
addParameter(p, 'flowOption', 'fill', @ischar);
addParameter(p, 'cave_ages', [5.5e6, 3.5e6, 1.5e6]);
addParameter(p, 'cave_heights', [100, 120, 250]);
addParameter(p, 'cave_err', [20, 20, 20]);
addParameter(p, 't_capture', 2.1e6, @isscalar);
addParameter(p, 't_capture_err', 1.0e6, @isscalar);
addParameter(p, 'output_dir', pwd, @ischar);

parse(p, DEM, K_range, varargin{:});
opts = p.Results;

nK = length(K_range);

%% Check that Gallen's code is on the path
if ~exist('linear_inversion_block_uplift_erodibility', 'file')
    % Try to find it
    block_uplift_path = fullfile(fileparts(mfilename('fullpath')), '..', ...
        'Block_Uplift_Linear_Inversion_Models', 'dimensional_block_uplift');
    if exist(block_uplift_path, 'dir')
        addpath(block_uplift_path);
        fprintf('Added Gallen code to path: %s\n', block_uplift_path);
    else
        error(['Cannot find linear_inversion_block_uplift_erodibility.m. ', ...
               'Add Block_Uplift_Linear_Inversion_Models to your MATLAB path.']);
    end
end

%% Run inversions for each K
results = struct();
colors = lines(nK);

for k = 1:nK
    K = K_range(k);
    fprintf('\n===== Running inversion with K = %.2e =====\n', K);

    % Adjust tau_inc based on K to get ~10 time steps
    % tau_max ~ chi_max / K, so scale tau_inc accordingly
    tau_inc_k = opts.tau_inc;

    try
        % Close any figures opened by Gallen's code
        close(gcf);

        [A, Umod, S, Stau, tau_steps] = linear_inversion_block_uplift_erodibility(...
            DEM, K, tau_inc_k, opts.crita, opts.mn, opts.flowOption);

        % Store results
        results(k).K = K;
        results(k).A = A;
        results(k).Umod = Umod;
        results(k).S = S;
        results(k).Stau = Stau;
        results(k).tau_steps = tau_steps;
        results(k).tau_max = max(tau_steps);

        % Check if capture signal falls in cave window
        t_low = opts.t_capture - opts.t_capture_err;
        t_high = opts.t_capture + opts.t_capture_err;

        % Find the peak uplift rate and its timing
        [U_peak, peak_idx] = max(Umod);
        t_peak = tau_steps(peak_idx);

        results(k).U_peak = U_peak;
        results(k).t_peak = t_peak;
        results(k).capture_in_window = (t_peak >= t_low && t_peak <= t_high);

        fprintf('  Max tau: %.1f Ma\n', max(tau_steps)/1e6);
        fprintf('  Peak uplift: %.2f mm/yr at tau = %.1f Ma\n', ...
            U_peak*1000, t_peak/1e6);
        fprintf('  Cave window: %.1f - %.1f Ma  --  %s\n', ...
            t_low/1e6, t_high/1e6, ...
            ternary(results(k).capture_in_window, 'CONSISTENT', 'OUTSIDE'));

    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        results(k).K = K;
        results(k).error = ME.message;
        continue;
    end
end

%% Create comparison figure
figure('Position', [100, 100, 1600, 900]);

for k = 1:nK
    if ~isfield(results(k), 'Umod') || isempty(results(k).Umod)
        continue;
    end

    % Uplift history
    subplot(2, nK, k)
    stairs(results(k).tau_steps / 1e6, results(k).Umod * 1000, ...
        'Color', colors(k,:), 'LineWidth', 2); hold on
    xlabel('\tau (Ma)'); ylabel('U (mm/yr)');
    title(sprintf('K = %.1e', K_range(k)));

    % Cave constraint overlay
    t_low  = (opts.t_capture - opts.t_capture_err) / 1e6;
    t_high = (opts.t_capture + opts.t_capture_err) / 1e6;
    ylims = ylim;
    fill([t_low, t_high, t_high, t_low], ...
         [ylims(1), ylims(1), ylims(2), ylims(2)], ...
         'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    xline(opts.t_capture / 1e6, 'r--', 'LineWidth', 1.5);

    % Cave incision rates
    yline(0.01, 'b--', 'U_{pre}', 'LineWidth', 1);
    yline(0.125, 'r--', 'U_{post}', 'LineWidth', 1);

    legend('Inversion', 'Cave window', 't_{cap} (cave)', ...
        'Location', 'northeast');

    % Cumulative incision
    subplot(2, nK, nK + k)
    cum_incision = cumtrapz(results(k).tau_steps, results(k).Umod);
    stairs(results(k).tau_steps / 1e6, cum_incision, ...
        'Color', colors(k,:), 'LineWidth', 2); hold on

    % Cave heights
    errorbar(opts.cave_ages / 1e6, opts.cave_heights, opts.cave_err, ...
        'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'LineWidth', 1.5);

    xlabel('\tau (Ma)'); ylabel('Cumulative Incision (m)');
    title(sprintf('Cumulative: K = %.1e', K_range(k)));
    legend('Model', 'Cave data', 'Location', 'northwest');
end

sgtitle('Linear Inversion Sensitivity to K (n=1)', 'FontSize', 14);
saveas(gcf, fullfile(opts.output_dir, 'linear_inversion_K_sensitivity.png'));

%% Summary figure: all K values overlaid
figure('Position', [200, 200, 800, 600]);

subplot(2,1,1)
for k = 1:nK
    if ~isfield(results(k), 'Umod') || isempty(results(k).Umod); continue; end
    stairs(results(k).tau_steps / 1e6, results(k).Umod * 1000, ...
        'Color', colors(k,:), 'LineWidth', 2); hold on
end
% Cave window
ylims = ylim;
fill([(opts.t_capture - opts.t_capture_err)/1e6, ...
      (opts.t_capture + opts.t_capture_err)/1e6, ...
      (opts.t_capture + opts.t_capture_err)/1e6, ...
      (opts.t_capture - opts.t_capture_err)/1e6], ...
     [ylims(1), ylims(1), ylims(2), ylims(2)], ...
     'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
xlabel('\tau (Ma)'); ylabel('U (mm/yr)');
legend_strs = arrayfun(@(k) sprintf('K=%.1e', K_range(k)), 1:nK, 'UniformOutput', false);
legend_strs{end+1} = 'Cave window';
legend(legend_strs, 'Location', 'northeast');
title('Uplift Rate History: Sensitivity to K');

subplot(2,1,2)
for k = 1:nK
    if ~isfield(results(k), 'Umod') || isempty(results(k).Umod); continue; end
    cum_inc = cumtrapz(results(k).tau_steps, results(k).Umod);
    stairs(results(k).tau_steps / 1e6, cum_inc, ...
        'Color', colors(k,:), 'LineWidth', 2); hold on
end
errorbar(opts.cave_ages / 1e6, opts.cave_heights, opts.cave_err, ...
    'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'LineWidth', 1.5);
xlabel('\tau (Ma)'); ylabel('Cumulative Incision (m)');
legend([legend_strs(1:nK), {'Cave data'}], 'Location', 'northwest');
title('Cumulative Incision vs Cave Data');

saveas(gcf, fullfile(opts.output_dir, 'linear_inversion_summary.png'));

fprintf('\nFigures saved to: %s\n', opts.output_dir);

end

function s = ternary(cond, a, b)
if cond; s = a; else; s = b; end
end
