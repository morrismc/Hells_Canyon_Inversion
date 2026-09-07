function out = hc_find_theta(varargin)
% HC_FIND_THETA  Independent estimate of the concavity theta (= m/n), and a
% test of whether a single theta collapses the tributaries at all.
%
%   out = hc_find_theta()
%   out = hc_find_theta('stream_data_file', '...', 'mn_posterior', 0.547)
%
% Uses TopoToolbox's mnoptimvar, which finds the mn-ratio that minimizes
% the variability of elevation within chi-distance bins -- i.e. the theta
% that best COLLAPSES every channel onto one chi-z curve.
%
% Why this matters here.  Under uniform K and a common steady state, all
% branches must fall on a single chi-z line at the correct theta.  In the
% Hells Canyon fit they visibly do not: the tributaries fan apart and each
% bends over at high chi, and 1e6 well-mixed samples cannot fix it, which
% points at the model rather than the sampler.  This function separates the
% two possibilities:
%
%   (a) the optimal theta collapses the network and is close to the
%       posterior m/n  -> concavity is fine, and the residual is a real
%       transient / relict signal;
%   (b) NO theta collapses the network (residual scatter stays high at the
%       optimum) -> uniform K across this network is refuted, and the
%       tributaries need their own K rather than a smooth spatial ramp.
%
% mnoptimvar accepts node-attribute lists directly -- mnoptimvar(S,z,a) --
% so hc_stream_data.mat is sufficient and no GRIDobj is needed.  Its
% default reference area a0 = 1e6 m^2 matches S_DA already being in m^2.
%
% See also: prepare_hc_stream_data, plot_hc_results

p = inputParser;
addParameter(p, 'stream_data_file', 'hc_stream_data.mat', @ischar);
addParameter(p, 'mn_posterior', [], @(x) isempty(x) || isscalar(x));
addParameter(p, 'varfun', @iqr);      % mnoptimvar default
addParameter(p, 'mn0', 0.5);
addParameter(p, 'plot', true);
parse(p, varargin{:});
o = p.Results;

if exist('mnoptimvar', 'file') ~= 2
    error('hc_find_theta:noTTB', ...
        ['mnoptimvar not found. Add TopoToolbox to the path:\n' ...
         '  addpath(genpath(''C:\\path\\to\\topotoolbox''))']);
end

%% --- data -------------------------------------------------------------
if evalin('base', 'exist(''S'',''var'')') && ...
        isa(evalin('base','S'), 'STREAMobj')
    S    = evalin('base', 'S');
    Sz   = double(evalin('base', 'Sz'));
    S_DA = double(evalin('base', 'S_DA'));
    fprintf('Using S, Sz, S_DA from the base workspace.\n');
else
    raw = load(o.stream_data_file);
    if isfield(raw, 'stream_data'), sd = raw.stream_data; else, sd = raw; end
    S    = sd.S;
    Sz   = double(sd.Sz(:));
    S_DA = double(sd.S_DA(:));
end
Sz = Sz(:); S_DA = S_DA(:);

%% --- optimize ---------------------------------------------------------
fprintf('\nOptimizing theta with mnoptimvar (minimum-variance)...\n');
[mn_opt, cm, zm, zsd] = mnoptimvar(S, Sz, S_DA, ...
    'varfun', o.varfun, 'mn0', o.mn0, 'plot', false);

fprintf('  Optimal theta (m/n) = %.4f\n', mn_opt);
if ~isempty(o.mn_posterior)
    fprintf('  Inversion posterior m/n = %.4f  (difference %+.4f)\n', ...
            o.mn_posterior, o.mn_posterior - mn_opt);
end

%% --- collapse quality -------------------------------------------------
% Residual scatter about the binned mean chi-z curve, at the OPTIMAL theta.
% If this stays large at the optimum, no single theta collapses the
% network and uniform K is refuted.
chi_opt = chitransform(S, S_DA, 'mn', mn_opt, 'a0', 1e6);
scatter_opt = local_scatter(chi_opt, Sz);

fprintf('\n  Collapse quality at the optimum:\n');
fprintf('    RMS scatter about binned chi-z curve : %.1f m\n', scatter_opt);
fprintf('    Total relief                          : %.1f m\n', ...
        max(Sz) - min(Sz));
fprintf('    Scatter as %% of relief                : %.1f %%\n', ...
        100 * scatter_opt / max(max(Sz) - min(Sz), eps));

if ~isempty(o.mn_posterior)
    chi_pos = chitransform(S, S_DA, 'mn', o.mn_posterior, 'a0', 1e6);
    scatter_pos = local_scatter(chi_pos, Sz);
    fprintf('    (same at posterior m/n = %.3f)         : %.1f m\n', ...
            o.mn_posterior, scatter_pos);
else
    chi_pos = []; scatter_pos = NaN;
end

fprintf(['\n  READ IT AS: scatter well under ~5%% of relief means a single\n' ...
         '  theta collapses the network, so the concavity is right and the\n' ...
         '  profile misfit is a real transient signal. Scatter still >10%%\n' ...
         '  of relief at the OPTIMUM means no theta collapses it, and a\n' ...
         '  single K for the whole network is refuted.\n']);

%% --- figure -----------------------------------------------------------
if o.plot
    fig = figure('Position', [100 100 1250 480], 'Color', 'w');

    subplot(1,2,1)
    plot(chi_opt, Sz, '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', 2);
    hold on
    plot(cm(cm>0), zm(cm>0), 'r-', 'LineWidth', 2);
    xlabel('\chi (m)'); ylabel('Elevation (m)');
    title(sprintf('Optimal \\theta = %.3f  (scatter %.0f m)', ...
                  mn_opt, scatter_opt));
    legend('All nodes', 'Binned mean', 'Location', 'northwest');
    grid on

    subplot(1,2,2)
    if ~isempty(chi_pos)
        plot(chi_pos, Sz, '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', 2);
        xlabel('\chi (m)'); ylabel('Elevation (m)');
        title(sprintf('Inversion m/n = %.3f  (scatter %.0f m)', ...
                      o.mn_posterior, scatter_pos));
        grid on
    else
        errorbar(cm(cm>0), zm(cm>0), zsd(cm>0), 'k-', 'LineWidth', 1.2);
        xlabel('\chi (m)'); ylabel('Elevation (m)');
        title('Binned chi-z with 1\sigma spread');
        grid on
    end

    sgtitle('Concavity check: does one \theta collapse the network?', ...
            'FontSize', 13);
end

%% --- package ----------------------------------------------------------
out = struct('mn_opt', mn_opt, 'chi_bins', cm, 'z_bins', zm, ...
             'z_sd', zsd, 'scatter_opt', scatter_opt, ...
             'scatter_posterior', scatter_pos, ...
             'mn_posterior', o.mn_posterior, ...
             'relief', max(Sz) - min(Sz));

end

%% ========================================================================
function s = local_scatter(chi, z)
% RMS deviation of elevation from the binned-mean chi-z curve.  This is the
% quantity that goes to zero only if every branch collapses onto one curve.
nb = 100;
edges = linspace(min(chi), max(chi), nb+1);
[~, ~, bin] = histcounts(chi, edges);
ok = bin > 0;
mu = accumarray(bin(ok), z(ok), [nb 1], @mean, NaN);
r  = z(ok) - mu(bin(ok));
s  = sqrt(mean(r(isfinite(r)).^2));
end
