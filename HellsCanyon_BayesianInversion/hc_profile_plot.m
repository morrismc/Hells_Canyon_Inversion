function out = hc_profile_plot(varargin)
% HC_PROFILE_PLOT  Long profile and chi profile of the drainage network.
%
%   out = hc_profile_plot()                 % theta from mnoptimvar
%   out = hc_profile_plot('mn', 0.6869)     % a specific concavity
%   out = hc_profile_plot('mn', 0.6869, 'mn_compare', 0.547)
%
% Data only -- no forward model, no MCMC.  Shows what the network actually
% looks like, and whether one concavity collapses it.
%
%   Panel 1  long profile: elevation vs distance from the outlet
%   Panel 2  chi profile at 'mn', coloured by log10 drainage area
%   Panel 3  chi profile at 'mn_compare', for a side-by-side (optional)
%
% In panel 2 a network in steady state with uniform K plots as a single
% straight line.  Branches that peel away from the trunk are the ones that
% a single-K model cannot reconcile, and colouring by drainage area shows
% whether the departure is systematic with catchment size (which would
% point at K varying with discharge or lithology) or scattered.
%
% Options:
%   'mn'               concavity for panel 2; default = mnoptimvar optimum
%   'mn_compare'       second concavity for panel 3; default [] (skip)
%   'stream_data_file' default 'hc_stream_data.mat'
%   'a0'               reference area, default 1e6 m^2 (TopoToolbox default)
%
% See also: hc_find_theta, prepare_hc_stream_data

p = inputParser;
addParameter(p, 'mn', [], @(x) isempty(x) || isscalar(x));
addParameter(p, 'mn_compare', [], @(x) isempty(x) || isscalar(x));
addParameter(p, 'stream_data_file', 'hc_stream_data.mat', @ischar);
addParameter(p, 'a0', 1e6, @isscalar);
parse(p, varargin{:});
o = p.Results;

%% --- data -------------------------------------------------------------
if evalin('base','exist(''S'',''var'')') && isa(evalin('base','S'),'STREAMobj')
    S    = evalin('base','S');
    Sz   = double(evalin('base','Sz'));
    S_DA = double(evalin('base','S_DA'));
else
    raw = load(o.stream_data_file);
    if isfield(raw,'stream_data'), sd = raw.stream_data; else, sd = raw; end
    S = sd.S;  Sz = double(sd.Sz(:));  S_DA = double(sd.S_DA(:));
end
Sz = Sz(:); S_DA = S_DA(:);
Sz_norm = Sz - min(Sz);

if isempty(o.mn)
    fprintf('No mn supplied; finding the optimum with mnoptimvar...\n');
    o.mn = mnoptimvar(S, Sz, S_DA, 'plot', false);
    fprintf('  optimal theta = %.4f\n', o.mn);
end

%% --- trunk ------------------------------------------------------------
try
    St = trunk(klargestconncomps(S,1));
    [~, ti] = ismember(St.IXgrid, S.IXgrid);
    ti = ti(ti>0);
catch
    ti = [];
end

d_km  = S.distance / 1e3;
logA  = log10(max(S_DA, 1));

%% --- figure -----------------------------------------------------------
npan = 2 + ~isempty(o.mn_compare);
fig  = figure('Position', [80 80 460*npan 520], 'Color', 'w');

% Panel 1: long profile
subplot(1, npan, 1)
plot(d_km, Sz_norm, '.', 'Color', [0.78 0.78 0.78], 'MarkerSize', 2);
hold on
if ~isempty(ti)
    [dt, is] = sort(d_km(ti), 'ascend');
    plot(dt, Sz_norm(ti(is)), '-', 'Color', [0.15 0.35 0.75], 'LineWidth', 2);
    legend('All channels', 'Trunk', 'Location', 'northwest');
end
xlabel('Distance from outlet (km)');
ylabel('Elevation above outlet (m)');
title('Long profile');
grid on; box on

% Panel 2: chi profile at mn
chi = chitransform(S, S_DA, 'mn', o.mn, 'a0', o.a0);
subplot(1, npan, 2)
scatter(chi, Sz_norm, 3, logA, 'filled');
hold on
if ~isempty(ti)
    [ct, is] = sort(chi(ti), 'ascend');
    plot(ct, Sz_norm(ti(is)), 'k-', 'LineWidth', 2);
end
cb = colorbar; cb.Label.String = 'log_{10} drainage area (m^2)';
colormap(gca, parula);
xlabel('\chi (m)'); ylabel('Elevation above outlet (m)');
title(sprintf('\\chi profile,  \\theta = %.4f', o.mn));
grid on; box on

out = struct('mn', o.mn, 'chi', chi, 'Sz_norm', Sz_norm, ...
             'trunk_idx', ti, 'scatter_mn', local_scatter(chi, Sz_norm));

% Panel 3: comparison concavity
if ~isempty(o.mn_compare)
    chi2 = chitransform(S, S_DA, 'mn', o.mn_compare, 'a0', o.a0);
    subplot(1, npan, 3)
    scatter(chi2, Sz_norm, 3, logA, 'filled');
    hold on
    if ~isempty(ti)
        [ct2, is2] = sort(chi2(ti), 'ascend');
        plot(ct2, Sz_norm(ti(is2)), 'k-', 'LineWidth', 2);
    end
    cb = colorbar; cb.Label.String = 'log_{10} drainage area (m^2)';
    xlabel('\chi (m)'); ylabel('Elevation above outlet (m)');
    title(sprintf('\\chi profile,  \\theta = %.4f', o.mn_compare));
    grid on; box on
    out.chi_compare  = chi2;
    out.scatter_compare = local_scatter(chi2, Sz_norm);
end

sgtitle('Drainage network: long profile and \chi profile', 'FontSize', 13);

%% --- report -----------------------------------------------------------
relief = max(Sz_norm);
fprintf('\nNetwork geometry:\n');
fprintf('  nodes                    : %d\n', numel(Sz_norm));
fprintf('  relief                   : %.0f m\n', relief);
fprintf('  max chi (theta = %.4f)  : %.2f\n', o.mn, max(chi));
fprintf('  chi-z scatter            : %.1f m  (%.1f%% of relief)\n', ...
        out.scatter_mn, 100*out.scatter_mn/max(relief,eps));
if ~isempty(o.mn_compare)
    fprintf('  chi-z scatter at %.4f  : %.1f m  (%.1f%% of relief)\n', ...
            o.mn_compare, out.scatter_compare, ...
            100*out.scatter_compare/max(relief,eps));
end

end

%% ========================================================================
function s = local_scatter(chi, z)
% RMS deviation from the binned-mean chi-z curve: goes to zero only if
% every branch collapses onto one curve.
nb = 100;
edges = linspace(min(chi), max(chi), nb+1);
[~,~,bin] = histcounts(chi, edges);
ok = bin > 0;
mu = accumarray(bin(ok), z(ok), [nb 1], @mean, NaN);
r  = z(ok) - mu(bin(ok));
s  = sqrt(mean(r(isfinite(r)).^2));
end
