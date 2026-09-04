function summary = hc_multibasin_summary(res, output_dir)
% HC_MULTIBASIN_SUMMARY  QA table and comparison figure across basins.
%
%   summary = hc_multibasin_summary(res, output_dir)
%
% res is a cell array of run_hc_inversion outputs.  Produces:
%   (1) a QA table you scan INSTEAD of eyeballing 4 figures per basin, and
%   (2) the cross-basin t_capture comparison -- the actual scientific
%       output of a multi-basin run.
%
% The QA gates below are the failure modes this project has actually hit,
% not hypothetical ones: parameters railed against prior bounds (three
% separate production runs lost to it), acceptance far from the efficient
% range, non-finite likelihoods, and forward-model failures.
%
% See also: run_all_basins, run_hc_inversion

if nargin < 2 || isempty(output_dir), output_dir = pwd; end
nb = numel(res);
if nb == 0, summary = struct([]); return; end

pn = {'U_pre','U_post','ksn_ref','n','m_over_n','t_capture','U_grad'};

summary = struct([]);
for b = 1:nb
    r = res{b};  d = r.diag;
    s = struct();
    s.name        = r.name;
    s.n_stream    = d.n_stream;
    s.has_caves   = d.has_caves;
    s.accept_post = d.accept_post;
    s.rms_map     = d.rms_map;
    s.iter_per_s  = d.iter_per_s;
    s.n_ffail     = d.n_ffail;
    s.n_nonfinite = d.n_nonfinite;

    for j = 1:min(numel(pn), size(r.params_post,2))
        s.([pn{j} '_med'])  = median(r.params_post(:,j));
        s.([pn{j} '_ci95']) = d.ci95(j,:);
        s.([pn{j} '_railed']) = d.railed(j);
    end
    s.railed_names = pn(d.railed(1:min(numel(pn), numel(d.railed))));

    % Overall verdict
    flags = {};
    if any(d.railed),                flags{end+1} = 'RAILED';    end %#ok<AGROW>
    if d.accept_post < 0.15 || d.accept_post > 0.50
                                     flags{end+1} = 'ACCEPT';    end %#ok<AGROW>
    if d.n_nonfinite > 0.01*numel(r.accepted)
                                     flags{end+1} = 'NONFINITE'; end %#ok<AGROW>
    if d.n_ffail > 0.01*numel(r.accepted)
                                     flags{end+1} = 'FWDFAIL';   end %#ok<AGROW>
    s.flags = flags;
    s.ok    = isempty(flags);
    summary(b).data = s; %#ok<AGROW>
end

%% ------------------------------------------------------------- table
fprintf('\n================ MULTI-BASIN QA SUMMARY ================\n');
fprintf('%-18s %7s %7s %8s %9s  %s\n', ...
        'Basin','accept','RMS(m)','t_cap','n','flags');
fprintf('%s\n', repmat('-', 1, 72));
for b = 1:nb
    s = summary(b).data;
    fl = 'ok';
    if ~isempty(s.flags), fl = strjoin(s.flags, ','); end
    fprintf('%-18s %6.1f%% %7.1f %8.2f %9.2f  %s\n', ...
            s.name, 100*s.accept_post, s.rms_map, ...
            s.t_capture_med/1e6, s.n_med, fl);
end
fprintf('%s\n', repmat('-', 1, 72));

for b = 1:nb
    s = summary(b).data;
    if ~isempty(s.railed_names)
        fprintf('  %s: RAILED against prior bounds -> %s\n', ...
                s.name, strjoin(s.railed_names, ', '));
    end
end
fprintf(['\nA railed parameter is set by the prior, not the data, and drags\n' ...
         'its correlates with it. Widen the bound and re-run that basin\n' ...
         'before interpreting ANY of its parameters.\n']);

%% ------------------------------------------------------- comparison
try
    fig = figure('Position', [100 100 1200 520]);

    % t_capture across basins -- the money plot
    subplot(1,2,1); hold on
    for b = 1:nb
        s = summary(b).data;
        ci = s.t_capture_ci95 / 1e6;
        plot(ci, [b b], '-', 'Color', [0.4 0.4 0.7], 'LineWidth', 2);
        plot(s.t_capture_med/1e6, b, 'o', 'MarkerSize', 9, ...
             'MarkerFaceColor', [0.85 0.33 0.1], 'MarkerEdgeColor', 'k');
    end
    set(gca, 'YTick', 1:nb, 'YTickLabel', {summary.data.name} , ...
             'YLim', [0.5 nb+0.5], 'TickLabelInterpreter', 'none');
    xlabel('t_{capture} (Ma)');
    title({'Capture timing by basin', ...
           'Agreement across independent basins = regional signal'});
    grid on

    % n across basins
    subplot(1,2,2); hold on
    for b = 1:nb
        s = summary(b).data;
        plot(s.n_ci95, [b b], '-', 'Color', [0.4 0.7 0.4], 'LineWidth', 2);
        plot(s.n_med, b, 'o', 'MarkerSize', 9, ...
             'MarkerFaceColor', [0.2 0.6 0.2], 'MarkerEdgeColor', 'k');
    end
    xline(1, 'k--', 'n=1');
    set(gca, 'YTick', 1:nb, 'YTickLabel', {summary.data.name}, ...
             'YLim', [0.5 nb+0.5], 'TickLabelInterpreter', 'none');
    xlabel('n (slope exponent)');
    title('Slope exponent by basin');
    grid on

    sgtitle('Multi-basin comparison', 'Interpreter', 'none');
    saveas(fig, fullfile(output_dir, 'multibasin_comparison.png'));
    fprintf('\nComparison figure: %s\n', ...
            fullfile(output_dir, 'multibasin_comparison.png'));
catch ME
    fprintf('(comparison figure skipped: %s)\n', ME.message);
end

end
