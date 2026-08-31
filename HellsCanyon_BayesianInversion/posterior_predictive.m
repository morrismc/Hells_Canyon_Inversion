%% posterior_predictive.m
load('params_HC_capture.mat');                 % params_post
sd = load('hc_stream_data.mat'); sd = sd.stream_data;
S = sd.S;  Sz = double(sd.Sz(:));  S_DA = double(sd.S_DA(:));
Sz_norm = Sz - min(Sz);
dt_forward = 25000;

fwd = @(p) hc_river_forward_model(S, S_DA, [p(1) p(2)], p(6), ...
                                  10^p(3), p(5)*p(4), p(4), dt_forward);

% --- median-parameter model ---
p_med = median(params_post, 1);
Z_med = fwd(p_med);

% --- ensemble, thinned evenly across the chain (not random: keeps spread) ---
n_draw = 100;
idx  = round(linspace(1, size(params_post,1), n_draw));
Zens = zeros(numel(Sz_norm), n_draw);
fprintf('Running %d forward models', n_draw);
for k = 1:n_draw
    Zens(:,k) = fwd(params_post(idx(k),:));
    if mod(k,10)==0, fprintf('.'); end
end
fprintf(' done\n');

% --- trunk nodes, ordered outlet -> head ---
St = trunk(klargestconncomps(S,1));
[~, ti] = ismember(St.IXgrid, S.IXgrid);   ti = ti(ti>0);
[xk, ord] = sort(S.distance(ti), 'ascend');
tn = ti(ord);  xk = xk/1e3;

band = prctile(Zens(tn,:), [2.5 16 50 84 97.5], 2);

figure('Position',[100 100 950 620]);
fill([xk; flipud(xk)], [band(:,1); flipud(band(:,5))], [0.88 0.88 0.96], ...
     'EdgeColor','none'); hold on
fill([xk; flipud(xk)], [band(:,2); flipud(band(:,4))], [0.72 0.72 0.90], ...
     'EdgeColor','none');
plot(xk, Sz_norm(tn), 'b-', 'LineWidth', 2);
plot(xk, Z_med(tn),   'r-', 'LineWidth', 2);
xlabel('Distance from outlet (km)'); ylabel('Elevation (m)');
legend('95% predictive','68% predictive','Observed trunk','Median model', ...
       'Location','northwest');
title('Trunk profile: posterior predictive', 'Interpreter','none');
grid on

% --- quantify: median vs MAP ---
load('mMAP_HC_capture.mat','Z_mod_map');
fprintf('\nRMS residual, all nodes:\n');
fprintf('  median-parameter model : %7.1f m\n', sqrt(mean((Sz_norm-Z_med).^2)));
fprintf('  MAP model              : %7.1f m\n', sqrt(mean((Sz_norm-Z_mod_map).^2)));