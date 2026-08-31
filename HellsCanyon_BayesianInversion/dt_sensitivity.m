%% dt_sensitivity.m
load('params_HC_capture.mat');
sd = load('hc_stream_data.mat'); sd = sd.stream_data;
S = sd.S; S_DA = double(sd.S_DA(:));

p   = median(params_post,1);
dts = [6250 12500 25000 50000 100000];
Zs  = zeros(numel(S_DA), numel(dts));  tt = zeros(size(dts));

for k = 1:numel(dts)
    tic;
    Zs(:,k) = hc_river_forward_model(S, S_DA, [p(1) p(2)], p(6), ...
                                     10^p(3), p(5)*p(4), p(4), dts(k));
    tt(k) = toc;
end

ref = Zs(:,1);   % finest dt = reference truth
fprintf('\n%8s %9s %14s %14s\n','dt (yr)','time (s)','RMS vs ref','max vs ref');
fprintf('%s\n', repmat('-',1,50));
for k = 1:numel(dts)
    fprintf('%8d %9.2f %12.2f m %12.2f m\n', dts(k), tt(k), ...
        sqrt(mean((Zs(:,k)-ref).^2)), max(abs(Zs(:,k)-ref)));
end