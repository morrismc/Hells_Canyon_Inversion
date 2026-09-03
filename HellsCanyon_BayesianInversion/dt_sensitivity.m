%% dt_sensitivity.m
%
% IMPORTANT -- dt sensitivity depends strongly on n, so the earlier result
% that cleared dt = 100,000 does NOT carry over.  That test ran at n ~ 0.6.
% For n > 1 the stream-power equation is shock-forming: knickpoints steepen
% as they migrate instead of spreading out, and the implicit scheme needs a
% finer step to place them correctly.  With the ksn analysis implying
% n ~ 3.7, RE-RUN this before adopting any coarse dt -- and repeat it at the
% high-K edge of the posterior, since knickpoint celerity scales with K.
load('params_HC_capture.mat');
sd = load('hc_stream_data.mat'); sd = sd.stream_data;
S = sd.S; S_DA = double(sd.S_DA(:));

p   = median(params_post,1);
dts = [6250 12500 25000 50000 100000];
Zs  = zeros(numel(S_DA), numel(dts));  tt = zeros(size(dts));

for k = 1:numel(dts)
    tic;
    if numel(p) >= 7, U_pat = hc_uplift_pattern(S, p(7));
    else,             U_pat = [];   % pre-U_grad chain: uniform uplift
    end
    Zs(:,k) = hc_river_forward_model(S, S_DA, [p(1) p(2)], p(6), ...
                                     hc_derive_K(p), p(5)*p(4), p(4), ...
                                     dts(k), U_pat);
    tt(k) = toc;
end

ref = Zs(:,1);   % finest dt = reference truth
fprintf('\n%8s %9s %14s %14s\n','dt (yr)','time (s)','RMS vs ref','max vs ref');
fprintf('%s\n', repmat('-',1,50));
for k = 1:numel(dts)
    fprintf('%8d %9.2f %12.2f m %12.2f m\n', dts(k), tt(k), ...
        sqrt(mean((Zs(:,k)-ref).^2)), max(abs(Zs(:,k)-ref)));
end