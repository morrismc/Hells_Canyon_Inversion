%% ess_check.m -- how many INDEPENDENT samples do you actually have?
load('params_HC_capture.mat');       % params_post, param_names

X = params_post;
[N, P] = size(X);
maxlag = min(20000, floor(N/5));     % ACF is unreliable beyond ~N/5
rate_per_sec = 0.4;                  % your measured throughput
target_ESS   = 200;

IACT = zeros(1,P);  ESS = zeros(1,P);

for p = 1:P
    x = X(:,p) - mean(X(:,p));
    if var(x,1) == 0, IACT(p)=Inf; ESS(p)=0; continue; end

    % autocorrelation via zero-padded FFT
    nfft = 2^nextpow2(2*N);
    f    = fft(x, nfft);
    acf  = real(ifft(f .* conj(f)));
    acf  = acf(1:maxlag+2);
    acf  = acf / acf(1);             % acf(1) = rho(0) = 1

    % Geyer initial-positive-sequence estimator
    s = 0; m = 0;
    while true
        i1 = 2*m + 1;  i2 = 2*m + 2;
        if i2 > numel(acf), break; end
        G = acf(i1) + acf(i2);
        if G <= 0, break; end
        s = s + G;  m = m + 1;
    end
    IACT(p) = max(2*s - 1, 1);       % tau = 2*sum(Gamma) - 1
    ESS(p)  = N / IACT(p);
end

fprintf('\nPost-burn-in samples: %d\n\n', N);
fprintf('%-20s %10s %10s %14s %12s\n', ...
        'Parameter','IACT','ESS','need (iters)','extra hours');
fprintf('%s\n', repmat('-',1,70));
for p = 1:P
    need  = target_ESS * IACT(p);
    hours = max(need - N, 0) / rate_per_sec / 3600;
    fprintf('%-20s %10.0f %10.1f %14.0f %12.1f\n', ...
            param_names{p}, IACT(p), ESS(p), need, hours);
end
fprintf('\nWorst parameter drives the run length.\n');