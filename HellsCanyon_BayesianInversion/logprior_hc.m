function [lp] = logprior_hc(candidate, prior_bounds, cave_prior)
% LOGPRIOR_HC  Log-prior for Hells Canyon Bayesian inversion.
%
% Combines uniform bounds with optional informative Gaussian priors
% from cave constraints on timing, incision rates, and the slope
% exponent n.
%
% 8-parameter ordering:
%   [U_pre, U_mid, U_post, log10_K, n, m_over_n, t1, t2]
%
% Inputs:
%   candidate    - [1 x 8] vector of candidate parameters
%   prior_bounds - [8 x 2] matrix of [lower, upper] bounds
%   cave_prior   - struct with optional informative priors from caves:
%                  .use_informative  - logical
%                  .t1_mean, .t1_std - Gaussian prior on t1
%                  .U_pre_mean, .U_pre_std
%                  .U_post_mean, .U_post_std
%                  .n_mean, .n_std   - Gaussian prior on slope exponent
%
% Outputs:
%   lp           - Log-prior probability

% Step 1: Check uniform bounds (hard walls)
for i = 1:length(candidate)
    if candidate(i) < prior_bounds(i,1) || candidate(i) > prior_bounds(i,2)
        lp = -Inf;
        return;
    end
end

% Step 2: Start with uniform prior within bounds
lp = 0;

% Step 3: Physical ordering constraints
% t1 must be strictly greater than t2 (older transition before younger)
if candidate(7) <= candidate(8)
    lp = -Inf;
    return;
end

% Step 4: Add informative Gaussian priors if requested
if nargin >= 3 && ~isempty(cave_prior) && cave_prior.use_informative

    % Prior on t1 (main capture event from cave burial ages)
    if isfield(cave_prior, 't1_mean') && isfield(cave_prior, 't1_std') ...
            && cave_prior.t1_std > 0
        lp = lp - 0.5 * ((candidate(7) - cave_prior.t1_mean) / ...
             cave_prior.t1_std)^2;
    end

    % Backwards-compatible: also check for old field name t_capture_mean
    if isfield(cave_prior, 't_capture_mean') && ~isfield(cave_prior, 't1_mean')
        if cave_prior.t_capture_std > 0
            lp = lp - 0.5 * ((candidate(7) - cave_prior.t_capture_mean) / ...
                 cave_prior.t_capture_std)^2;
        end
    end

    % Prior on U_pre (background rate from relict landscape / upper caves)
    if isfield(cave_prior, 'U_pre_mean') && cave_prior.U_pre_std > 0
        lp = lp - 0.5 * ((candidate(1) - cave_prior.U_pre_mean) / ...
             cave_prior.U_pre_std)^2;
    end

    % Prior on U_post (most recent incision rate)
    if isfield(cave_prior, 'U_post_mean') && cave_prior.U_post_std > 0
        lp = lp - 0.5 * ((candidate(3) - cave_prior.U_post_mean) / ...
             cave_prior.U_post_std)^2;
    end

    % Gaussian prior on n (slope exponent) to break the K-n trade-off.
    % Most geomorphological studies find n in [0.5, 2] with n=1 standard.
    if isfield(cave_prior, 'n_mean') && isfield(cave_prior, 'n_std') ...
            && cave_prior.n_std > 0
        lp = lp - 0.5 * ((candidate(5) - cave_prior.n_mean) / ...
             cave_prior.n_std)^2;
    end

    % Physical constraint: post-capture rates should be faster than
    % pre-capture (drainage capture accelerates incision).
    % U_mid > U_pre  AND  U_post > U_pre
    if candidate(2) <= candidate(1) || candidate(3) <= candidate(1)
        lp = -Inf;
        return;
    end
end

end
