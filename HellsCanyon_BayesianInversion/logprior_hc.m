function [lp] = logprior_hc(candidate, prior_bounds, cave_prior)
% LOGPRIOR_HC  Log-prior for Hells Canyon Bayesian inversion.
%
% Combines uniform bounds with optional informative Gaussian priors
% from cave constraints on timing, incision rates, and the slope
% exponent n.
%
% 6-parameter ordering:
%   [U_pre, U_post, log10_K, n, m_over_n, t_capture]
%
% Inputs:
%   candidate    - [1 x 6] vector of candidate parameters
%   prior_bounds - [6 x 2] matrix of [lower, upper] bounds
%   cave_prior   - struct with optional informative priors:
%                  .use_informative  - logical
%                  .t_capture_mean, .t_capture_std
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

% Step 3: Add informative Gaussian priors if requested
if nargin >= 3 && ~isempty(cave_prior) && cave_prior.use_informative

    % Parameter indices:
    % 1=U_pre, 2=U_post, 3=log10_K, 4=n, 5=m/n, 6=t_capture

    % Informative prior on t_capture from cave burial ages
    if isfield(cave_prior, 't_capture_mean') && cave_prior.t_capture_std > 0
        lp = lp - 0.5 * ((candidate(6) - cave_prior.t_capture_mean) / ...
             cave_prior.t_capture_std)^2;
    end

    % Informative prior on U_pre from relict landscape / upper caves
    if isfield(cave_prior, 'U_pre_mean') && cave_prior.U_pre_std > 0
        lp = lp - 0.5 * ((candidate(1) - cave_prior.U_pre_mean) / ...
             cave_prior.U_pre_std)^2;
    end

    % Informative prior on U_post from lower caves / adjusted landscape
    if isfield(cave_prior, 'U_post_mean') && cave_prior.U_post_std > 0
        lp = lp - 0.5 * ((candidate(2) - cave_prior.U_post_mean) / ...
             cave_prior.U_post_std)^2;
    end

    % Gaussian prior on n (slope exponent) to break the K-n trade-off.
    % Most geomorphological studies find n in [0.5, 2] with n=1 standard.
    if isfield(cave_prior, 'n_mean') && isfield(cave_prior, 'n_std') ...
            && cave_prior.n_std > 0
        lp = lp - 0.5 * ((candidate(4) - cave_prior.n_mean) / ...
             cave_prior.n_std)^2;
    end

    % Physical constraint: U_post > U_pre (capture increases incision)
    if candidate(2) <= candidate(1)
        lp = -Inf;
        return;
    end
end

end
