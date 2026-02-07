function [lp] = logprior_hc(candidate, prior_bounds, cave_prior)
% LOGPRIOR_HC  Log-prior for Hells Canyon Bayesian inversion.
%
% Combines uniform bounds with optional informative Gaussian priors
% from cave constraints on timing and incision rates.
%
% Inputs:
%   candidate    - [1 x nparams] vector of candidate parameters
%                  [U_pre, U_post, log10_K, n, m_over_n, t_capture]
%   prior_bounds - [nparams x 2] matrix of [lower, upper] bounds
%   cave_prior   - struct with optional informative priors from caves:
%                  .use_informative  - logical, use informative priors?
%                  .t_capture_mean   - Mean capture time (years)
%                  .t_capture_std    - Std of capture time (years)
%                  .U_pre_mean       - Mean pre-capture rate (m/yr)
%                  .U_pre_std        - Std of pre-capture rate (m/yr)
%                  .U_post_mean      - Mean post-capture rate (m/yr)
%                  .U_post_std       - Std of post-capture rate (m/yr)
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

% Step 2: Start with uniform (uninformative) prior within bounds
lp = 0;

% Step 3: Add informative Gaussian priors from cave constraints if requested
if nargin >= 3 && ~isempty(cave_prior) && cave_prior.use_informative

    % Parameter indices (must match main script ordering):
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

    % Physical constraint: U_post > U_pre (capture increases incision)
    if candidate(2) <= candidate(1)
        lp = -Inf;
        return;
    end
end

end
