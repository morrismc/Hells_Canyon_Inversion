function [logL, logL_stream, logL_cave] = hc_loglikelihood(obs_stream, ...
    mod_stream, sigma_stream, obs_cave, mod_cave, sigma_cave)
% HC_LOGLIKELIHOOD  Combined log-likelihood for river profiles + cave data.
%
% Computes Gaussian log-likelihood for two data types with balancing
% weights to prevent the more numerous stream data from dominating.
%
% Inputs:
%   obs_stream   - Observed stream elevations (m)
%   mod_stream   - Modeled stream elevations (m)
%   sigma_stream - Uncertainty on stream elevations (m), scalar or vector
%   obs_cave     - Observed cave heights above river (m)
%   mod_cave     - Modeled cave heights above river (m)
%   sigma_cave   - Uncertainty on cave heights (m), scalar or vector
%
% Outputs:
%   logL         - Total log-likelihood
%   logL_stream  - Stream component
%   logL_cave    - Cave component

% Force column vectors to prevent broadcast/size mismatch
obs_stream   = obs_stream(:);
mod_stream   = mod_stream(:);
sigma_stream = sigma_stream(:);

% Diagnostic: check sizes match
if length(obs_stream) ~= length(mod_stream) || length(obs_stream) ~= length(sigma_stream)
    error('hc_loglikelihood:sizeMismatch', ...
        'Size mismatch: obs_stream=%d, mod_stream=%d, sigma_stream=%d', ...
        length(obs_stream), length(mod_stream), length(sigma_stream));
end

% Stream profile likelihood
resid_stream = (obs_stream - mod_stream) ./ sigma_stream;
logL_stream = -0.5 * sum(resid_stream.^2);

% Cave data likelihood
if ~isempty(obs_cave) && ~isempty(mod_cave)
    resid_cave = (obs_cave - mod_cave) ./ sigma_cave;
    logL_cave = -0.5 * sum(resid_cave.^2);

    % Stream and cave likelihoods are combined with equal weight per data
    % point (Ws = 1).  Cave constraints on timing and rates already enter
    % through informative Gaussian priors (logprior_hc), so there is no
    % need to down-weight the stream data here.  The previous weighting
    % (Ws = n_cave/n_stream ≈ 0.006) effectively silenced the stream
    % profile, producing poor river-profile fits.
    logL = logL_stream + logL_cave;
else
    logL_cave = 0;
    logL = logL_stream;
end

end
