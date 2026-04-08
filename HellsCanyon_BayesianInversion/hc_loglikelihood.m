function [logL, logL_stream, logL_cave] = hc_loglikelihood(obs_stream, ...
    mod_stream, sigma_stream, obs_cave, mod_cave, sigma_cave)
% HC_LOGLIKELIHOOD  Combined log-likelihood for river profiles + cave data.
%
% Uses a moderate balancing approach:  the Gallen & Fernandez-Blanco (2021)
% sqrt-scaling is capped so that with very unbalanced datasets (e.g. 5000
% stream nodes vs. 3 caves) the stream data retains enough weight to
% constrain K, n, and m/n, while cave data still contributes meaningfully
% through both the likelihood and informative priors.
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

    % Moderate balancing: scale stream logL so that the *effective* number
    % of stream data points equals a user-tunable target.  This preserves
    % enough stream weight to constrain K/n/m-n while letting cave data
    % contribute.  Default N_eff = 50 (i.e., treat the stream profile as
    % ~50 independent constraints rather than 5000+ correlated nodes).
    n_stream = length(obs_stream);
    N_eff_stream = 50;  % effective independent stream constraints
    Ws = N_eff_stream / n_stream;

    logL = Ws * logL_stream + logL_cave;
else
    logL_cave = 0;
    logL = logL_stream;
end

end
