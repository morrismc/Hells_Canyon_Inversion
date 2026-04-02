function [logL, logL_stream, logL_cave] = hc_loglikelihood(obs_stream, ...
    mod_stream, sigma_stream, obs_cave, mod_cave, sigma_cave)
% HC_LOGLIKELIHOOD  Combined log-likelihood for river profiles + cave data.
%
% Balances the two datasets using the error-inflation approach of Gallen &
% Fernandez-Blanco (2021): stream-node uncertainties are inflated by
% sqrt(n_stream / n_cave) so that the *per-dataset* contribution to the
% likelihood is comparable, despite the large difference in data counts.
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
%   logL_stream  - Stream component (after inflation)
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

% Cave data likelihood
if ~isempty(obs_cave) && ~isempty(mod_cave)
    n_stream = length(obs_stream);
    n_cave   = length(obs_cave);

    % Gallen-style error inflation: multiply stream sigma by
    % sqrt(n_stream / n_cave).  This is equivalent to weighting the stream
    % log-likelihood by n_cave / n_stream, but it operates on the errors
    % directly (matching the published implementation) so the combined
    % likelihood can be summed without an explicit Ws factor.
    sigma_stream_eff = sigma_stream * sqrt(n_stream / n_cave);

    resid_stream = (obs_stream - mod_stream) ./ sigma_stream_eff;
    logL_stream  = -0.5 * sum(resid_stream.^2);

    resid_cave = (obs_cave - mod_cave) ./ sigma_cave;
    logL_cave  = -0.5 * sum(resid_cave.^2);

    logL = logL_stream + logL_cave;
else
    % No cave data — use stream likelihood unweighted
    resid_stream = (obs_stream - mod_stream) ./ sigma_stream;
    logL_stream  = -0.5 * sum(resid_stream.^2);
    logL_cave    = 0;
    logL         = logL_stream;
end

end
