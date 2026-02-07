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

% Stream profile likelihood
resid_stream = (obs_stream - mod_stream) ./ sigma_stream;
logL_stream = -0.5 * sum(resid_stream.^2);

% Cave data likelihood
if ~isempty(obs_cave) && ~isempty(mod_cave)
    resid_cave = (obs_cave - mod_cave) ./ sigma_cave;
    logL_cave = -0.5 * sum(resid_cave.^2);

    % Balance weights: scale stream likelihood so cave data has equal influence
    n_stream = length(obs_stream);
    n_cave   = length(obs_cave);
    Ws = n_cave / n_stream;

    logL = Ws * logL_stream + logL_cave;
else
    logL_cave = 0;
    logL = logL_stream;
end

end
