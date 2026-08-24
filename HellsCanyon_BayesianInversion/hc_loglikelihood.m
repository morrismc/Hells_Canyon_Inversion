function [logL, logL_stream, logL_cave] = hc_loglikelihood(obs_stream, ...
    mod_stream, sigma_stream, obs_cave, mod_cave, sigma_cave)
% HC_LOGLIKELIHOOD  Combined log-likelihood for river profiles + cave data.
%
% Dataset weighting.  Gallen & Fernandez-Blanco (2021) balance their two
% datasets by inflating the stream errors by sqrt(n_stream/n_terrace),
% which is equivalent to down-weighting the stream log-likelihood so its
% EFFECTIVE number of data points equals the terrace count (equal total
% weight per dataset).  Applied literally here, that would give the
% ~5000-node stream network an effective weight of only n_cave = 3 points,
% leaving K, n, and m/n essentially unconstrained by the profile.  We
% therefore use the same mechanism (a multiplicative weight Ws on the
% stream log-likelihood, identical algebra to Gallen's error inflation)
% but set the effective stream count to N_eff_stream instead of n_cave —
% a deliberate, documented departure motivated by the extreme 5000:3
% imbalance and by spatial autocorrelation of profile residuals.
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

    % Scale stream logL so that the *effective* number of stream data
    % points equals a tunable target (see header note on how this relates
    % to Gallen's equal-weight error inflation; Ws here plays the same
    % role as his 1/Ws^2 error scaling).
    %
    % For a river profile with strong spatial autocorrelation (node
    % spacing ~30 m, autocorrelation length of hundreds of meters), a
    % reasonable effective count is n_stream / 30 to n_stream / 60.  For
    % ~5000 nodes that is ~80-170 independent constraints.  We use
    % N_eff_stream = 100 which pushes a bit harder on the stream-power
    % parameters than the previous value of 50 without over-weighting
    % correlated noise.
    n_stream = length(obs_stream);
    N_eff_stream = 100;  % effective independent stream constraints
    Ws = N_eff_stream / n_stream;

    logL = Ws * logL_stream + logL_cave;
else
    logL_cave = 0;
    logL = logL_stream;
end

end
