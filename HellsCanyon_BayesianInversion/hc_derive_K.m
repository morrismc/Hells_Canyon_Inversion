function K = hc_derive_K(p)
% HC_DERIVE_K  Erodibility implied by the sampled parameters.
%
%   K = U_pre / ksn_ref^n
%
% with the 6-parameter row layout used throughout this project:
%   p = [U_pre, U_post, ksn_ref, n, m_over_n, t_capture]
%
% Accepts a single row or an [N x 6] matrix of posterior samples and
% returns a scalar or [N x 1] column of K values.  Evaluating K per
% sample (rather than from summary statistics) correctly propagates the
% U_pre / ksn_ref / n covariance into the K posterior.
%
% Why K is derived rather than sampled: at steady state the relict
% channel has ksn = (U_pre/K)^(1/n), so the K consistent with the data
% shifts by ~6 orders of magnitude as n varies.  Sampling log10(K) and n
% as independent coordinates proposes across that narrow curved ridge and
% mixes very poorly (IACT ~5400 in the previous run).  ksn_ref is
% directly measurable and only weakly correlated with n.
%
% GUARD: chains produced before the reparameterization stored log10(K) in
% column 3, which is NEGATIVE.  Raising a negative base to a fractional
% power returns a COMPLEX number in MATLAB without raising an error, so
% an old file would silently poison the forward model.  Caught here.

if size(p, 2) < 4
    error('hc_derive_K:badSize', ...
          'Expected at least 6 columns [U_pre U_post ksn_ref n m/n t_cap], got %d.', ...
          size(p, 2));
end

if any(p(:,3) <= 0)
    error('hc_derive_K:oldParameterization', ...
         ['Column 3 contains non-positive values, so this looks like a chain ' ...
          'from the OLD log10(K) parameterization.\nColumn 3 must now be ' ...
          'ksn_ref (positive, order 100). Re-run the inversion, or check out ' ...
          'the pre-reparameterization code from git history to analyse that file.']);
end

K = p(:,1) ./ p(:,3).^p(:,4);

end
