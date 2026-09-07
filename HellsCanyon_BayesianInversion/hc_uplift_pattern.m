function U_pat = hc_uplift_pattern(S, U_grad)
% HC_UPLIFT_PATTERN  Normalized spatial uplift/incision pattern.
%
%   U(x) = U_outlet * (1 + U_grad * x/x_max)
%
% where x is along-channel distance from the outlet.  The pattern is
% normalized so that U_pat = 1 AT THE OUTLET.  That choice matters:
% it means the sampled U_pre and U_post remain the rates at the outlet,
% which is where the caves are, so cave_forward_model needs no change and
% the rate parameters stay directly comparable to the cave data.
%
% Motivation: Gallen & Fernandez-Blanco (2021) allow spatially varying
% uplift -- their S_U is a vector, a flexural profile
% Uf*exp(-x/lambda)*cos(x/lambda) for the Corinth rift.  The version here
% assumes uniform uplift, which cannot reproduce a channel whose relict
% upper reach is systematically lower-relief than a single steepness
% predicts.  A linear gradient is the minimal extension: one parameter,
% and it applies to BOTH phases because the misfit it targets lies in the
% relict (pre-capture) part of the profile.
%
% Inputs:
%   S      - TopoToolbox STREAMobj
%   U_grad - fractional change in uplift from outlet to the farthest
%            channel head.  0 = spatially uniform (the previous model);
%            -0.5 = head uplifts at half the outlet rate; +0.5 = 1.5x.
%
% Outputs:
%   U_pat  - [N x 1] multiplier on the phase uplift rate, one per node
%
% Reference:
%   Gallen & Fernandez-Blanco (2021), JGR-Earth Surface

x    = double(S.distance(:));
xmax = max(x);

if ~isfinite(xmax) || xmax <= 0
    U_pat = ones(size(x));
    return;
end

U_pat = 1 + U_grad * (x / xmax);

% Keep strictly positive.  A negative uplift field would drive the
% steady-state solve to complex values via (U/K)^(1/n); the prior bound on
% U_grad should prevent this, but the clamp makes it structural.
U_pat = max(U_pat, 1e-6);

end
