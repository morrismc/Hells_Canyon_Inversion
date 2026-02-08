function lr = logproposal_hc(p_1, p_2, step)
% LOGPROPOSAL_HC  Gaussian transition probability for Metropolis-Hastings.
%
% Since the proposal is symmetric (Gaussian centered on current state),
% the forward and reverse proposals are equal, and this term cancels in
% the MH acceptance ratio. Included for completeness.
%
% Inputs:
%   p_1   - [1 x nparams] starting parameter vector
%   p_2   - [1 x nparams] proposed parameter vector
%   step  - [1 x nparams] step sizes for each parameter
%
% Output:
%   lr    - Log of the transition probability q(p_2 | p_1)

lr = sum(-0.5 * ((p_2 - p_1) ./ step).^2);

end
