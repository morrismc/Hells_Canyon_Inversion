function [cave_heights_pred] = cave_forward_model(cave_ages, U_rates, ...
    t_transitions)
% CAVE_FORWARD_MODEL  Predict the height of cave sediment deposits above
% the modern river level given a multi-phase incision history.
%
% General model for N phases:
%   Phase 1 (oldest):  t > t_transitions(1),              rate = U_rates(1)
%   Phase 2:           t_transitions(1) > t > t_transitions(2), rate = U_rates(2)
%   ...
%   Phase N (present): t < t_transitions(end),             rate = U_rates(end)
%
% For each cave, the predicted height above the modern river is the total
% incision accumulated from the cave's abandonment age to the present,
% i.e.  sum over all phases of  U_phase * (time spent in that phase).
%
% Inputs:
%   cave_ages       - Vector of cave burial ages (years before present)
%   U_rates         - Vector of incision rates from oldest to youngest (m/yr)
%                     2-phase: [U_pre, U_post]
%                     3-phase: [U_pre, U_mid, U_post]
%   t_transitions   - Transition times in years BP, strictly decreasing.
%                     2-phase: [t_capture]
%                     3-phase: [t1, t2]  where t1 > t2.
%
% Outputs:
%   cave_heights_pred - Predicted height above modern river (m)
%
% Reference:
%   Morriss et al. (2025), PNAS: Cave burial dating of Hells Canyon

cave_heights_pred = zeros(size(cave_ages));

% Phase boundaries: [Inf, t1, t2, ..., 0]
t_bounds = [Inf, t_transitions(:)', 0];
n_phases = length(U_rates);

for i = 1:length(cave_ages)
    A = cave_ages(i);
    h = 0;
    for p = 1:n_phases
        upper = min(A, t_bounds(p));
        lower = t_bounds(p + 1);
        if upper > lower
            h = h + U_rates(p) * (upper - lower);
        end
    end
    cave_heights_pred(i) = h;
end

end
