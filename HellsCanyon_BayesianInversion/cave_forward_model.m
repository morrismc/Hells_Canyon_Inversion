function [cave_heights_pred] = cave_forward_model(cave_ages, t_capture, ...
    U_pre, U_post)
% CAVE_FORWARD_MODEL  Predict the height of cave sediment deposits above
% the modern river level given a two-phase incision history.
%
% Model:
%   Before capture (cave_age > t_capture):
%     height = U_post * t_capture + U_pre * (cave_age - t_capture)
%
%   After capture (cave_age <= t_capture):
%     height = U_post * cave_age
%
% The logic: a cave at a given age records the river's position when
% the cave was active. The total drop below that position is the
% cumulative incision since the cave was abandoned.
%
% Inputs:
%   cave_ages   - Vector of cave burial ages (years before present)
%   t_capture   - Time of drainage capture (years before present)
%   U_pre       - Pre-capture incision rate (m/yr)
%   U_post      - Post-capture incision rate (m/yr)
%
% Outputs:
%   cave_heights_pred - Predicted height above modern river (m)
%
% Reference:
%   Morriss et al. (2025), PNAS: Cave burial dating of Hells Canyon

cave_heights_pred = zeros(size(cave_ages));

for i = 1:length(cave_ages)
    if cave_ages(i) > t_capture
        % Cave formed before capture event
        % Total incision = fast phase (t_capture duration) + slow phase (remainder)
        cave_heights_pred(i) = U_post * t_capture + ...
                               U_pre * (cave_ages(i) - t_capture);
    else
        % Cave formed after capture event - only fast phase
        cave_heights_pred(i) = U_post * cave_ages(i);
    end
end

end
