function [Z_mod] = hc_river_forward_model(S, S_DA, U_pre, U_post, ...
    t_capture, K, m, n, run_time, dt)
% HC_RIVER_FORWARD_MODEL  Forward model for river incision following
% drainage capture in Hells Canyon.
%
% Solves the stream power incision model:
%   dz/dt = U(t) - K * A^m * |dz/dx|^n
%
% Two-phase tectonic scenario:
%   Phase 1 (t < t_capture): Uniform incision at rate U_pre
%   Phase 2 (t >= t_capture): Uniform incision at rate U_post
%
% Uses the implicit finite difference scheme of Braun & Willett (2013)
% with Newton-Raphson iteration for n != 1.
%
% Inputs:
%   S          - TopoToolbox STREAMobj
%   S_DA       - Drainage area at each stream node (m^2), vector
%   U_pre      - Pre-capture incision/uplift rate (m/yr), scalar
%   U_post     - Post-capture incision/uplift rate (m/yr), scalar
%   t_capture  - Time of capture event (years before present)
%   K          - Erodibility coefficient (m^(1-2m)/yr), scalar or vector
%   m          - Drainage area exponent, scalar
%   n          - Slope exponent, scalar
%   run_time   - Total model run time (years), should be >= t_capture
%   dt         - Time step (years)
%
% Outputs:
%   Z_mod      - Modeled elevations at each stream node (m)
%
% Adapted from Gallen & Fernandez-Blanco (2021) bayes_profiler code.
% Modified for Hells Canyon drainage capture scenario.
%
% Reference:
%   Braun & Willett (2013), Geomorphology
%   Gallen & Fernandez-Blanco (2021), JGR-Earth Surface

%% Compute initial steady-state profile with U_pre
% At steady state: U_pre = K * A^m * S^n
% => S = (U_pre / K)^(1/n) * A^(-m/n)
% => z = integral of S dx from outlet upstream
mn = m / n;

Z_mod = calculate_z(S, S_DA, U_pre, K, mn, n);
Z_mod = check_z(S, Z_mod);

%% Prepare the implicit finite-difference solver
% Precompute the "velocity field" = K * A^m for the stream power term
if isscalar(K)
    S_K = K * ones(size(S_DA));
else
    S_K = K;
end

[d, r, Af, dx, outlet_nodes] = fastscape_eroder_data_prep(S, S_DA, S_K, m);

%% Time-stepping forward model
n_steps = ceil(run_time / dt);

for step = 1:n_steps
    current_time = step * dt;  % time elapsed from start of model

    % Determine which uplift rate applies
    % t_capture is measured as "years before present"
    % run_time should equal or exceed t_capture
    % time_before_present = run_time - current_time
    time_bp = run_time - current_time;

    if time_bp >= t_capture
        % We're in the pre-capture phase (before the capture event)
        U_current = U_pre;
    else
        % Post-capture phase
        U_current = U_post;
    end

    % Apply uplift to all nodes
    Z_mod = Z_mod + U_current * dt;

    % Apply erosion using implicit scheme
    Z_mod = fastscape_eroder_outlets(Z_mod, n, dt, Af, d, r, dx, ...
        U_current, outlet_nodes);
end

end

%% ===================== Helper Functions =====================

function z = calculate_z(S, S_DA, U, K, mn, n)
% Calculate steady-state elevation profile using chi integration.
% z = integral from outlet to x of: (U/K)^(1/n) * (1/A)^(m/n) dx
%
% For a scalar K:
%   z(x) = (U/K)^(1/n) * integral of (1/A)^(m/n) dx

if isscalar(K)
    Sa = (U / K)^(1/n) * (1 ./ S_DA).^mn;
else
    Sa = (U ./ K).^(1/n) .* (1 ./ S_DA).^mn;
end

z = zeros(size(S.distance));
Six  = S.ix;
Sixc = S.ixc;
Sx   = S.distance;

% Integrate upstream from outlet
for lp = numel(Six):-1:1
    z(Six(lp)) = z(Sixc(lp)) + ...
        (Sa(Sixc(lp)) + (Sa(Six(lp)) - Sa(Sixc(lp))) / 2) * ...
        abs(Sx(Sixc(lp)) - Sx(Six(lp)));
end

end

function z = check_z(S, z)
% Ensure no node is lower than its downstream neighbor.
% Prevents backward-flowing rivers from initial conditions.
Six  = S.ix;
Sixc = S.ixc;

for lp = numel(Six):-1:1
    if z(Six(lp)) < z(Sixc(lp))
        z(Six(lp)) = z(Sixc(lp)) + 0.001;
    end
end

end

function [d, r, A, dx, outlet_nodes] = fastscape_eroder_data_prep(S, S_DA, S_K, m)
% Prepare data structures for the implicit finite difference solver.
%
% Outputs:
%   d   - donor (upstream) node indices
%   r   - receiver (downstream) node indices
%   A   - K * DA^m at each node (erosion velocity field)
%   dx  - distance between donor-receiver pairs
%   outlet_nodes - indices of outlet nodes

d  = S.ix;
r  = S.ixc;
dx = abs(S.distance(d) - S.distance(r));

% Erosion velocity field: K * A^m
A = S_K .* S_DA.^m;

% Find outlet nodes
outlet_nodes = streampoi(S, 'outlets', 'ix');

end

function S_Z = fastscape_eroder_outlets(S_Z, n, dt, A, d, r, dx, U, outlets)
% Implicit finite-difference erosion following Braun & Willett (2013).
%
% Fixes outlet elevations (base level = 0) and solves upstream.
% For n=1: direct algebraic solution.
% For n!=1: Newton-Raphson iterative solution at each node.

% Fix outlet elevations at base level = 0
for k = 1:length(outlets)
    S_Z(outlets(k)) = 0;
end

% Erode from downstream to upstream (implicit scheme)
for j = 1:length(d)
    % tt is the dimensionless erosion coefficient
    tt = A(d(j)) * dt / dx(j);

    if abs(n - 1) < 1e-6
        % n = 1: Direct implicit solution
        S_Z(d(j)) = (S_Z(d(j)) + S_Z(r(j)) * tt) / (1 + tt);
    else
        % n != 1: Newton-Raphson iteration
        S_Z(d(j)) = newtonraphson(S_Z(d(j)), S_Z(r(j)), dx(j), tt, n);
    end
end

end

function ztp1 = newtonraphson(zt, ztp1d, dx, tt, n)
% Newton-Raphson solver for the implicit stream power equation when n != 1.
%
% Solves: f(z) = z - zt + tt * dx * ((z - ztp1d) / dx)^n = 0
%
% Inputs:
%   zt    - elevation at current time step (before erosion)
%   ztp1d - downstream node elevation (already updated)
%   dx    - distance between nodes
%   tt    - K * A^m * dt / dx
%   n     - slope exponent
%
% Output:
%   ztp1  - updated elevation at next time step

ztp1 = zt;  % initial guess
tol  = 1e-3;  % convergence tolerance (meters)
max_iter = 100;

for iter = 1:max_iter
    slope = (ztp1 - ztp1d) / dx;

    if slope <= 0
        % If slope is negative or zero, set to small positive value
        ztp1 = ztp1d + 0.001;
        slope = 0.001 / dx;
    end

    f  = ztp1 - zt + tt * dx * slope^n;
    fp = 1 + n * tt * slope^(n-1);

    ztp1_new = ztp1 - f / fp;

    if abs(ztp1_new - ztp1) < tol
        ztp1 = ztp1_new;
        return;
    end

    ztp1 = ztp1_new;

    % Ensure elevation stays above downstream node
    if ztp1 < ztp1d
        ztp1 = ztp1d + 0.001;
    end
end

% If not converged, use last estimate
end
