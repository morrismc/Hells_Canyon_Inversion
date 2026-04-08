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
% All internal arrays are in NODE space (length = numel(S.IXgrid)),
% matching the TopoToolbox convention where S.ix and S.ixc are indices
% into the stream-node list, not into the DEM grid.
%
% Inputs:
%   S          - TopoToolbox STREAMobj
%   S_DA       - Drainage area at each stream node (m^2), vector (N nodes)
%   U_pre      - Pre-capture incision/uplift rate (m/yr), scalar
%   U_post     - Post-capture incision/uplift rate (m/yr), scalar
%   t_capture  - Time of capture event (years before present)
%   K          - Erodibility coefficient, scalar or vector (N nodes)
%   m          - Drainage area exponent, scalar
%   n          - Slope exponent, scalar
%   run_time   - Total model run time (years), should be >= t_capture
%   dt         - Time step (years)
%
% Outputs:
%   Z_mod      - Modeled elevations at each stream node (m), column vector
%                length = numel(S.IXgrid)
%
% Closely follows Gallen & Fernandez-Blanco (2021) bayes_profiler
% river_incision_forward_model.m, adapted for Hells Canyon two-phase
% capture scenario.
%
% Reference:
%   Braun & Willett (2013), Geomorphology
%   Gallen & Fernandez-Blanco (2021), JGR-Earth Surface

S_DA = S_DA(:);

%% Compute initial steady-state profile with U_pre
mn = m / n;
Z = calculate_z(S, S_DA, U_pre, K, mn, n);
Z = check_z(S, Z);

%% Prepare the implicit finite-difference solver
[d, r, Af, dx, outlet_nodes] = fastscape_eroder_data_prep(S, S_DA, K, m);

%% Diagnostic output (first call only)
persistent diag_printed;
if isempty(diag_printed)
    diag_printed = true;
    fprintf('\n--- Forward Model Diagnostics ---\n');
    fprintf('  Stream nodes: %d\n', numel(S.IXgrid));
    fprintf('  Outlet nodes found: %d  (node indices: %s)\n', ...
        numel(outlet_nodes), mat2str(outlet_nodes(:)'));
    fprintf('  Outlet elevations in Z: %s\n', ...
        mat2str(Z(outlet_nodes)', 4));
    fprintf('  Initial steady-state relief: %.1f m\n', max(Z) - min(Z));
    fprintf('  dx range: [%.1f, %.1f] m\n', min(dx(dx>0)), max(dx));
    fprintf('  max(K*A^m): %.2e   min(K*A^m): %.2e\n', max(Af), min(Af));
    fprintf('--- End Diagnostics ---\n\n');
end

%% Time-stepping forward model
n_steps = ceil(run_time / dt);

for step = 1:n_steps
    current_time = step * dt;
    time_bp = run_time - current_time;

    if time_bp >= t_capture
        U_current = U_pre;
    else
        U_current = U_post;
    end

    % Uplift + erosion in one call (matching Gallen's implementation,
    % which applies uplift and then removes it at outlets to keep them
    % at base level).
    Z = fastscape_eroder_outlets(Z, n, dt, Af, d, r, dx, ...
        U_current, outlet_nodes);

    % Ensure river doesn't flow backwards
    Z = check_z(S, Z);
end

%% Return node-ordered elevations
Z_mod = Z(:);

end

%% ===================== Helper Functions =====================

function z = calculate_z(S, S_DA, U, K, mn, n)
% Calculate steady-state elevation profile using chi integration.
% All arrays in node space (same length as S.distance).

z  = zeros(size(S.distance));
Sa = (U / K)^(1/n) * (1 ./ S_DA).^mn;

Six  = S.ix;
Sixc = S.ixc;
Sd   = S.distance;

% Integrate upstream from outlet
for lp = numel(Six):-1:1
    z(Six(lp)) = z(Sixc(lp)) + ...
        (Sa(Sixc(lp)) + (Sa(Six(lp)) - Sa(Sixc(lp))) / 2) * ...
        abs(Sd(Sixc(lp)) - Sd(Six(lp)));
end

end

function z = check_z(S, z)
% Ensure no node is lower than its downstream neighbor.
Six  = S.ix;
Sixc = S.ixc;

for lp = numel(Six):-1:1
    if z(Six(lp)) <= z(Sixc(lp))
        z(Six(lp)) = z(Sixc(lp)) + 0.1;
    end
end

end

function [d, r, A, dx, outlet_nodes] = fastscape_eroder_data_prep(S, S_DA, K, m)
% Prepare data structures for the implicit finite difference solver.
% All arrays in node space.

d = S.ix;
r = S.ixc;

% Erosion velocity field: K * A^m (node-indexed)
if isscalar(K)
    A = K .* S_DA.^m;
else
    A = K(:) .* S_DA.^m;
end

% Find outlet nodes and convert to node-list indices.
% S.outlets (or streampoi 'outlets') returns grid-linear indices;
% convert to positions within S.IXgrid (node indices) following
% Gallen's fastscape_eroder_data_prep.
outlet_grid = streampoi(S, 'outlets', 'ix');
outlet_nodes = zeros(size(outlet_grid));
for i = 1:length(outlet_grid)
    outlet_nodes(i) = find(S.IXgrid == outlet_grid(i));
end

% Distance between donor-receiver pairs
dx = abs(S.distance(d) - S.distance(r));

end

function S_Z = fastscape_eroder_outlets(S_Z, n, dt, A, d, r, dx, U, outlets)
% Implicit finite-difference erosion following Braun & Willett (2013).
% Uplift is applied inside this function (matching Gallen's implementation).
% Outlets are kept at base level by removing the uplift increment.

% Apply uplift to all nodes
S_Z = S_Z + dt * U;

% Remove uplift at outlets to maintain base level
if isscalar(U)
    S_Z(outlets) = S_Z(outlets) - dt * U;
else
    S_Z(outlets) = S_Z(outlets) - dt * U(outlets);
end

% Implicit erosion from downstream to upstream
for j = numel(d):-1:1
    tt = A(d(j)) * dt / dx(j);

    if abs(n - 1) < 1e-6
        % n = 1: Direct implicit solution
        S_Z(d(j)) = (S_Z(d(j)) + S_Z(r(j)) * tt) / (1 + tt);
    else
        % n != 1: Newton-Raphson iteration
        zt    = S_Z(d(j));
        ztp1d = S_Z(r(j));
        if ztp1d < zt
            S_Z(d(j)) = newtonraphson(zt, ztp1d, dx(j), tt, n);
        end
    end
end

end

function ztp1 = newtonraphson(zt, ztp1d, dx, tt, n)
% Newton-Raphson solver for the implicit stream power equation when n != 1.
% Matches Gallen's implementation.

tempz = zt;
tol   = inf;

while tol > 1e-3
    ztp1 = tempz - (tempz - zt + (tt * dx) * ((tempz - ztp1d) / dx)^n) / ...
        (1 + n * tt * ((tempz - ztp1d) / dx)^(n-1));
    tol   = abs(ztp1 - tempz);
    tempz = ztp1;
end

if ~isreal(ztp1) || isnan(ztp1)
    ztp1 = real(ztp1);
end

end
