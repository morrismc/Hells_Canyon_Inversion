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
% IMPORTANT: TopoToolbox STREAMobj properties S.ix and S.ixc are linear
% indices into the DEM grid, NOT indices into the stream node list.
% All internal arrays are therefore allocated in grid-index space, and
% stream-node values are extracted at the end via S.IXgrid.
%
% Inputs:
%   S          - TopoToolbox STREAMobj
%   S_DA       - Drainage area at each stream node (m^2), vector (N nodes)
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
%   Z_mod      - Modeled elevations at each stream node (m), column vector
%                length = numel(S.IXgrid)
%
% Adapted from Gallen & Fernandez-Blanco (2021) bayes_profiler code.
% Modified for Hells Canyon drainage capture scenario.
%
% Reference:
%   Braun & Willett (2013), Geomorphology
%   Gallen & Fernandez-Blanco (2021), JGR-Earth Surface

%% Map node-indexed inputs to grid-index space
% S.ix, S.ixc, and S.IXgrid are all linear indices into the DEM grid.
% We must work in grid space so that indexing with S.ix/S.ixc is valid.
S_DA = S_DA(:);
n_grid = max(S.IXgrid);  % size of grid-index space

% Map drainage area from node space to grid space
DA_grid = zeros(n_grid, 1);
DA_grid(S.IXgrid) = S_DA;

% Map distances from node space to grid space
dist_grid = zeros(n_grid, 1);
dist_grid(S.IXgrid) = S.distance;

% Map K to grid space
if isscalar(K)
    K_grid = K * ones(n_grid, 1);
else
    K_grid = zeros(n_grid, 1);
    K_grid(S.IXgrid) = K(:);
end

%% Compute initial steady-state profile with U_pre
mn = m / n;
Z_grid = calculate_z(S, DA_grid, dist_grid, U_pre, K, mn, n);
Z_grid = check_z(S, Z_grid);

%% Prepare the implicit finite-difference solver
[d, r, Af, dx, outlet_nodes] = fastscape_eroder_data_prep(S, DA_grid, K_grid, m, dist_grid);

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

    % Apply uplift to stream nodes only
    Z_grid(S.IXgrid) = Z_grid(S.IXgrid) + U_current * dt;

    % Apply erosion using implicit scheme
    Z_grid = fastscape_eroder_outlets(Z_grid, n, dt, Af, d, r, dx, ...
        U_current, outlet_nodes);
end

%% Extract stream node elevations from grid space
Z_mod = Z_grid(S.IXgrid);
Z_mod = Z_mod(:);

end

%% ===================== Helper Functions =====================

function z = calculate_z(S, DA_grid, dist_grid, U, K, mn, n)
% Calculate steady-state elevation profile using chi integration.
% Works in grid-index space.

if isscalar(K)
    Sa = zeros(size(DA_grid));
    Sa(S.IXgrid) = (U / K)^(1/n) * (1 ./ DA_grid(S.IXgrid)).^mn;
else
    Sa = zeros(size(DA_grid));
    Sa(S.IXgrid) = (U ./ DA_grid(S.IXgrid)).^(1/n) .* (1 ./ DA_grid(S.IXgrid)).^mn;
end

z   = zeros(size(DA_grid));
Six  = S.ix;
Sixc = S.ixc;

% Integrate upstream from outlet
for lp = numel(Six):-1:1
    z(Six(lp)) = z(Sixc(lp)) + ...
        (Sa(Sixc(lp)) + (Sa(Six(lp)) - Sa(Sixc(lp))) / 2) * ...
        abs(dist_grid(Sixc(lp)) - dist_grid(Six(lp)));
end

end

function z = check_z(S, z)
% Ensure no node is lower than its downstream neighbor.
Six  = S.ix;
Sixc = S.ixc;

for lp = numel(Six):-1:1
    if z(Six(lp)) < z(Sixc(lp))
        z(Six(lp)) = z(Sixc(lp)) + 0.001;
    end
end

end

function [d, r, A, dx, outlet_nodes] = fastscape_eroder_data_prep(S, DA_grid, K_grid, m, dist_grid)
% Prepare data structures for the implicit finite difference solver.
% All arrays are in grid-index space.

d  = S.ix;
r  = S.ixc;
dx = abs(dist_grid(d) - dist_grid(r));

% Erosion velocity field: K * A^m (grid-indexed)
A = K_grid .* DA_grid.^m;

% Find outlet nodes (returns grid linear indices)
outlet_nodes = streampoi(S, 'outlets', 'ix');

end

function S_Z = fastscape_eroder_outlets(S_Z, n, dt, A, d, r, dx, U, outlets)
% Implicit finite-difference erosion following Braun & Willett (2013).
% All arrays are in grid-index space.

% Fix outlet elevations at base level = 0
for k = 1:length(outlets)
    S_Z(outlets(k)) = 0;
end

% Erode from downstream to upstream (implicit scheme).
% TopoToolbox orders S.ix/S.ixc from headwaters (index 1) to outlet
% (index end), so we iterate backwards to process outlet-to-headwaters,
% ensuring each downstream node is updated before its upstream donors.
for j = length(d):-1:1
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

ztp1 = zt;
tol  = 1e-3;
max_iter = 100;

for iter = 1:max_iter
    slope = (ztp1 - ztp1d) / dx;

    if slope <= 0
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

    if ztp1 < ztp1d
        ztp1 = ztp1d + 0.001;
    end
end

end
