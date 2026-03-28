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
% TopoToolbox note: S.ix and S.ixc are linear indices into the DEM grid.
% We build a grid-to-node mapping so all internal arrays use compact
% node-indexed vectors (N nodes) rather than full grid arrays (~7x fewer
% elements), significantly improving speed and memory.
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
%
% Reference:
%   Braun & Willett (2013), Geomorphology
%   Gallen & Fernandez-Blanco (2021), JGR-Earth Surface

%% Build grid-to-node index mapping
% S.ix/S.ixc are grid linear indices. We convert them to node indices
% (1..N) so all arrays can be compact N-element vectors.
S_DA = S_DA(:);
N = numel(S.IXgrid);

grid2node = zeros(max(S.IXgrid), 1, 'int32');
grid2node(S.IXgrid) = int32(1:N)';

% Convert edge lists from grid indices to node indices
d_node = grid2node(S.ix);   % donor (upstream) node indices
r_node = grid2node(S.ixc);  % receiver (downstream) node indices

% Convert outlet grid indices to node indices
outlet_grid = streampoi(S, 'outlets', 'ix');
outlet_node = grid2node(outlet_grid);

% Node-indexed distance
dist_node = S.distance(:);  % already N elements, node-indexed

% Edge distances
dx = abs(dist_node(d_node) - dist_node(r_node));

% Number of edges
n_edges = numel(d_node);

%% Compute initial steady-state profile with U_pre
mn = m / n;
Z_mod = calculate_z_node(S_DA, dist_node, d_node, r_node, n_edges, ...
    U_pre, K, mn, n, N);
Z_mod = check_z_node(Z_mod, d_node, r_node, n_edges);

%% Prepare erosion velocity field: K * A^m at each node
if isscalar(K)
    Af = K * S_DA.^m;
else
    Af = K(:) .* S_DA.^m;
end

%% Precompute Af * dt / dx for each edge (constant across time steps)
tt_base = Af(d_node) ./ dx;  % K*A^m/dx at each edge

%% Time-stepping forward model
n_steps = ceil(run_time / dt);
n_outlets = length(outlet_node);

for step = 1:n_steps
    time_bp = run_time - step * dt;

    if time_bp >= t_capture
        U_current = U_pre;
    else
        U_current = U_post;
    end

    % Apply uplift
    Z_mod = Z_mod + U_current * dt;

    % Fix outlet elevations at base level
    for k = 1:n_outlets
        Z_mod(outlet_node(k)) = 0;
    end

    % Erode: implicit scheme, downstream to upstream
    tt_edge = tt_base * dt;  % K*A^m*dt/dx for each edge

    if abs(n - 1) < 1e-6
        % n = 1: vectorized implicit solution
        Z_mod(d_node) = (Z_mod(d_node) + Z_mod(r_node) .* tt_edge) ./ ...
            (1 + tt_edge);
    else
        % n != 1: Newton-Raphson at each edge
        for j = 1:n_edges
            Z_mod(d_node(j)) = newtonraphson(Z_mod(d_node(j)), ...
                Z_mod(r_node(j)), dx(j), tt_edge(j), n);
        end
    end
end

Z_mod = Z_mod(:);

end

%% ===================== Helper Functions =====================

function z = calculate_z_node(S_DA, dist_node, d_node, r_node, n_edges, ...
    U, K, mn, n_exp, N)
% Steady-state elevation via chi integration, node-indexed.

if isscalar(K)
    Sa = (U / K)^(1/n_exp) * (1 ./ S_DA).^mn;
else
    Sa = (U ./ K(:)).^(1/n_exp) .* (1 ./ S_DA).^mn;
end

z = zeros(N, 1);

for lp = n_edges:-1:1
    z(d_node(lp)) = z(r_node(lp)) + ...
        (Sa(r_node(lp)) + (Sa(d_node(lp)) - Sa(r_node(lp))) / 2) * ...
        abs(dist_node(r_node(lp)) - dist_node(d_node(lp)));
end

end

function z = check_z_node(z, d_node, r_node, n_edges)
% Ensure no node is lower than its downstream neighbor.

for lp = n_edges:-1:1
    if z(d_node(lp)) < z(r_node(lp))
        z(d_node(lp)) = z(r_node(lp)) + 0.001;
    end
end

end

function ztp1 = newtonraphson(zt, ztp1d, dx, tt, n)
% Newton-Raphson solver for implicit stream power equation when n != 1.

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
