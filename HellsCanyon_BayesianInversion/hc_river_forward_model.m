function [Z_mod] = hc_river_forward_model(S, S_DA, U_rates, ...
    t_transitions, K, m, n, dt)
% HC_RIVER_FORWARD_MODEL  Forward model for river incision following
% drainage capture in Hells Canyon.
%
% Solves the stream power incision model:
%   dz/dt = U(t) - K * A^m * |dz/dx|^n
%
% Multi-phase tectonic scenario (general):
%   Phase 1 (oldest):  t > t_transitions(1),  rate = U_rates(1)
%   Phase 2:           t_transitions(1) > t > t_transitions(2), rate = U_rates(2)
%   ...
%   Phase N (present):  t < t_transitions(end), rate = U_rates(end)
%
% The model starts from a steady-state profile under U_rates(1), then
% steps through each subsequent phase to the present.
%
% Uses the implicit finite difference scheme of Braun & Willett (2013)
% with Newton-Raphson iteration for n != 1.
%
% All internal arrays are in NODE space (length = numel(S.IXgrid)),
% matching the TopoToolbox convention where S.ix and S.ixc are indices
% into the stream-node list, not into the DEM grid.
%
% Inputs:
%   S              - TopoToolbox STREAMobj
%   S_DA           - Drainage area at each stream node (m^2), vector
%   U_rates        - Vector of uplift rates from oldest to youngest (m/yr).
%                    For 2-phase: [U_pre, U_post].
%                    For 3-phase: [U_pre, U_mid, U_post].
%   t_transitions  - Transition times in years BP, strictly decreasing.
%                    For 2-phase: [t_capture].
%                    For 3-phase: [t1, t2]  where t1 > t2.
%   K              - Erodibility coefficient, scalar or vector (N nodes)
%   m              - Drainage area exponent, scalar
%   n              - Slope exponent, scalar
%   dt             - Time step (years)
%
% Outputs:
%   Z_mod          - Modeled elevations at each stream node (m), column vector
%
% Closely follows Gallen & Fernandez-Blanco (2021) bayes_profiler
% river_incision_forward_model.m, extended for multi-phase capture
% scenario.
%
% Reference:
%   Braun & Willett (2013), Geomorphology
%   Gallen & Fernandez-Blanco (2021), JGR-Earth Surface

S_DA = S_DA(:);
U_rates = U_rates(:)';
t_transitions = t_transitions(:)';

%% Compute initial steady-state profile with the oldest (background) rate
mn = m / n;
Z = calculate_z(S, S_DA, U_rates(1), K, mn, n);
Z = check_z(S, Z);

%% Prepare the implicit finite-difference solver
[d, r, Af, dx, outlet_nodes] = fastscape_eroder_data_prep(S, S_DA, K, m);

%% Diagnostic output (first call only)
persistent diag_printed;
if isempty(diag_printed)
    diag_printed = true;
    fprintf('\n--- Forward Model Diagnostics ---\n');
    fprintf('  Stream nodes: %d\n', numel(S.IXgrid));
    fprintf('  Uplift phases: %d  (rates: %s m/yr)\n', ...
        length(U_rates), mat2str(U_rates, 4));
    fprintf('  Transition times: %s yr BP\n', mat2str(t_transitions));
    fprintf('  Outlet nodes found: %d  (node indices: %s)\n', ...
        numel(outlet_nodes), mat2str(outlet_nodes(:)'));
    fprintf('  Outlet elevations in Z: %s\n', ...
        mat2str(Z(outlet_nodes)', 4));
    fprintf('  Initial steady-state relief: %.1f m\n', max(Z) - min(Z));
    fprintf('  dx range: [%.1f, %.1f] m\n', min(dx(dx>0)), max(dx));
    fprintf('  max(K*A^m): %.2e   min(K*A^m): %.2e\n', max(Af), min(Af));
    fprintf('--- End Diagnostics ---\n\n');
end

%% Time-stepping: march through each phase from oldest to present
% t_bounds marks the interval endpoints: [t1, t2, ..., 0]
t_bounds = [t_transitions, 0];

for phase = 2:length(U_rates)
    duration = t_bounds(phase - 1) - t_bounds(phase);
    n_steps  = max(1, ceil(duration / dt));
    dt_phase = duration / n_steps;   % keep exact duration

    for step = 1:n_steps
        Z = fastscape_eroder_outlets(Z, n, dt_phase, Af, d, r, dx, ...
            U_rates(phase), outlet_nodes);
    end
end

% Final safety check (no-op if the implicit solver behaved correctly)
Z = check_z(S, Z);

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
outlet_grid = streampoi(S, 'outlets', 'ix');
outlet_nodes = zeros(size(outlet_grid));
for i = 1:length(outlet_grid)
    outlet_nodes(i) = find(S.IXgrid == outlet_grid(i));
end

% Distance between donor-receiver pairs.  Guard against the rare case
% where two consecutive stream nodes end up at identical distance
% values (which would produce tt = A*dt/dx = Inf in the solver).
dx = abs(S.distance(d) - S.distance(r));
dx_min = eps(max(S.distance)) * 10;
dx(dx < dx_min) = dx_min;

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
% Matches Gallen's implementation but with a max-iteration safeguard and
% clamping to prevent tempz from dropping below the receiver (which would
% make ((tempz-ztp1d)/dx)^n complex for non-integer n).

tempz   = zt;
tol     = inf;
max_it  = 50;
it      = 0;

while tol > 1e-3 && it < max_it
    it = it + 1;

    dz_dx = (tempz - ztp1d) / dx;
    if dz_dx <= 0
        tempz = ztp1d + 1e-3;
        dz_dx = 1e-3 / dx;
    end

    f_val   = tempz - zt + (tt * dx) * dz_dx^n;
    f_prime = 1 + n * tt * dz_dx^(n - 1);

    ztp1  = tempz - f_val / f_prime;
    tol   = abs(ztp1 - tempz);
    tempz = ztp1;
end

if ~isreal(ztp1) || isnan(ztp1) || isinf(ztp1) || ztp1 < ztp1d
    ztp1 = max(real(ztp1), ztp1d + 1e-3);
end

end
