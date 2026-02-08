function [stream_data] = prepare_hc_stream_data(dem_path, varargin)
% PREPARE_HC_STREAM_DATA  Extract stream network data from a DEM for use
% in the Hells Canyon Bayesian inversion.
%
% Can load either:
%   (1) A .mat file containing pre-extracted TopoToolbox objects
%   (2) A GeoTIFF DEM file
%
% Inputs:
%   dem_path    - Path to DEM (.tif) or pre-extracted data (.mat)
%
%   Optional name-value pairs:
%     'crita'       - Critical drainage area for channel heads (m^2, default: 1e6)
%     'mn'          - Reference concavity m/n (default: 0.45)
%     'outlet_xy'   - [x, y] coordinates of desired outlet (default: auto)
%     'dem_varname'  - Variable name for DEM in .mat file (default: 'DEMoc')
%     'S_varname'    - Variable name for STREAMobj (default: 'Sc')
%     'DA_varname'   - Variable name for drainage area (default: 'Ac')
%     'FD_varname'   - Variable name for FLOWobj (default: 'FDc')
%     'output_file'  - Path to save extracted data (default: 'hc_stream_data.mat')
%
% Outputs:
%   stream_data - struct with fields:
%     .S       - TopoToolbox STREAMobj
%     .Sz      - Elevations at stream nodes (m)
%     .S_DA    - Drainage area at stream nodes (m^2)
%     .S_dist  - Distance from outlet at each node (m)
%     .Schi    - Chi coordinate at each node
%     .mn      - m/n ratio used
%     .crita   - Critical area used
%
% Requires: TopoToolbox 2.4+

%% Parse inputs
p = inputParser;
addRequired(p, 'dem_path', @ischar);
addParameter(p, 'crita', 1e6, @isscalar);
addParameter(p, 'mn', 0.45, @isscalar);
addParameter(p, 'outlet_xy', [], @(x) isempty(x) || numel(x)==2);
addParameter(p, 'dem_varname', 'DEMoc', @ischar);
addParameter(p, 'S_varname', 'Sc', @ischar);
addParameter(p, 'DA_varname', 'Ac', @ischar);
addParameter(p, 'FD_varname', 'FDc', @ischar);
addParameter(p, 'output_file', 'hc_stream_data.mat', @ischar);

parse(p, dem_path, varargin{:});
opts = p.Results;

%% Load data
[~, ~, ext] = fileparts(dem_path);

if strcmpi(ext, '.mat')
    fprintf('Loading pre-extracted data from: %s\n', dem_path);
    data = load(dem_path);

    % Try to extract objects from .mat file
    if isfield(data, opts.S_varname) && isfield(data, opts.dem_varname)
        % Pre-extracted TopoToolbox objects available
        S   = data.(opts.S_varname);
        DEM = data.(opts.dem_varname);

        if isfield(data, opts.DA_varname)
            DA_grid = data.(opts.DA_varname);
            % DA might be a GRIDobj or a flow accumulation grid
            if isa(DA_grid, 'GRIDobj')
                S_DA = DA_grid.Z(S.IXgrid);
            else
                S_DA = DA_grid(S.IXgrid);
            end
        else
            % Need to compute DA from FD
            if isfield(data, opts.FD_varname)
                FD = data.(opts.FD_varname);
                DA = flowacc(FD) .* (FD.cellsize^2);
                S_DA = DA.Z(S.IXgrid);
            else
                error('Cannot compute drainage area: no FLOWobj or DA grid found.');
            end
        end

        % Get elevations
        if isa(DEM, 'GRIDobj')
            Sz = DEM.Z(S.IXgrid);
        else
            Sz = DEM(S.IXgrid);
        end

    elseif isfield(data, 'tDEM')
        % Sean Gallen's format
        DEM = data.tDEM;
        fprintf('Found tDEM format, extracting stream network...\n');
        DEM = fillsinks(DEM);
        FD = FLOWobj(DEM);
        DA = flowacc(FD) .* (FD.cellsize^2);
        S  = STREAMobj(FD, 'minarea', opts.crita / (DEM.cellsize)^2);
        S  = klargestconncomps(S, 1);
        Sz   = DEM.Z(S.IXgrid);
        S_DA = DA.Z(S.IXgrid);
    else
        error('Could not find expected variables in .mat file.');
    end

elseif strcmpi(ext, '.tif') || strcmpi(ext, '.tiff')
    fprintf('Loading GeoTIFF DEM from: %s\n', dem_path);
    DEM = GRIDobj(dem_path);
    DEM.Z(DEM.Z <= -9999) = NaN;

    DEM = fillsinks(DEM);
    FD  = FLOWobj(DEM);
    DA  = flowacc(FD) .* (FD.cellsize^2);
    S   = STREAMobj(FD, 'minarea', opts.crita / (DEM.cellsize)^2);
    S   = klargestconncomps(S, 1);

    % Optionally snap to a specific outlet
    if ~isempty(opts.outlet_xy)
        Sout = snap2stream(S, opts.outlet_xy(1), opts.outlet_xy(2));
        DB   = drainagebasins(FD, Sout);
        DEM.Z(DB.Z == 0) = NaN;
        DEM  = crop(DEM);
        FD   = FLOWobj(DEM);
        DA   = flowacc(FD) .* (FD.cellsize^2);
        S    = STREAMobj(FD, 'minarea', opts.crita / (DEM.cellsize)^2);
        S    = klargestconncomps(S, 1);
    end

    Sz   = DEM.Z(S.IXgrid);
    S_DA = DA.Z(S.IXgrid);
else
    error('Unsupported file format: %s. Use .mat or .tif', ext);
end

%% Compute chi
fprintf('Computing chi coordinates (m/n = %.2f)...\n', opts.mn);
Schi = compute_chi(S, S_DA, opts.mn);

%% Package output
stream_data = struct();
stream_data.S      = S;
stream_data.Sz     = Sz;
stream_data.S_DA   = S_DA;
stream_data.S_dist = S.distance;
stream_data.Schi   = Schi;
stream_data.mn     = opts.mn;
stream_data.crita  = opts.crita;

% Summary statistics
fprintf('Stream network extracted:\n');
fprintf('  Nodes: %d\n', length(Sz));
fprintf('  Elevation range: %.0f - %.0f m\n', min(Sz), max(Sz));
fprintf('  Chi range: %.1f - %.1f\n', min(Schi), max(Schi));
fprintf('  Max drainage area: %.1f km^2\n', max(S_DA) / 1e6);

%% Save
if ~isempty(opts.output_file)
    save(opts.output_file, 'stream_data');
    fprintf('Stream data saved to: %s\n', opts.output_file);
end

end

%% ===================== Helper Functions =====================

function Schi = compute_chi(S, S_DA, mn)
% Compute chi coordinate for the stream network.
% chi = integral from outlet to x of (1/A)^(m/n) dx

Schi = zeros(size(S.distance));
Six  = S.ix;
Sixc = S.ixc;
Sx   = S.distance;
Sa   = (1 ./ S_DA).^mn;

for lp = numel(Six):-1:1
    Schi(Six(lp)) = Schi(Sixc(lp)) + ...
        (Sa(Sixc(lp)) + (Sa(Six(lp)) - Sa(Sixc(lp))) / 2) * ...
        abs(Sx(Sixc(lp)) - Sx(Six(lp)));
end

end
