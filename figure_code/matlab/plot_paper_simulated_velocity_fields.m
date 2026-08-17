function manifest = plot_paper_simulated_velocity_fields(varargin)
%PLOT_PAPER_SIMULATED_VELOCITY_FIELDS Rebuild simulated velocity snapshots.
%
% The catalog marks are stored in the scaled units used by simulation.m.
% This routine converts them back to grid and velocity units, applies linear
% amplitude decay over each eddy lifetime, and superposes all active eddies
% using the compactly supported basic-eddy profile defined in the paper.
%
% Name-value options:
%   CatalogPath    Constant-background simulated catalog .mat file.
%   OutputDir      Destination directory for the snapshot PDFs.
%   Snapshots      Time-block values to render (default [252 253]).
%   FilePattern    sprintf pattern for snapshot filenames.
%   ManifestFile  Optional CSV filename written inside OutputDir.
%   ProgressEvery Print progress every N snapshots.

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(scriptDir, '..', '..');

defaultCatalogPath = fullfile(repoRoot, 'exp_kernel', 'results', ...
    'Lifetime-first14days', 'validation', 'sim_catalog.mat');
defaultFigureDir = fullfile(repoRoot, 'figures');

p = inputParser;
addParameter(p, 'CatalogPath', defaultCatalogPath, ...
    @(v) ischar(v) || (isstring(v) && isscalar(v)));
addParameter(p, 'OutputDir', defaultFigureDir, ...
    @(v) ischar(v) || (isstring(v) && isscalar(v)));
addParameter(p, 'Snapshots', [252 253], ...
    @(v) isnumeric(v) && isvector(v) && all(isfinite(v)));
addParameter(p, 'FilePattern', 'simulated_191-%d.pdf', ...
    @(v) ischar(v) || (isstring(v) && isscalar(v)));
addParameter(p, 'ManifestFile', '', ...
    @(v) ischar(v) || (isstring(v) && isscalar(v)));
addParameter(p, 'ProgressEvery', 1, ...
    @(v) isnumeric(v) && isscalar(v) && isfinite(v) && v >= 1);
parse(p, varargin{:});

catalogPath = char(p.Results.CatalogPath);
figureDir = char(p.Results.OutputDir);
snapshots = double(p.Results.Snapshots(:)');
filePattern = char(p.Results.FilePattern);
manifestFile = char(p.Results.ManifestFile);
progressEvery = max(1, round(p.Results.ProgressEvery));

gridValues = 2:2:90; % Matches the observed 191-*.txt radar grid.
axisPadding = 25;
velocityScale = 1;   % One plot unit equals one velocity unit.

if ~isfile(catalogPath)
    error('Simulated catalog not found: %s', catalogPath);
end
if ~isfolder(figureDir)
    mkdir(figureDir);
end

loaded = load(catalogPath, 'catalog');
catalog = loaded.catalog;
[X, Y] = meshgrid(gridValues, gridValues);
nSnapshots = numel(snapshots);
nActiveValues = zeros(nSnapshots, 1);
maxSpeedValues = zeros(nSnapshots, 1);
filePaths = strings(nSnapshots, 1);

fprintf('Source catalog: %s\n', catalogPath);
for k = 1:numel(snapshots)
    snapshot = snapshots(k);
    [U, V, nActive] = reconstructVelocityField(catalog, snapshot, X, Y);

    figurePath = fullfile(figureDir, sprintf(filePattern, snapshot));
    exportVelocityFigure(figurePath, X, Y, U, V, ...
        velocityScale, axisPadding);

    nActiveValues(k) = nActive;
    maxSpeedValues(k) = max(hypot(U(:), V(:)));
    filePaths(k) = string(figurePath);
    if mod(k, progressEvery) == 0 || k == 1 || k == nSnapshots
        fprintf(['Rendered %d/%d (block %g): %d active eddies, ' ...
            'maximum speed %.6g\n'], k, nSnapshots, snapshot, nActive, ...
            maxSpeedValues(k));
    end
end

manifest = table(snapshots(:), snapshots(:) / 96, snapshots(:) * 0.25, ...
    nActiveValues, maxSpeedValues, filePaths, 'VariableNames', ...
    {'snapshot_block','day','time_hours','active_eddies','max_speed','pdf_file'});
if ~isempty(manifestFile)
    writetable(manifest, fullfile(figureDir, manifestFile));
end
end

function [U, V, nActive] = reconstructVelocityField(catalog, snapshot, X, Y)
time = catalog.time(:);
lifetime = max(catalog.lifetime(:), 1);
xCenter = catalog.x(:);
yCenter = catalog.y(:);

% simulation.m stores amplitude multiplied by 9 and radius in metres.
amplitude = catalog.aa(:) / 9;
radius = max(catalog.bb(:) / 125, 0.5);

active = find(time <= snapshot & snapshot < time + lifetime);
nActive = numel(active);
U = zeros(size(X));
V = zeros(size(Y));

for k = 1:nActive
    i = active(k);
    age = snapshot - time(i);
    amplitudeNow = amplitude(i) * max(1 - age / lifetime(i), 0);
    radiusNow = radius(i);

    xiX = (X - xCenter(i)) / radiusNow;
    xiY = (Y - yCenter(i)) / radiusNow;
    rho = hypot(xiX, xiY);
    inside = rho > 0 & rho <= 1;

    magnitude = zeros(size(rho));
    magnitude(inside) = (1 - cos(2 * pi * rho(inside))) / 2;
    rotation = zeros(size(rho));
    rotation(inside) = amplitudeNow .* magnitude(inside) ./ rho(inside);

    U = U - rotation .* xiY;
    V = V + rotation .* xiX;
end
end

function exportVelocityFigure(filePath, X, Y, U, V, velocityScale, axisPadding)
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Units', 'pixels', 'Position', [100 100 760 720]);
ax = axes(fig);
if isprop(ax, 'Toolbar')
    ax.Toolbar.Visible = 'off';
end

quiver(ax, X, Y, U / velocityScale, V / velocityScale, 0, ...
    'Color', [0 0 0], ...
    'LineWidth', 0.35, ...
    'MaxHeadSize', 0.45);
axis(ax, 'equal');
axis(ax, [-axisPadding, 90 + axisPadding, ...
    -axisPadding, 90 + axisPadding]);
axis(ax, 'off');

exportgraphics(fig, filePath, 'ContentType', 'image', 'Resolution', 220);
close(fig);
end
