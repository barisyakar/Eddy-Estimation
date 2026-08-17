function plot_velocity_snapshot_pdfs(inputDir, outputDir, axisPadding)
%PLOT_VELOCITY_SNAPSHOT_PDFS Render all numbered x-y-u-v snapshot files.

if nargin < 1 || isempty(inputDir)
    scriptDir = fileparts(mfilename('fullpath'));
    repoRoot = fullfile(scriptDir, '..', '..');
    inputDir = fullfile(repoRoot, 'output', 'velocity_fields');
end
if nargin < 2 || isempty(outputDir)
    outputDir = [char(inputDir) '_pdfs'];
end
if nargin < 3 || isempty(axisPadding)
    axisPadding = 25;
end
if isstring(inputDir)
    inputDir = char(inputDir);
end
if isstring(outputDir)
    outputDir = char(outputDir);
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

files = dir(fullfile(inputDir, '191-*.txt'));
files = files(arrayfun(@(f) ...
    ~isempty(regexp(f.name, '^191-\d+\.txt$', 'once')), files));
if isempty(files)
    error('No numbered 191-*.txt files found in %s', inputDir);
end
files = sortSnapshotFiles(files);

for k = 1:numel(files)
    txtPath = fullfile(inputDir, files(k).name);
    raw = readmatrix(txtPath);
    validateSnapshot(raw, txtPath);

    valid = all(isfinite(raw(:, 1:4)), 2);
    x = raw(valid, 1);
    y = raw(valid, 2);
    u = raw(valid, 3);
    v = raw(valid, 4);

    fig = makeVelocityFigure(x, y, u, v, axisPadding, false);
    [~, baseName] = fileparts(files(k).name);
    exportgraphics(fig, fullfile(outputDir, [baseName '.pdf']), ...
        'ContentType', 'image', 'Resolution', 220);
    close(fig);

    if mod(k, 25) == 0 || k == numel(files)
        fprintf('Exported %d/%d PDFs\n', k, numel(files));
    end
end

fprintf('PDF snapshots written to %s\n', outputDir);
end

function files = sortSnapshotFiles(files)
numbers = zeros(numel(files), 1);
for k = 1:numel(files)
    token = regexp(files(k).name, '^191-(\d+)\.txt$', 'tokens', 'once');
    numbers(k) = str2double(token{1});
end
[~, order] = sort(numbers);
files = files(order);
end

function validateSnapshot(raw, filePath)
if size(raw, 2) < 4
    error('Expected x, y, u, v columns in %s', filePath);
end
end

function fig = makeVelocityFigure(x, y, u, v, axisPadding, visible)
visibility = 'off';
if visible
    visibility = 'on';
end
fig = figure('Visible', visibility, 'Color', 'w', ...
    'Units', 'pixels', 'Position', [100 100 760 720]);
ax = axes(fig);
if isprop(ax, 'Toolbar')
    ax.Toolbar.Visible = 'off';
end
quiver(ax, x, y, u, v, 0, 'Color', [0 0 0], ...
    'LineWidth', 0.35, 'MaxHeadSize', 0.45);
axis(ax, 'equal');
axis(ax, [-axisPadding 90 + axisPadding -axisPadding 90 + axisPadding]);
axis(ax, 'off');
end
