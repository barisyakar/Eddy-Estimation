function make_overleaf_figures()
%MAKE_OVERLEAF_FIGURES Rebuild the paper figures supported by this package.
%
% The two simulated velocity snapshots can be rebuilt from the included
% catalogs alone. Figures derived from observations require the Lifetime-*.txt
% source catalogs under data/. Existing PDFs remain unchanged when those
% inputs are unavailable.

rootDir = fileparts(mfilename('fullpath'));
figureCodeDir = fullfile(rootDir, 'figure_code', 'matlab');
constantModelDir = fullfile(rootDir, 'exp_kernel');
outDir = fullfile(rootDir, 'figures');

addpath(figureCodeDir);
addpath(constantModelDir);
if ~isfolder(outDir)
    mkdir(outDir);
end

% These two panels use only the included first-window simulation.
plot_paper_simulated_velocity_fields();

dataNames = {'Lifetime-first14days', 'Lifetime-last14days'};
dataFiles = cellfun(@(name) fullfile(rootDir, 'data', [name '.txt']), ...
    dataNames, 'UniformOutput', false);
missingData = dataFiles(~cellfun(@isfile, dataFiles));
if ~isempty(missingData)
    warning('make_overleaf_figures:sourceDataMissing', ...
        ['The simulated velocity panels were rebuilt. The remaining ' ...
         'observation-derived figures were retained because these source ' ...
         'inputs are absent:\n  %s'], strjoin(missingData, '\n  '));
    return;
end

eddyplot.PlotStyle.apply();
datasets = cell(1, numel(dataNames));
for k = 1:numel(dataNames)
    dataName = dataNames{k};
    data = eddyplot.EddyData(dataFiles{k}).load();
    datasets{k} = data;

    saveFigure(eddyplot.FigACF(data).draw(), outDir, ...
        sprintf('acf15min-%s', dataName));
    saveFigure(eddyplot.FigACFAggregated(data).draw(), outDir, ...
        sprintf('acf1_of_aggregatedwindows-%s', dataName));
    saveFigure(eddyplot.FigSpatialDensity(data).draw(), outDir, ...
        sprintf('fig_eddy_dist-%s', dataName));
    saveFigure( ...
        eddyplot.FigSpatialDensity(data, 60, 0, true, 10).draw(), ...
        outDir, sprintf('fig_eddy_dist_3d-%s', dataName));
    makeConstantIntensity(rootDir, dataName, data, outDir);
end

makeLifetimeSurvival(datasets, outDir);
copyTimeVaryingFigures(rootDir, outDir);

fprintf(['Observed snapshots real_191-41.pdf and real_191-42.pdf require ' ...
    'the source x-y-u-v radar fields and were not overwritten.\n']);
fprintf('[make_overleaf_figures] Finished writing %s\n', outDir);
end

function makeConstantIntensity(rootDir, dataName, data, outDir)
catalogFile = fullfile(rootDir, 'exp_kernel', 'results', dataName, ...
    'validation', 'sim_catalog.mat');
loaded = load(catalogFile, 'catalog');
catalog = loaded.catalog;
params = catalog.summary.params;

[tVolterra, qVolterra] = volterra_solver( ...
    fullfile(rootDir, 'data', [dataName '.txt']), ...
    params.lambda0, params.K0, params.omega, params.rho, ...
    catalog.summary.T, 'dt', 0.25, 'Verbose', false);
volterra = struct('t_volterra', tVolterra, 'Q_volterra', qVolterra);

comparison = eddyplot.FigSimValidation(data, catalog, volterra);
saveFigure(comparison.drawIntensityComparison(), outDir, ...
    sprintf('sim_intensity_comparison-%s', dataName));
end

function makeLifetimeSurvival(datasets, outDir)
ps = eddyplot.PlotStyle;
colors = ps.colors();
fig = ps.newFig(ps.FIG_WIDTH_FULL * 1.15, ps.FIG_HEIGHT * 1.2);
ax = axes(fig);
hold(ax, 'on');

lineStyles = {'-', '--'};
lineColors = {colors.blue, colors.orange};
labels = {'First 14 Days', 'Last 14 Days'};
for k = 1:numel(datasets)
    lifetime = sort(double(datasets{k}.lifetime(:)));
    thresholds = unique(lifetime);
    survival = arrayfun(@(value) mean(lifetime > value), thresholds);
    stairs(ax, thresholds, max(survival, 1e-3), lineStyles{k}, ...
        'Color', lineColors{k}, 'LineWidth', 1.5, ...
        'DisplayName', labels{k});
end

ax.YScale = 'log';
ylim(ax, [1e-3, 1]);
xlim(ax, [0, 360]);
legend(ax, 'Location', 'northeast', 'Box', 'off');
ps.polishAxes(ax, 'Empirical Lifetime Survival $\bar{\gamma}(\tau)$', ...
    'Lifetime threshold $\tau$ (min)', '$P(L>\tau)$');
hold(ax, 'off');
saveFigure(fig, outDir, 'lifetime_survival_comparison');
end

function copyTimeVaryingFigures(rootDir, outDir)
specs = {
    84, 'paper_timevarying_background_volterra_bins84.pdf', ...
        'timevarying_background_volterra_4h.pdf';
    84, 'paper_timevarying_arrivals_observed_simulated_bins84.pdf', ...
        'timevarying_arrivals_observed_simulated_4h.pdf';
    42, 'paper_timevarying_background_volterra_bins42.pdf', ...
        'timevarying_background_volterra_8h.pdf';
    42, 'paper_timevarying_arrivals_observed_simulated_bins42.pdf', ...
        'timevarying_arrivals_observed_simulated_8h.pdf'
};

for k = 1:size(specs, 1)
    binCount = specs{k, 1};
    source = fullfile(rootDir, 'time_varying_intensity', ...
        sprintf('results_bins%d', binCount), 'figures', specs{k, 2});
    destination = fullfile(outDir, specs{k, 3});
    if isfile(source)
        copyfile(source, destination);
        fprintf('  Copied: %s\n', destination);
    else
        warning('make_overleaf_figures:timeVaryingResultMissing', ...
            ['Retained existing %s. To rebuild it with authorized data, ' ...
             'run time_varying_intensity/run_timevarying_bin_analysis.m.'], ...
            specs{k, 3});
    end
end
end

function saveFigure(fig, outDir, baseName)
outputFile = fullfile(outDir, [baseName '.pdf']);
exportgraphics(fig, outputFile, 'ContentType', 'vector');
close(fig);
fprintf('  Saved: %s\n', outputFile);
end
