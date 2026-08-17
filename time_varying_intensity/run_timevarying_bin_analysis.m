% run_timevarying_bin_analysis.m
% =========================================================================
% Fit time-varying baseline Hawkes models for selected bin counts, simulate
% one representative catalog, and create paper-ready comparison figures.
%
% Paper bin counts:
%   84 bins = 4-hour bins over 14 days
%   42 bins = 8-hour bins over 14 days
%
% Outputs for each B:
%   time_varying_intensity/results_binsB/
%       Lifetime-first14days_binsB.mat
%       Lifetime-last14days_binsB.mat
%       paper_summary_binsB.csv
%       paper_baseline_timeseries_binsB.csv
%       simulation/<dataset>/sim_catalog_seed42.mat
%       figures/paper_timevarying_background_volterra_binsB.{pdf,png}
%       figures/paper_timevarying_arrivals_observed_simulated_binsB.{pdf,png}
%
% The two figure PDFs are also copied into figures/ with names that
% encode the physical bin width:
%   timevarying_background_volterra_4h.pdf
%   timevarying_arrivals_observed_simulated_4h.pdf
%   timevarying_background_volterra_8h.pdf
%   timevarying_arrivals_observed_simulated_8h.pdf
% =========================================================================

clear; close all; clc;

SCRIPT_DIR = fileparts(mfilename('fullpath'));
ROOT_DIR = fileparts(SCRIPT_DIR);
DATA_DIR = fullfile(ROOT_DIR, 'data');
PAPER_FIG_DIR = fullfile(ROOT_DIR, 'figures');

BIN_COUNTS = [84 42];
REP_SEED = 42;
GRID_SIDE = 91;
GRID_AREA = GRID_SIDE^2;
WINDOW_STEPS = 14 * 24 * 4;
VOLTERRA_DT = 0.25;

DATASETS = {
    'Lifetime-first14days', 'First 14 days', fullfile(DATA_DIR, 'Lifetime-first14days.txt');
    'Lifetime-last14days',  'Last 14 days',  fullfile(DATA_DIR, 'Lifetime-last14days.txt')
};

colors = struct();
colors.blue = [0.000, 0.447, 0.741];
colors.orange = [0.850, 0.325, 0.098];
colors.green = [0.466, 0.674, 0.188];
colors.gray = [0.420, 0.420, 0.420];
colors.dark = [0.100, 0.100, 0.100];

if ~isfolder(PAPER_FIG_DIR); mkdir(PAPER_FIG_DIR); end

for bIdx = 1:numel(BIN_COUNTS)
    B = BIN_COUNTS(bIdx);
    binHours = (WINDOW_STEPS / B) * 0.25;
    resultsDir = fullfile(SCRIPT_DIR, sprintf('results_bins%d', B));
    simDir = fullfile(resultsDir, 'simulation');
    figDir = fullfile(resultsDir, 'figures');
    if ~isfolder(resultsDir); mkdir(resultsDir); end
    if ~isfolder(simDir); mkdir(simDir); end
    if ~isfolder(figDir); mkdir(figDir); end

    fprintf('\n========================================================\n');
    fprintf('[run_timevarying_bin_analysis] B=%d (%.1f-hour bins)\n', B, binHours);
    fprintf('========================================================\n');

    results = struct([]);
    for d = 1:size(DATASETS, 1)
        dataName = DATASETS{d, 1};
        displayLabel = DATASETS{d, 2};
        dataFile = DATASETS{d, 3};
        fitFile = fullfile(resultsDir, sprintf('%s_bins%d.mat', dataName, B));

        fprintf('\n--- %s, B=%d ---\n', dataName, B);
        fitWasUpdated = false;
        fitIsCurrent = false;
        if isfile(fitFile)
            fitMeta = load(fitFile, 'T', 'eventTimeBlocks');
            fitIsCurrent = isfield(fitMeta, 'T') && fitMeta.T == WINDOW_STEPS && ...
                           isfield(fitMeta, 'eventTimeBlocks');
        end
        if fitIsCurrent
            fprintf('  Fit exists, loading: %s\n', fitFile);
        else
            [initLambda0, initK0, initOmega, initRho] = initialValues(dataFile, dataName, SCRIPT_DIR, GRID_AREA, WINDOW_STEPS);
            fprintf('  Initial values: lambda0=%.6e, K0=%.6e, omega=%.4f, rho=%.6f\n', ...
                initLambda0, initK0, initOmega, initRho);
            ExpKernelEstimation(dataFile, initLambda0, initK0, initOmega, initRho, ...
                fitFile, B, WINDOW_STEPS);
            fitWasUpdated = true;
        end

        fit = load(fitFile, 'final_lambda0', 'final_K0', 'final_omega', ...
            'final_rho', 'T', 'GridArea', 'n', 'G_sum', 'Pr0');
        lambda0 = fit.final_lambda0(:);
        K0 = fit.final_K0;
        omega = fit.final_omega;
        rho = fit.final_rho;
        T = fit.T;
        binDuration = T / B;
        binEdges = linspace(0, T, B + 1);
        binCenters = 0.5 * (binEdges(1:end-1) + binEdges(2:end));
        binCentersDay = binCenters(:) / 96;
        baselineRate = lambda0 * fit.GridArea;

        outSimDir = fullfile(simDir, dataName);
        if ~isfolder(outSimDir); mkdir(outSimDir); end

        spatialBgFile = fullfile(outSimDir, 'spatial_bg.mat');
        if ~isfile(spatialBgFile)
            refSpatialBg = fullfile(SCRIPT_DIR, 'results_bins84', 'simulation', dataName, 'spatial_bg.mat');
            if isfile(refSpatialBg)
                copyfile(refSpatialBg, spatialBgFile);
                fprintf('  Reused spatial background from bins84.\n');
            else
                learn_spatial_background(dataFile, 'OutputFile', spatialBgFile, 'Verbose', true);
            end
        end
        tmpBg = load(spatialBgFile, 'bgMap');

        simFile = fullfile(outSimDir, sprintf('sim_catalog_seed%d.mat', REP_SEED));
        if isfile(simFile) && ~fitWasUpdated
            tmpSim = load(simFile, 'catalog');
            simIsCurrent = isfield(tmpSim.catalog, 'summary') && ...
                           isfield(tmpSim.catalog.summary, 'T') && ...
                           tmpSim.catalog.summary.T == T;
        else
            simIsCurrent = false;
        end
        if simIsCurrent
            catalog = tmpSim.catalog;
            fprintf('  Representative simulation exists, loading seed=%d.\n', REP_SEED);
        else
            catalog = simulation(dataFile, lambda0, K0, omega, rho, T, ...
                'Seed', REP_SEED, 'SpatialBgMap', tmpBg.bgMap, ...
                'OutputFile', simFile, 'Verbose', true);
        end

        real = loadRealCatalog(dataFile);
        realCounts = countByEdges(real.timeStep, binEdges);
        simCounts = countByEdges(catalog.time, binEdges);
        realRate = realCounts(:) / binDuration;
        simRate = simCounts(:) / binDuration;

        [tVolterra, QVolterra] = volterra_solver(dataFile, lambda0, K0, omega, rho, ...
            T, 'dt', VOLTERRA_DT, 'Verbose', false);
        volterraAtBins = interp1(tVolterra, QVolterra, binCenters(:), 'linear', 'extrap');
        results(d).dataName = dataName;
        results(d).displayLabel = displayLabel;
        results(d).B = B;
        results(d).binHours = binHours;
        results(d).T = T;
        results(d).binDuration = binDuration;
        results(d).binEdgesDay = binEdges(:) / 96;
        results(d).binCentersTime = binCenters(:);
        results(d).binCentersDay = binCentersDay;
        results(d).baselineRate = baselineRate;
        results(d).volterraRate = volterraAtBins;
        results(d).realRate = realRate;
        results(d).simRate = simRate;
        results(d).nReal = fit.n;
        results(d).nSim = catalog.summary.N_total;
        results(d).nSimBackground = catalog.summary.N_background;
        results(d).nSimOffspring = catalog.summary.N_offspring;
        results(d).fitTriggeredPct = 100 * (1 - mean(fit.Pr0));
        results(d).fitBranching = K0 * fit.G_sum / fit.n;
        results(d).K0 = K0;
        results(d).omega = omega;
        results(d).rho = rho;
        results(d).lambda0Mean = mean(lambda0);
        results(d).baselineMean = mean(baselineRate);
        results(d).baselineMin = min(baselineRate);
        results(d).baselineMax = max(baselineRate);
        results(d).baselineCV = std(baselineRate) / mean(baselineRate);

        fprintf('  Real N=%d, sim N=%d, branching=%.4f, EM triggered=%.2f%%\n', ...
            fit.n, catalog.summary.N_total, results(d).fitBranching, results(d).fitTriggeredPct);
    end

    backgroundFig = plotBackgroundVolterra(results, figDir, colors);
    arrivalsFig = plotObservedSimulated(results, figDir, colors);
    copyPaperFigures(backgroundFig, arrivalsFig, PAPER_FIG_DIR, binHours);
    writeSummaryFiles(results, resultsDir);
end

fprintf('\n[run_timevarying_bin_analysis] Done.\n\n');

%% ---- Helpers ------------------------------------------------------------

function [lambda0, K0, omega, rho] = initialValues(dataFile, dataName, scriptDir, gridArea, windowSteps)
    refFile = fullfile(scriptDir, 'results_bins84', sprintf('%s_bins84.mat', dataName));
    if isfile(refFile)
        ref = load(refFile, 'final_lambda0', 'final_K0', 'final_omega', 'final_rho');
        lambda0 = mean(ref.final_lambda0(:));
        K0 = ref.final_K0;
        omega = ref.final_omega;
        rho = ref.final_rho;
        return;
    end

    kappa = computeKappa(dataFile, gridArea, windowSteps);
    lambda0 = 0.5 * kappa;
    K0 = 1;
    omega = 0.5;
    rho = 0.1;
end

function kappa = computeKappa(dataFile, gridArea, windowSteps)
    fid = fopen(dataFile, 'r');
    if fid < 0
        error('run_timevarying_bin_analysis:fileNotFound', 'Cannot open: %s', dataFile);
    end
    raw = textscan(fid, '%s %d %d %d %f %f %d %d');
    fclose(fid);
    nTotal = numel(raw{1});
    kappa = nTotal / (gridArea * windowSteps);
end

function real = loadRealCatalog(dataFile)
    [s, ~, ~, ~, ~, ~, ~, ~] = textread(dataFile, '%s%d%d%d%f%f%d%d'); %#ok<DTXTRD>
    real.timeStep = timestamp_blocks(s);
end

function counts = countByEdges(values, edges)
    counts = histcounts(values, edges);
    if numel(values) > 0
        counts(end) = counts(end) + sum(values == edges(end));
    end
end

function backgroundPdf = plotBackgroundVolterra(results, figDir, colors)
    B = results(1).B;
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 920 430]);
    tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for d = 1:numel(results)
        r = results(d);
        ax = nexttile;
        hold(ax, 'on');
        stairs(ax, r.binEdgesDay, [r.baselineRate; r.baselineRate(end)], ...
            'Color', colors.orange, 'LineWidth', 1.8);
        plot(ax, r.binCentersDay, r.volterraRate, '--', ...
            'Color', colors.blue, 'LineWidth', 1.7);
        hold(ax, 'off');
        title(ax, sprintf('%s (%g-hour bins)', r.displayLabel, r.binHours));
        xlabel(ax, 'Day within period');
        ylabel(ax, 'Events per 15-min step');
        grid(ax, 'on'); box(ax, 'off');
        legend(ax, {'Estimated background $\hat{\lambda}_0(t)|S|$', ...
            'Volterra $Q(t)$'}, ...
            'Interpreter', 'latex', 'Location', 'northwest', 'Box', 'off');
    end
    outBase = sprintf('paper_timevarying_background_volterra_bins%d', B);
    backgroundPdf = saveFigure(fig, figDir, outBase);
end

function arrivalsPdf = plotObservedSimulated(results, figDir, colors)
    B = results(1).B;
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 460]);
    tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'loose');
    for d = 1:numel(results)
        r = results(d);
        ax = nexttile;
        hold(ax, 'on');
        plot(ax, r.binCentersTime, r.simRate, 'Color', colors.blue, 'LineWidth', 1.05);
        plot(ax, r.binCentersTime, r.volterraRate, '--', ...
            'Color', colors.orange, 'LineWidth', 1.7);
        plot(ax, r.binCentersTime, r.realRate, '--', ...
            'Color', [0.72 0.72 0.72], 'LineWidth', 1.05);
        hold(ax, 'off');
        title(ax, sprintf('Intensity Comparison -- %s (%g-hour bins)', ...
            r.displayLabel, r.binHours), 'FontName', 'Times', 'FontWeight', 'normal');
        xlabel(ax, 'Time (time steps)', 'FontName', 'Times');
        ylabel(ax, 'Event rate per time step', 'FontName', 'Times');
        xlim(ax, [0 r.T]);
        ymax = max([r.simRate(:); r.volterraRate(:); r.realRate(:)]);
        ylim(ax, [0 max(1, 1.08 * ymax)]);
        grid(ax, 'on');
        ax.GridAlpha = 0.25;
        ax.MinorGridAlpha = 0.10;
        ax.XMinorTick = 'on';
        ax.YMinorTick = 'on';
        ax.TickDir = 'both';
        ax.Box = 'on';
        ax.FontName = 'Times';
        legend(ax, {'Simulation', ...
            'Volterra $Q(t)$', 'Real data'}, ...
            'Interpreter', 'latex', 'Location', 'northeast', 'Box', 'off');
    end
    outBase = sprintf('paper_timevarying_arrivals_observed_simulated_bins%d', B);
    arrivalsPdf = saveFigure(fig, figDir, outBase);
end

function pdfFile = saveFigure(fig, outDir, name)
    pngFile = fullfile(outDir, [name '.png']);
    pdfFile = fullfile(outDir, [name '.pdf']);
    exportgraphics(fig, pngFile, 'Resolution', 220, 'Padding', 'figure');
    exportgraphics(fig, pdfFile, 'ContentType', 'vector', 'Padding', 'figure');
    close(fig);
    fprintf('  Saved %s\n', pdfFile);
end

function copyPaperFigures(backgroundPdf, arrivalsPdf, paperFigDir, binHours)
    suffix = sprintf('%gh', binHours);
    suffix = strrep(suffix, '.', 'p');
    copyfile(backgroundPdf, fullfile(paperFigDir, sprintf('timevarying_background_volterra_%s.pdf', suffix)));
    copyfile(arrivalsPdf, fullfile(paperFigDir, sprintf('timevarying_arrivals_observed_simulated_%s.pdf', suffix)));
end

function writeSummaryFiles(results, resultsDir)
    B = results(1).B;
    summaryRows = cell(numel(results), 1);
    baselineRows = cell(numel(results), 1);
    for d = 1:numel(results)
        r = results(d);
        summaryRows{d} = table(string(r.dataName), r.B, r.binHours, r.nReal, r.nSim, ...
            r.nSimBackground, r.nSimOffspring, r.lambda0Mean, r.baselineMean, ...
            r.baselineMin, r.baselineMax, r.baselineCV, r.fitBranching, ...
            r.fitTriggeredPct, r.K0, r.omega, r.rho, ...
            'VariableNames', {'dataset', 'num_bins', 'bin_hours', 'real_events', ...
            'sim_events_seed42', 'sim_background_seed42', 'sim_offspring_seed42', ...
            'lambda0_mean', 'baseline_mean_events_per_step', ...
            'baseline_min_events_per_step', 'baseline_max_events_per_step', ...
            'baseline_cv', 'fit_branching_ratio', 'fit_triggered_pct', ...
            'K0', 'omega', 'rho'});

        nBins = numel(r.baselineRate);
        baselineRows{d} = table(repmat(string(r.dataName), nBins, 1), ...
            (1:nBins)', r.binEdgesDay(1:end-1), r.binEdgesDay(2:end), ...
            r.baselineRate, r.volterraRate, r.realRate, r.simRate, ...
            'VariableNames', {'dataset', 'bin', 'bin_start_day', 'bin_end_day', ...
            'background_events_per_step', 'volterra_events_per_step', ...
            'observed_events_per_step', 'simulated_events_per_step'});
    end

    writetable(vertcat(summaryRows{:}), fullfile(resultsDir, sprintf('paper_summary_bins%d.csv', B)));
    writetable(vertcat(baselineRows{:}), fullfile(resultsDir, sprintf('paper_baseline_timeseries_bins%d.csv', B)));
end
