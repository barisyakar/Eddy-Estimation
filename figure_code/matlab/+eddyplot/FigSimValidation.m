classdef FigSimValidation
    % FigSimValidation  Comparison figures: real data vs. simulated data.
    %
    %   Methods:
    %     drawACFComparison       – side-by-side ACF stem plots
    %     drawSpatialComparison   – side-by-side spatial heatmaps
    %     drawIntensityComparison – empirical Q(t) vs Volterra solution
    %     drawParamRecovery       – bar chart of true vs recovered parameters

    properties (SetAccess = private)
        realData      eddyplot.EddyData
        simCatalog    struct
        volterra      struct
        recovery      struct
    end

    methods
        %------------------------------------------------------------------
        function obj = FigSimValidation(realData, simCatalog, volterra, recovery)
            arguments
                realData      eddyplot.EddyData
                simCatalog    struct
                volterra      struct = struct()
                recovery      struct = struct()
            end
            obj.realData   = realData;
            obj.simCatalog = simCatalog;
            obj.volterra   = volterra;
            obj.recovery   = recovery;
        end

        %------------------------------------------------------------------
        function fig = drawACFComparison(obj, maxLag)
            arguments
                obj
                maxLag (1,1) double = 40
            end
            ps = eddyplot.PlotStyle;
            c = eddyplot.PlotStyle.colors();
            fig = ps.newFig(ps.FIG_WIDTH_FULL * 2, ps.FIG_HEIGHT * 1.2);

            % Real ACF
            ax1 = subplot(1, 2, 1, 'Parent', fig);
            windowSteps = obj.realData.totalDuration() / 15;
            realCounts = obj.binCounts(obj.realData.time_min / 15, 1.0, windowSteps);
            realACF = obj.computeACF(realCounts, maxLag);
            stem(ax1, 1:maxLag, realACF(2:end), 'filled', 'Color', c.blue);
            hold(ax1, 'on');
            yline(ax1, 0, '--', 'Color', [0.5 0.5 0.5]);
            hold(ax1, 'off');
            ps.polishAxes(ax1, ...
                sprintf('ACF --- Real Data (%s)', obj.realData.formalLabel), ...
                'Lag (time steps)', 'Autocorrelation');

            % Simulated ACF
            ax2 = subplot(1, 2, 2, 'Parent', fig);
            simTime = obj.simCatalog.time(:);
            simCounts = obj.binCounts(simTime, 1.0, windowSteps);
            simACF = obj.computeACF(simCounts, maxLag);
            stem(ax2, 1:maxLag, simACF(2:end), 'filled', 'Color', c.red);
            hold(ax2, 'on');
            yline(ax2, 0, '--', 'Color', [0.5 0.5 0.5]);
            hold(ax2, 'off');
            ps.polishAxes(ax2, 'ACF --- Simulated Data', ...
                'Lag (time steps)', 'Autocorrelation');
        end

        %------------------------------------------------------------------
        function fig = drawSpatialComparison(obj, nBins)
            arguments
                obj
                nBins (1,1) double = 9
            end
            ps = eddyplot.PlotStyle;
            fig = ps.newFig(ps.FIG_WIDTH_FULL * 2.2, ps.FIG_HEIGHT * 1.2);

            % Real
            ax1 = subplot(1, 2, 1, 'Parent', fig);
            N_real = obj.realData.numEvents();
            histogram2(ax1, obj.realData.x, obj.realData.y, nBins, ...
                'DisplayStyle', 'tile', 'ShowEmptyBins', 'on');
            colormap(ax1, 'hot');
            colorbar(ax1);
            axis(ax1, 'equal');
            ps.polishAxes(ax1, ...
                sprintf('Spatial Density --- Real (N=%d)', N_real), 'X', 'Y');

            % Simulated
            ax2 = subplot(1, 2, 2, 'Parent', fig);
            simX = obj.simCatalog.x(:);
            simY = obj.simCatalog.y(:);
            N_sim = length(simX);
            histogram2(ax2, simX, simY, nBins, ...
                'DisplayStyle', 'tile', 'ShowEmptyBins', 'on');
            colormap(ax2, 'hot');
            colorbar(ax2);
            axis(ax2, 'equal');
            ps.polishAxes(ax2, ...
                sprintf('Spatial Density --- Simulated (N=%d)', N_sim), 'X', 'Y');
        end

        %------------------------------------------------------------------
        function fig = drawIntensityComparison(obj, binWidth)
            arguments
                obj
                binWidth (1,1) double = 5.0
            end
            ps = eddyplot.PlotStyle;
            c = eddyplot.PlotStyle.colors();
            fig = ps.newFig(ps.FIG_WIDTH_FULL * 1.5, ps.FIG_HEIGHT * 1.2);
            ax = axes(fig);
            hold(ax, 'on');

            % Simulated empirical rate
            simTime = obj.simCatalog.time(:);
            if ~isempty(simTime)
                Tmax = obj.realData.totalDuration() / 15;
                edges = 0:binWidth:(ceil(Tmax / binWidth) * binWidth);
                counts = histcounts(simTime, edges);
                centers = 0.5 * (edges(1:end-1) + edges(2:end));
                rate = counts / binWidth;
                plot(ax, centers, rate, 'Color', c.blue, 'LineWidth', 1.0, ...
                     'DisplayName', 'Simulated (empirical)');
            end

            % Volterra solution
            if isfield(obj.volterra, 't_volterra')
                tV = obj.volterra.t_volterra(:);
                QV = obj.volterra.Q_volterra(:);
                plot(ax, tV, QV, 'Color', c.red, 'LineWidth', 1.5, ...
                     'LineStyle', '--', 'DisplayName', 'Volterra Q(t)');
            end

            % Real data empirical rate
            realTimeSteps = obj.realData.time_min / 15;
            if ~isempty(realTimeSteps)
                TmaxR = obj.realData.totalDuration() / 15;
                edgesR = 0:binWidth:(ceil(TmaxR / binWidth) * binWidth);
                countsR = histcounts(realTimeSteps, edgesR);
                centersR = 0.5 * (edgesR(1:end-1) + edgesR(2:end));
                rateR = countsR / binWidth;
                plot(ax, centersR, rateR, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, ...
                     'DisplayName', 'Real data (empirical)');
            end

            legend(ax, 'Location', 'northeast', 'FontSize', ps.LEGEND_SIZE);
            ps.polishAxes(ax, ...
                sprintf('Intensity Comparison --- %s', obj.realData.formalLabel), ...
                'Time (time steps)', 'Event rate per time step');
            hold(ax, 'off');
        end

        %------------------------------------------------------------------
        function fig = drawParamRecovery(obj)
            ps = eddyplot.PlotStyle;
            c = eddyplot.PlotStyle.colors();
            fig = ps.newFig(ps.FIG_WIDTH_FULL, ps.FIG_HEIGHT * 1.2);
            ax = axes(fig);

            if ~isfield(obj.recovery, 'true_params')
                text(ax, 0.5, 0.5, 'No recovery data', ...
                     'HorizontalAlignment', 'center', 'FontSize', 14);
                return;
            end

            trueP = obj.recovery.true_params(:)';
            recP  = obj.recovery.recovered_params(:)';

            ratios = recP ./ max(abs(trueP), 1e-15);

            x = 1:4;
            w = 0.35;

            hold(ax, 'on');
            bar(ax, x - w/2, ones(1,4), w, 'FaceColor', c.blue, 'FaceAlpha', 0.8, ...
                'DisplayName', 'True (=1.0)');
            bar(ax, x + w/2, ratios, w, 'FaceColor', c.red, 'FaceAlpha', 0.8, ...
                'DisplayName', 'Recovered / True');

            % Value labels
            for i = 1:4
                text(ax, x(i) - w/2, 1.02, sprintf('%.2e', trueP(i)), ...
                     'HorizontalAlignment', 'center', 'FontSize', 6, 'Rotation', 45);
                text(ax, x(i) + w/2, ratios(i) + 0.02, sprintf('%.2e', recP(i)), ...
                     'HorizontalAlignment', 'center', 'FontSize', 6, 'Rotation', 45);
            end

            yline(ax, 1.0, '--', 'Color', [0.5 0.5 0.5], 'Alpha', 0.5);
            set(ax, 'XTick', x, 'XTickLabel', {'\lambda_0', 'K_0', '\omega', '\rho'});
            legend(ax, 'Location', 'northeast', 'FontSize', ps.LEGEND_SIZE);
            ps.polishAxes(ax, ...
                sprintf('Parameter Recovery --- %s', obj.realData.formalLabel), ...
                'Parameter', 'Ratio (Recovered / True)');
            hold(ax, 'off');
        end
    end

    methods (Static, Access = private)
        function counts = binCounts(times, binWidth, maxTime)
            if isempty(times)
                counts = 0;
                return;
            end
            if nargin < 3
                maxTime = max(times) + binWidth;
            end
            edges = 0:binWidth:(ceil(maxTime / binWidth) * binWidth);
            counts = histcounts(times, edges);
        end

        function acf = computeACF(series, maxLag)
            n = length(series);
            if n < 2
                acf = zeros(maxLag + 1, 1);
                return;
            end
            m = mean(series);
            v = var(series);
            if v < 1e-15
                acf = zeros(maxLag + 1, 1);
                return;
            end
            centered = series(:) - m;
            acf = zeros(maxLag + 1, 1);
            for lag = 0:min(maxLag, n-1)
                acf(lag+1) = sum(centered(1:n-lag) .* centered(lag+1:n)) / (n * v);
            end
        end
    end
end
