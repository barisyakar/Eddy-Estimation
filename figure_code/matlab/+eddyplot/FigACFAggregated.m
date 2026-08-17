classdef FigACFAggregated
% FIGACFAGGREGATED  Lag-1 ACF vs aggregation window size.

    properties (SetAccess = private)
        data        eddyplot.EddyData
        binWidth    (1,1) double = 15
        windowBins  (1,:) double = [1 2 4 8 16 32]
    end

    methods
        function obj = FigACFAggregated(data, binWidth, windowBins)
            arguments
                data        eddyplot.EddyData
                binWidth    (1,1) double  = 15
                windowBins  (1,:) double  = [1 2 4 8 16 32]
            end
            obj.data       = data;
            obj.binWidth   = binWidth;
            obj.windowBins = windowBins;
        end

        function fig = draw(obj)
            [counts, ~]              = obj.baseCounts();
            [windowMin, lagOneACF]   = obj.computeLagOne(counts);

            ps  = eddyplot.PlotStyle;
            c   = ps.colors();
            fig = ps.newFig(ps.FIG_WIDTH_HALF * 1.3, ps.FIG_HEIGHT);
            ax  = axes(fig);

            plot(ax, windowMin, lagOneACF, 'o-', ...
                 'Color', c.blue, 'MarkerFaceColor', c.blue, ...
                 'MarkerEdgeColor', c.blue, ...
                 'LineWidth', ps.LINE_WIDTH, 'MarkerSize', ps.MARKER_SIZE);

            ax.XScale = 'log';
            xticks(ax, windowMin);
            xticklabels(ax, arrayfun(@(v) sprintf('%g',v), windowMin, 'UniformOutput', false));
            ylim(ax, [0, 1]);

            ps.polishAxes(ax, ...
                sprintf('Lag-1 ACF vs Window --- %s', obj.data.formalLabel), ...
                'Aggregation window (minutes)', 'Lag-1 ACF');
        end
    end

    methods (Access = private)
        function [counts, edges] = baseCounts(obj)
            t = obj.data.time_min;
            edges = 0:obj.binWidth:obj.data.totalDuration();
            counts = histcounts(t, edges);
        end

        function [windowMin, lagOneACF] = computeLagOne(obj, counts)
            windowMin = zeros(1, numel(obj.windowBins));
            lagOneACF = zeros(1, numel(obj.windowBins));
            for k = 1:numel(obj.windowBins)
                w   = obj.windowBins(k);
                m   = floor(numel(counts)/w)*w;
                agg = sum(reshape(counts(1:m), w, []), 1);
                if numel(agg) > 2
                    R = corrcoef(agg(1:end-1), agg(2:end));
                    lagOneACF(k) = R(1,2);
                else
                    lagOneACF(k) = NaN;
                end
                windowMin(k) = w * obj.binWidth;
            end
            valid     = ~isnan(lagOneACF);
            windowMin = windowMin(valid);
            lagOneACF = lagOneACF(valid);
        end
    end
end
