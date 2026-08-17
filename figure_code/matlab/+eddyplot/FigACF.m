classdef FigACF
% FIGACF  ACF stem plot of binned eddy event counts (15-min bins).

    properties (SetAccess = private)
        data      eddyplot.EddyData
        binWidth  (1,1) double = 15
        nLags     (1,1) double = 48
    end

    methods
        function obj = FigACF(data, binWidth, nLags)
            arguments
                data      eddyplot.EddyData
                binWidth  (1,1) double = 15
                nLags     (1,1) double = 48
            end
            obj.data     = data;
            obj.binWidth = binWidth;
            obj.nLags    = nLags;
        end

        function fig = draw(obj)
            [counts, ~]  = obj.binnedCounts();
            [acfVals, lagHours] = obj.computeACF(counts);

            ps  = eddyplot.PlotStyle;
            c   = ps.colors();
            fig = ps.newFig(ps.FIG_WIDTH_FULL, ps.FIG_HEIGHT);
            ax  = axes(fig);

            stem(ax, lagHours, acfVals, ...
                 'Color', c.blue, 'MarkerFaceColor', c.blue, ...
                 'MarkerEdgeColor', c.blue, 'MarkerSize', 3, ...
                 'LineWidth', ps.STEM_LINE_WIDTH);

            ps.polishAxes(ax, ...
                sprintf('ACF of Event Counts (%d-min bins) --- %s', obj.binWidth, obj.data.formalLabel), ...
                'Lag (hours)', 'ACF');
            ylim(ax, [-0.3, 1.05]);
            xlim(ax, [0, lagHours(end)]);
        end
    end

    methods (Access = private)
        function [counts, edges] = binnedCounts(obj)
            t = obj.data.time_min;
            edges = 0:obj.binWidth:obj.data.totalDuration();
            counts = histcounts(t, edges);
        end

        function [acfVals, lagHours] = computeACF(obj, counts)
            x   = double(counts(:)) - mean(double(counts));
            r   = xcorr(x, obj.nLags, 'unbiased');
            mid = obj.nLags + 1;
            acov = r(mid:end);
            acfVals  = acov / acov(1);
            lagHours = (0:obj.nLags) * (obj.binWidth / 60);
        end
    end
end
