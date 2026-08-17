classdef FigSpatialDensity
% FIGSPATIALDENSITY  Spatial density of eddy events.
%
%   Two views:
%     2D contour (Gaussian KDE)  — view3d = false (default)
%     3D bar chart (histogram2)  — view3d = true  (matches paper Figure 3)
%
%   Usage:
%       d   = eddyplot.EddyData('...').load();
%       fig = eddyplot.FigSpatialDensity(d).draw();                % 2D contour
%       fig = eddyplot.FigSpatialDensity(d, 60, 0, true, 10).draw(); % 3D bars

    properties (SetAccess = private)
        data       eddyplot.EddyData
        nGrid      (1,1) double = 60
        bandwidth  (1,1) double = 0    % 0 = Scott's rule (2D only)
        view3d     (1,1) logical = false
        binSize3d  (1,1) double = 10   % spatial bin width for 3D bars
    end

    methods
        function obj = FigSpatialDensity(data, nGrid, bandwidth, view3d, binSize3d)
            arguments
                data      eddyplot.EddyData
                nGrid     (1,1) double = 60
                bandwidth (1,1) double = 0
                view3d    (1,1) logical = false
                binSize3d (1,1) double = 10
            end
            obj.data      = data;
            obj.nGrid     = nGrid;
            obj.bandwidth = bandwidth;
            obj.view3d    = view3d;
            obj.binSize3d = binSize3d;
        end

        function fig = draw(obj)
            if obj.view3d
                fig = obj.draw3d();
            else
                fig = obj.draw2d();
            end
        end
    end

    methods (Access = private)
        %------------------------------------------------------------------
        function fig = draw2d(obj)
            % Gaussian KDE contour with scatter overlay
            x = obj.data.x;
            y = obj.data.y;

            xv = linspace(min(x)-1, max(x)+1, obj.nGrid);
            yv = linspace(min(y)-1, max(y)+1, obj.nGrid);
            [Xg, Yg] = meshgrid(xv, yv);

            if obj.bandwidth == 0
                n  = numel(x);
                hx = n^(-1/6) * std(x);
                hy = n^(-1/6) * std(y);
            else
                hx = obj.bandwidth;
                hy = obj.bandwidth;
            end
            Z = eddyplot.FigSpatialDensity.gaussianKDE(x, y, Xg, Yg, hx, hy);

            ps  = eddyplot.PlotStyle;
            c   = ps.colors();
            fig = ps.newFig(ps.FIG_WIDTH_HALF * 1.5, ps.FIG_HEIGHT * 1.3);
            ax  = axes(fig);

            contourf(ax, Xg, Yg, Z, 20, 'LineColor','none');
            cmap = flipud(hot(256));
            cmap(1,:) = [0.95 0.95 0.95];
            colormap(ax, cmap);
            cb = colorbar(ax);
            cb.Label.String  = 'Density';
            cb.Label.FontSize = ps.LABEL_SIZE;

            hold(ax,'on');
            scatter(ax, x, y, 4, c.blue, 'filled', 'MarkerFaceAlpha', 0.25);
            hold(ax,'off');
            axis(ax,'equal','tight');

            ps.polishAxes(ax, sprintf('Spatial Density --- %s', obj.data.formalLabel), ...
                'X (grid units)', 'Y (grid units)');
        end

        %------------------------------------------------------------------
        function fig = draw3d(obj)
            % 3-D bar chart of event counts (matches paper Figure 3)
            x  = obj.data.x;
            y  = obj.data.y;
            bw = obj.binSize3d;

            x_min = floor(min(x));   x_max = ceil(max(x));
            y_min = floor(min(y));   y_max = ceil(max(y));
            x_edges = x_min : bw : x_max;
            y_edges = y_min : bw : y_max;

            ps  = eddyplot.PlotStyle;
            fig = ps.newFig(ps.FIG_WIDTH_FULL * 1.3, ps.FIG_HEIGHT * 1.8);
            ax  = axes(fig);

            histogram2(ax, x, y, x_edges, y_edges, ...
                       'FaceColor', 'flat', ...
                       'EdgeColor', 'k', ...
                       'LineWidth', 0.25, ...
                       'FaceAlpha', 0.95);

            colormap(ax, parula);
            cb = colorbar(ax);
            cb.Label.String  = 'Event Count';
            cb.Label.FontSize = ps.LABEL_SIZE;

            view(ax, -55, 30);
            xlabel(ax, 'X Coordinate', 'FontSize', ps.LABEL_SIZE, 'Interpreter', 'latex');
            ylabel(ax, 'Y Coordinate', 'FontSize', ps.LABEL_SIZE, 'Interpreter', 'latex');
            zlabel(ax, 'Count',         'FontSize', ps.LABEL_SIZE, 'Interpreter', 'latex');
            title(ax, sprintf('Spatial Density of Events --- %s', obj.data.formalLabel), ...
                  'FontSize', ps.TITLE_SIZE, 'FontWeight', 'bold', 'Interpreter', 'latex');
            grid(ax, 'on');
        end
    end

    methods (Static, Access = private)
        function Z = gaussianKDE(x, y, Xg, Yg, hx, hy)
            N = numel(x);
            Z = zeros(size(Xg));
            for k = 1:N
                Z = Z + exp(-0.5*((Xg-x(k)).^2/hx^2 + (Yg-y(k)).^2/hy^2));
            end
            Z = Z / (N * 2*pi * hx * hy);
        end
    end
end
