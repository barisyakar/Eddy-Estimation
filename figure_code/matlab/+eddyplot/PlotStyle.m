classdef PlotStyle
% PLOTSTYLE  Centralised style constants for all eddy figures.
%
%   All methods are Static:
%       eddyplot.PlotStyle.apply()       % set global defaults
%       c = eddyplot.PlotStyle.colors()  % colorblind-safe palette

    properties (Constant)
        FONT_SIZE        = 10
        TITLE_SIZE       = 11
        LABEL_SIZE       = 10
        TICK_SIZE        = 9
        LEGEND_SIZE      = 9
        LINE_WIDTH       = 1.4
        MARKER_SIZE      = 6
        STEM_LINE_WIDTH  = 0.9
        FIG_WIDTH_FULL   = 6.3    % inches (~16 cm, full text width)
        FIG_WIDTH_HALF   = 3.1    % inches (~8 cm, half column)
        FIG_HEIGHT       = 3.2
        EXPORT_DPI       = 300
    end

    methods (Static)
        %------------------------------------------------------------------
        function apply()
            % Set global MATLAB defaults for publication-quality plots with
            % LaTeX-rendered text.
            set(groot, 'defaultAxesFontSize',               eddyplot.PlotStyle.FONT_SIZE);
            set(groot, 'defaultAxesFontName',               'Times New Roman');
            set(groot, 'defaultAxesLineWidth',              0.8);
            set(groot, 'defaultAxesBox',                    'off');
            set(groot, 'defaultAxesXGrid',                  'on');
            set(groot, 'defaultAxesYGrid',                  'on');
            set(groot, 'defaultAxesGridAlpha',              0.25);
            set(groot, 'defaultAxesGridLineStyle',          '--');
            set(groot, 'defaultAxesTickLabelInterpreter',   'latex');
            set(groot, 'defaultLineLineWidth',              eddyplot.PlotStyle.LINE_WIDTH);
            set(groot, 'defaultLineMarkerSize',             eddyplot.PlotStyle.MARKER_SIZE);
            set(groot, 'defaultTextFontName',               'Times New Roman');
            set(groot, 'defaultTextFontSize',               eddyplot.PlotStyle.FONT_SIZE);
            set(groot, 'defaultTextInterpreter',            'latex');
            set(groot, 'defaultLegendInterpreter',          'latex');
            set(groot, 'defaultColorbarTickLabelInterpreter','latex');
            set(groot, 'defaultFigureColor',                'w');
        end

        %------------------------------------------------------------------
        function c = colors()
            c.blue   = [0  114 178] / 255;
            c.orange = [230 159   0] / 255;
            c.green  = [  0 158 115] / 255;
            c.red    = [213  94   0] / 255;
            c.purple = [204 121 167] / 255;
            c.sky    = [ 86 180 233] / 255;
            c.gray   = [0.55 0.55 0.55];
        end

        %------------------------------------------------------------------
        function fig = newFig(width, height)
            arguments
                width  (1,1) double = eddyplot.PlotStyle.FIG_WIDTH_FULL
                height (1,1) double = eddyplot.PlotStyle.FIG_HEIGHT
            end
            fig = figure('Units','inches','Position',[1 1 width height], ...
                         'Color','w','PaperPositionMode','auto');
        end

        %------------------------------------------------------------------
        function polishAxes(ax, titleStr, xlabelStr, ylabelStr)
            arguments
                ax         matlab.graphics.axis.Axes
                titleStr   (1,1) string = ""
                xlabelStr  (1,1) string = ""
                ylabelStr  (1,1) string = ""
            end
            ps = eddyplot.PlotStyle;
            if titleStr  ~= "", title(ax,  titleStr,  'FontSize', ps.TITLE_SIZE, 'FontWeight','bold'); end
            if xlabelStr ~= "", xlabel(ax, xlabelStr, 'FontSize', ps.LABEL_SIZE); end
            if ylabelStr ~= "", ylabel(ax, ylabelStr, 'FontSize', ps.LABEL_SIZE); end
            ax.FontSize = ps.FONT_SIZE;
            ax.FontName = 'Times New Roman';
            ax.TickDir  = 'out';
            box(ax, 'off');
        end

        %------------------------------------------------------------------
        function saveFig(fig, outDir, baseName)
            arguments
                fig     matlab.ui.Figure
                outDir  (1,1) string
                baseName(1,1) string
            end
            if ~isfolder(outDir), mkdir(outDir); end
            pdfPath = fullfile(outDir, baseName + ".pdf");
            exportgraphics(fig, pdfPath, 'ContentType','vector');
            fprintf('  Saved: %s\n', pdfPath);
        end
    end
end
