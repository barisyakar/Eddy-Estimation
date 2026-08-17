function bgMap = learn_spatial_background(dataFilePath, varargin)
%LEARN_SPATIAL_BACKGROUND  Learn the spatial density of eddy formations
%   from empirical data using 2D Gaussian kernel density estimation.
%
%   bgMap = learn_spatial_background(dataFilePath)
%   bgMap = learn_spatial_background(..., 'Name', Value)
%
%   This function reads the empirical eddy catalog, computes a 2D Gaussian
%   KDE of the spatial locations on a discrete grid, and returns a
%   normalized probability map. This map is used in simulation to sample
%   background event locations from the empirical spatial distribution
%   instead of from a uniform distribution.
%
%   The estimation model retains a scalar lambda_0 (average background rate
%   per unit time per unit area); this offline spatial learning is used ONLY
%   in simulation to produce spatially realistic synthetic catalogs.
%
%   Inputs:
%     dataFilePath  - path to a Lifetime-*.txt data file
%
%   Name-Value Options:
%     'GridSide'    - side length of the square spatial domain (default 91)
%     'Bandwidth'   - KDE bandwidth: 'silverman' (default) or a positive scalar
%     'OutputFile'  - path to save the bgMap struct as .mat (default: no save)
%     'Verbose'     - print progress (default true)
%
%   Output:
%     bgMap - struct with fields:
%       .density     - GridSide x GridSide normalized density matrix
%       .prob        - GridSide x GridSide probability matrix (sums to 1)
%       .cumProb     - flattened cumulative probability for sampling
%       .gridX       - x grid centers (0:GridSide-1)
%       .gridY       - y grid centers (0:GridSide-1)
%       .bandwidth   - bandwidth used
%       .N_data      - number of data points used
%       .GridSide    - grid side length
%
%   Example:
%     bgMap = learn_spatial_background('../data/Lifetime-first14days.txt', ...
%                                      'OutputFile', 'results/spatial_bg.mat');
%
% See also: simulation, ExpKernelEstimation

%% ---- Parse inputs -------------------------------------------------------
p = inputParser;
addRequired(p, 'dataFilePath', @ischar);
addParameter(p, 'GridSide', 91, @(v) isfinite(v) && v > 0);
addParameter(p, 'Bandwidth', 'silverman');
addParameter(p, 'OutputFile', '', @ischar);
addParameter(p, 'Verbose', true, @islogical);
parse(p, dataFilePath, varargin{:});

GridSide   = p.Results.GridSide;
bwInput    = p.Results.Bandwidth;
outputFile = p.Results.OutputFile;
verbose    = p.Results.Verbose;

%% ---- Load empirical data ------------------------------------------------
if verbose; fprintf('[learn_spatial_background] Loading data from: %s\n', dataFilePath); end

[~, ~, firstx, firsty, ~, ~, ~, ~] = ...
    textread(dataFilePath, '%s%d%d%d%f%f%d%d');  %#ok<DTXTRD>

x_data = double(firstx);
y_data = double(firsty);
N_data = length(x_data);

if verbose; fprintf('  Loaded %d event locations\n', N_data); end

%% ---- Compute bandwidth --------------------------------------------------
if ischar(bwInput) && strcmpi(bwInput, 'silverman')
    % Silverman's rule of thumb for 2D: h = N^(-1/6) * sigma_pooled
    sigma_x = std(x_data);
    sigma_y = std(y_data);
    sigma_pooled = sqrt(0.5 * (sigma_x^2 + sigma_y^2));
    h = N_data^(-1/6) * sigma_pooled;
    if verbose
        fprintf('  Silverman bandwidth: h = %.3f (sigma_x=%.2f, sigma_y=%.2f)\n', ...
                h, sigma_x, sigma_y);
    end
elseif isnumeric(bwInput) && isscalar(bwInput) && bwInput > 0
    h = bwInput;
    if verbose; fprintf('  User-specified bandwidth: h = %.3f\n', h); end
else
    error('learn_spatial_background:invalidBandwidth', ...
          'Bandwidth must be ''silverman'' or a positive scalar.');
end

%% ---- Compute 2D KDE on grid --------------------------------------------
% Grid centers at integer coordinates 0, 1, ..., GridSide-1
% (matching the coordinate system used in ExpKernelEstimation and simulation)
gx = (0:GridSide-1)';
gy = (0:GridSide-1)';
[GX, GY] = meshgrid(gx, gy);

% Evaluate KDE: f(x,y) = (1/N) sum_i (1/(2*pi*h^2)) exp(-((x-xi)^2+(y-yi)^2)/(2*h^2))
density = zeros(GridSide, GridSide);
coeff = 1 / (2 * pi * h^2);

if verbose; fprintf('  Computing 2D KDE on %d x %d grid...\n', GridSide, GridSide); end

for i = 1:N_data
    dx2 = (GX - x_data(i)).^2;
    dy2 = (GY - y_data(i)).^2;
    density = density + coeff * exp(-(dx2 + dy2) / (2 * h^2));
end
density = density / N_data;

%% ---- Normalize to probability map --------------------------------------
% Each grid cell is 1x1, so probability ~ density (up to normalization)
prob = density / sum(density(:));  % Normalize to sum to 1

% Compute cumulative probability for efficient inverse-CDF sampling
% (column-major linearization: prob(:) reads column by column)
cumProb = cumsum(prob(:));
cumProb(end) = 1.0;  % Ensure exact normalization

%% ---- Build output struct ------------------------------------------------
bgMap = struct();
bgMap.density   = density;
bgMap.prob      = prob;
bgMap.cumProb   = cumProb;
bgMap.gridX     = gx;
bgMap.gridY     = gy;
bgMap.bandwidth = h;
bgMap.N_data    = N_data;
bgMap.GridSide  = GridSide;

%% ---- Summary ------------------------------------------------------------
if verbose
    fprintf('\n=== SPATIAL BACKGROUND MAP ===\n');
    fprintf('  Grid: %d x %d\n', GridSide, GridSide);
    fprintf('  Bandwidth: h = %.3f grid units\n', h);
    fprintf('  Peak density: %.6e\n', max(density(:)));
    fprintf('  Min density:  %.6e\n', min(density(:)));
    
    % Find peak location
    [~, peakIdx] = max(prob(:));
    [peakRow, peakCol] = ind2sub([GridSide, GridSide], peakIdx);
    fprintf('  Peak location: (x=%d, y=%d)\n', gx(peakCol), gy(peakRow));
    
    % Effective support (cells with >1% of peak probability)
    threshold = 0.01 * max(prob(:));
    nActive = sum(prob(:) > threshold);
    fprintf('  Active cells (>1%% of peak): %d / %d (%.1f%%)\n', ...
            nActive, GridSide^2, 100*nActive/GridSide^2);
    fprintf('==============================\n\n');
end

%% ---- Save ---------------------------------------------------------------
if ~isempty(outputFile)
    if verbose; fprintf('[learn_spatial_background] Saving to: %s\n', outputFile); end
    save(outputFile, 'bgMap');
end

end
