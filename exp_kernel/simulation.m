function catalog = simulation(dataFilePath, lambda0, K0, omega, rho, T, varargin)
%SIMULATION  Simulate the spatio-temporal Hawkes eddy process via the
%            branching (cluster) representation (Algorithm 2 in the paper).
%
%   catalog = simulation(dataFilePath, lambda0, K0, omega, rho, T)
%   catalog = simulation(..., 'Name', Value)
%
%   This function implements the rejection-free cluster representation
%   described in Hawkes & Oakes (1974) and detailed in §Simulation of
%   the paper. The process is built generation-by-generation:
%
%     Generation 0 (background):
%       N_bg ~ Poisson(lambda0 * |S| * T)
%       Times are uniform on [0, T]. Locations follow the supplied empirical
%       spatial density, with a uniform distribution used only as a fallback.
%
%     Generation l+1 (offspring of generation l):
%       For each parent i with strain |a_i|/b_i and lifetime l_i:
%         N_i  ~ Poisson( K0 * (|a_i|/b_i) * (1 - exp(-omega*l_i)) )
%         delay tau ~ TruncatedExponential(rate=omega, bound=l_i)
%         location  ~ N(z_i, (1/(2*rho)) * I_2)
%
%     Marks (a, b, l) are bootstrapped from the empirical data file.
%
%   Inputs:
%     dataFilePath  – path to a Lifetime-*.txt data file (for mark bootstrap)
%     lambda0       – background rate per unit space per unit time
%     K0            – productivity parameter
%     omega         – temporal decay rate
%     rho           – spatial decay rate
%     T             – observation window length (in time-step units)
%
%   Name-Value Options:
%     'GridSide'    – side length of the square spatial domain (default 91)
%     'Seed'        – random seed for reproducibility (default 42)
%     'MaxGen'      – safety cap on generation depth (default 50)
%     'OutputFile'  – path to save catalog as .mat (default: no save)
%     'Verbose'     – print progress (default true)
%     'SpatialBgMap'– spatial-density struct from learn_spatial_background;
%                     empty selects the uniform fallback (default empty)
%
%   Output:
%     catalog – scalar struct containing arrays with fields:
%       .time     – event time (in time-step units)
%       .x, .y    – spatial coordinates
%       .aa       – amplitude (mark, already scaled as in estimation)
%       .bb       – radius    (mark, already scaled as in estimation)
%       .lifetime – lifetime  (mark, in time-step blocks)
%       .gen      – generation index (0 = background)
%       .parent   – parent index (0 = background)
%
%   Example (using estimated parameters from first 14 days):
%     cat = simulation('data/Lifetime-first14days.txt', ...
%                       0.00049833, 2.1154, 0.77657, 0.020512, 1344, ...
%                       'SpatialBgMap', bgMap);
%
% See also: ExpKernelEstimation, run_sweep

%% ---- Parse inputs -------------------------------------------------------
p = inputParser;
addRequired(p, 'dataFilePath', @ischar);
addRequired(p, 'lambda0',     @(v) isfinite(v) && v > 0);
addRequired(p, 'K0',          @(v) isfinite(v) && v > 0);
addRequired(p, 'omega',       @(v) isfinite(v) && v > 0);
addRequired(p, 'rho',         @(v) isfinite(v) && v > 0);
addRequired(p, 'T',           @(v) isfinite(v) && v > 0);
addParameter(p, 'GridSide', 91,   @(v) isfinite(v) && v > 0);
addParameter(p, 'Seed',    42,    @(v) isfinite(v) && v >= 0);
addParameter(p, 'MaxGen',  50,    @(v) isfinite(v) && v > 0);
addParameter(p, 'OutputFile', '', @ischar);
addParameter(p, 'Verbose', true,  @islogical);
addParameter(p, 'SpatialBgMap', [], @(v) isempty(v) || isstruct(v) || ischar(v));
parse(p, dataFilePath, lambda0, K0, omega, rho, T, varargin{:});

GridSide   = p.Results.GridSide;
GridArea   = GridSide^2;
seed       = p.Results.Seed;
maxGen     = p.Results.MaxGen;
outputFile = p.Results.OutputFile;
verbose    = p.Results.Verbose;
spatialBgMap = p.Results.SpatialBgMap;

% Load spatial background map if a file path is provided
if ischar(spatialBgMap) && ~isempty(spatialBgMap)
    if verbose; fprintf('[simulation] Loading spatial background map from: %s\n', spatialBgMap); end
    tmp = load(spatialBgMap, 'bgMap');
    spatialBgMap = tmp.bgMap;
end
useSpatialBg = isstruct(spatialBgMap) && isfield(spatialBgMap, 'cumProb');
if useSpatialBg && verbose
    fprintf('[simulation] Using learned spatial background distribution\n');
elseif verbose
    fprintf('[simulation] Using uniform spatial background\n');
end

rng(seed, 'twister');

%% ---- Load empirical mark distributions ----------------------------------
if verbose; fprintf('[simulation] Loading empirical marks from: %s\n', dataFilePath); end

[~, t_raw, ~, ~, aa_raw, bb_raw, ~, ~] = ...
    textread(dataFilePath, '%s%d%d%d%f%f%d%d');  %#ok<DTXTRD>

% Apply the same scaling as ExpKernelEstimation.m (lines 69-71)
time_scalar = 1;
aa_emp = aa_raw * (60/100) * 15;       % amplitude scaling
bb_emp = bb_raw * 125;                  % radius scaling
lt_emp = double(t_raw) * time_scalar / 15;  % lifetime in time-step blocks

% Remove any entries with non-positive radius (guard for strain = |a|/b)
valid = bb_emp > 0;
aa_emp = aa_emp(valid);
bb_emp = bb_emp(valid);
lt_emp = lt_emp(valid);

nEmp = numel(aa_emp);
if verbose; fprintf('  Empirical mark pool: %d events\n', nEmp); end

%% ---- Helper: sample marks from empirical pool ---------------------------
    function [a, b, l] = sampleMarks(nSamp)
        idx = randi(nEmp, nSamp, 1);
        a = aa_emp(idx);
        b = bb_emp(idx);
        l = lt_emp(idx);
    end

%% ---- Helper: truncated exponential draw ---------------------------------
    function tau = truncExp(rate, bound, nSamp)
        % Draw from Exp(rate) truncated to [0, bound]
        % CDF inversion: tau = -log(1 - U*(1 - exp(-rate*bound))) / rate
        U = rand(nSamp, 1);
        tau = -log(1 - U .* (1 - exp(-rate * bound))) / rate;
    end

%% ---- Generation 0: Background events -----------------------------------
expectedBg = lambda0 * GridArea * T;
N_bg = poissrnd(expectedBg);

if verbose
    fprintf('[simulation] Background: E[N_bg] = %.1f, sampled N_bg = %d\n', expectedBg, N_bg);
end

% Pre-allocate the catalog as arrays that expand with each generation.
allTime     = zeros(N_bg, 1);
allX        = zeros(N_bg, 1);
allY        = zeros(N_bg, 1);
allAA       = zeros(N_bg, 1);
allBB       = zeros(N_bg, 1);
allLifetime = zeros(N_bg, 1);
allGen      = zeros(N_bg, 1);
allParent   = zeros(N_bg, 1);

if N_bg > 0
    allTime(1:N_bg)     = rand(N_bg, 1) * T;               % Uniform(0, T)
    if useSpatialBg
        % Sample locations from learned spatial density via inverse-CDF
        U = rand(N_bg, 1);
        edges = [0; spatialBgMap.cumProb(:)];
        [~, ~, flatIdx] = histcounts(U, edges);
        flatIdx = max(flatIdx, 1);  % Guard against edge cases
        gs = spatialBgMap.GridSide;
        rowIdx = mod(flatIdx - 1, gs) + 1;
        colIdx = floor((flatIdx - 1) / gs) + 1;
        % Grid center + uniform sub-cell jitter
        allX(1:N_bg) = double(spatialBgMap.gridX(colIdx)) + rand(N_bg, 1) - 0.5;
        allY(1:N_bg) = double(spatialBgMap.gridY(rowIdx)) + rand(N_bg, 1) - 0.5;
        % Clamp to domain boundaries
        allX(1:N_bg) = max(0, min(allX(1:N_bg), GridSide));
        allY(1:N_bg) = max(0, min(allY(1:N_bg), GridSide));
    else
        allX(1:N_bg)    = rand(N_bg, 1) * GridSide;         % Uniform(0, GridSide)
        allY(1:N_bg)    = rand(N_bg, 1) * GridSide;
    end
    [a0, b0, l0]        = sampleMarks(N_bg);
    allAA(1:N_bg)       = a0;
    allBB(1:N_bg)       = b0;
    allLifetime(1:N_bg) = l0;
    allGen(1:N_bg)      = 0;                                % generation 0
    allParent(1:N_bg)   = 0;                                % no parent
end

%% ---- Branching: generate offspring generation-by-generation -------------
genStart = 1;             % index of first event in current generation
genEnd   = N_bg;          % index of last event in current generation
totalEvents = N_bg;
gen = 0;

while genStart <= genEnd && gen < maxGen
    gen = gen + 1;
    nParents = genEnd - genStart + 1;
    
    % For each parent in the current generation, compute expected offspring
    newTime     = [];
    newX        = [];
    newY        = [];
    newAA       = [];
    newBB       = [];
    newLifetime = [];
    newParent   = [];
    
    for pidx = genStart:genEnd
        ai = allAA(pidx);
        bi = allBB(pidx);
        li = allLifetime(pidx);
        ti = allTime(pidx);
        xi = allX(pidx);
        yi = allY(pidx);
        
        if bi <= 0; continue; end
        
        strain_i = abs(ai) / bi;
        
        % Expected offspring: K0 * strain * (1 - exp(-omega * l))
        mi = K0 * strain_i * (1 - exp(-omega * li));
        
        % Draw number of offspring
        Ni = poissrnd(mi);
        if Ni == 0; continue; end
        
        % Draw temporal delays from TruncatedExponential(omega, li)
        tau = truncExp(omega, li, Ni);
        t_offspring = ti + tau;
        
        % Keep only offspring within [0, T]
        keep = t_offspring < T & t_offspring >= 0;
        Ni_kept = sum(keep);
        if Ni_kept == 0; continue; end
        
        t_offspring = t_offspring(keep);
        
        % Draw spatial offsets from N(0, 1/(2*rho) * I_2)
        sigma = sqrt(1 / (2 * rho));
        dx = sigma * randn(Ni_kept, 1);
        dy = sigma * randn(Ni_kept, 1);
        x_offspring = xi + dx;
        y_offspring = yi + dy;
        
        % Sample marks for offspring
        [a_off, b_off, l_off] = sampleMarks(Ni_kept);
        
        newTime     = [newTime;     t_offspring];       %#ok<AGROW>
        newX        = [newX;        x_offspring];       %#ok<AGROW>
        newY        = [newY;        y_offspring];       %#ok<AGROW>
        newAA       = [newAA;       a_off];             %#ok<AGROW>
        newBB       = [newBB;       b_off];             %#ok<AGROW>
        newLifetime = [newLifetime; l_off];             %#ok<AGROW>
        newParent   = [newParent;   repmat(pidx, Ni_kept, 1)]; %#ok<AGROW>
    end
    
    nNew = numel(newTime);
    
    if verbose
        fprintf('  Generation %d: %d parents → %d offspring\n', gen, nParents, nNew);
    end
    
    if nNew == 0; break; end
    
    % Append to catalog
    idx1 = totalEvents + 1;
    idx2 = totalEvents + nNew;
    
    allTime(idx1:idx2)     = newTime;
    allX(idx1:idx2)        = newX;
    allY(idx1:idx2)        = newY;
    allAA(idx1:idx2)       = newAA;
    allBB(idx1:idx2)       = newBB;
    allLifetime(idx1:idx2) = newLifetime;
    allGen(idx1:idx2)      = gen;
    allParent(idx1:idx2)   = newParent;
    
    genStart    = idx1;
    genEnd      = idx2;
    totalEvents = idx2;
end

%% ---- Trim and sort chronologically ------------------------------------
allTime     = allTime(1:totalEvents);
allX        = allX(1:totalEvents);
allY        = allY(1:totalEvents);
allAA       = allAA(1:totalEvents);
allBB       = allBB(1:totalEvents);
allLifetime = allLifetime(1:totalEvents);
allGen      = allGen(1:totalEvents);
allParent   = allParent(1:totalEvents);

[allTime, sortIdx] = sort(allTime);
allX        = allX(sortIdx);
allY        = allY(sortIdx);
allAA       = allAA(sortIdx);
allBB       = allBB(sortIdx);
allLifetime = allLifetime(sortIdx);
allGen      = allGen(sortIdx);
allParent   = allParent(sortIdx);

% Parent indices were assigned in generation order. After chronological
% sorting, remap every non-background parent to its new catalog index.
oldToNew = zeros(totalEvents, 1);
oldToNew(sortIdx) = (1:totalEvents)';
hasParent = allParent > 0;
allParent(hasParent) = oldToNew(allParent(hasParent));

%% ---- Build output struct -----------------------------------------------
catalog = struct();
catalog.time     = allTime;
catalog.x        = allX;
catalog.y        = allY;
catalog.aa       = allAA;
catalog.bb       = allBB;
catalog.lifetime = allLifetime;
catalog.gen      = allGen;
catalog.parent   = allParent;
catalog.N        = totalEvents;

% Summary statistics
nBg  = sum(allGen == 0);
nOff = totalEvents - nBg;
maxGenReached = max(allGen);

catalog.summary.N_total      = totalEvents;
catalog.summary.N_background = nBg;
catalog.summary.N_offspring  = nOff;
catalog.summary.frac_bg      = nBg / max(totalEvents, 1);
catalog.summary.max_gen      = maxGenReached;
catalog.summary.params       = struct('lambda0', lambda0, 'K0', K0, ...
                                       'omega', omega, 'rho', rho);
catalog.summary.T            = T;
catalog.summary.GridSide     = GridSide;
catalog.summary.useSpatialBg = useSpatialBg;

%% ---- Print summary -----------------------------------------------------
if verbose
    fprintf('\n=== SIMULATION SUMMARY ===\n');
    fprintf('  Parameters: lambda0=%.4e, K0=%.4e, omega=%.4f, rho=%.6f\n', ...
            lambda0, K0, omega, rho);
    fprintf('  Window:     T=%d time steps, GridArea=%d\n', T, GridArea);
    fprintf('  Total events:     %d\n', totalEvents);
    fprintf('  Background:       %d (%.1f%%)\n', nBg, 100*nBg/max(totalEvents,1));
    fprintf('  Offspring:        %d (%.1f%%)\n', nOff, 100*nOff/max(totalEvents,1));
    fprintf('  Max generation:   %d\n', maxGenReached);
    
    % Empirical branching ratio
    if nBg > 0
        empBranching = nOff / totalEvents;
        fprintf('  Empirical branching ratio: %.4f\n', empBranching);
    end
    
    % Theoretical branching ratio
    meanStrain = mean(abs(aa_emp) ./ bb_emp);
    meanTimeInt = mean(1 - exp(-omega * lt_emp));
    theoBranching = K0 * meanStrain * meanTimeInt;
    fprintf('  Theoretical branching (K0*E[|a|/b]*E[1-exp(-wl)]): %.4f\n', theoBranching);
    fprintf('==========================\n\n');
end

%% ---- Save ---------------------------------------------------------------
if ~isempty(outputFile)
    if verbose; fprintf('[simulation] Saving catalog to: %s\n', outputFile); end
    save(outputFile, 'catalog', 'lambda0', 'K0', 'omega', 'rho', 'T', ...
         'GridSide', 'seed');
end

end
