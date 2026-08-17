% run_sweep.m
% =========================================================================
% Parameter sweep over random initial conditions for ExpKernelEstimation.
%
% For every .txt file found in DATA_DIR the script:
%   1. Parses timestamps to compute the mean arrival rate kappa.
%   2. Draws N_RUNS random initial parameter vectors within the bounds:
%         lambda0 in [kappa * LAMBDA0_KAPPA_LO, kappa]  (uniform, kappa-relative)
%         K0      in [K0_LO, K0_HI]                     (log-uniform)
%         omega   in [OMEGA_LO, OMEGA_HI]               (uniform)
%         rho     in [RHO_LO, RHO_HI]                   (uniform)
%   3. Calls ExpKernelEstimation for each draw.
%   4. Saves the full workspace to:
%         RESULTS_DIR/<dataName>/
%           result_<dataName>_l<lambda0>_K<K0>_w<omega>_r<rho>.mat
%      and a summary CSV:
%         RESULTS_DIR/<dataName>/summary.csv
%
% Configuration:
%   Edit the CONFIGURATION block below before running.
%
% Usage (from exp_kernel/ folder):
%   >> run_sweep
%   >> run_sweep   % re-run adds to existing results directory
% =========================================================================

%% ---- CONFIGURATION (edit here) ----------------------------------------
N_RUNS     = 10;       % number of random starts per data file

% Parameter bounds used for the paper fits
%   lambda0: sampled as a fraction of kappa (the per-data mean arrival
%   rate). This keeps the sampling interval valid for each data set.
LAMBDA0_KAPPA_LO = 0.01;   % lower bound = 1 % of kappa
LAMBDA0_KAPPA_HI = 1.00;   % upper bound = 100 % of kappa  (= kappa itself)

K0_LO      = 0.01;    % lower bound for K0  (log-uniform sampling)
K0_HI      = 1e6;     % upper bound for K0

OMEGA_LO   = 0.1;     % lower / upper bounds for omega
OMEGA_HI   = 2.0;

RHO_LO     = 0.1;     % lower / upper bounds for rho
RHO_HI     = 2.0;

GRID_AREA  = 91^2;    % spatial domain area (must match ExpKernelEstimation)
WINDOW_STEPS = 14 * 24 * 4;  % 14 days in 15-minute blocks
RAND_SEED  = 42;      % set to [] to disable seeding (non-reproducible)

% Paths (relative to this file's directory, i.e. exp_kernel/)
SCRIPT_DIR  = fileparts(mfilename('fullpath'));
DATA_DIR    = fullfile(SCRIPT_DIR, '..', 'data');
RESULTS_DIR = fullfile(SCRIPT_DIR, 'results');
%% -----------------------------------------------------------------------

%% ---- Reproducibility --------------------------------------------------
if ~isempty(RAND_SEED)
    rng(RAND_SEED, 'twister');
end

%% ---- Discover data files ----------------------------------------------
listing = dir(fullfile(DATA_DIR, '*.txt'));
if isempty(listing)
    error('run_sweep:noData', 'No .txt files found in: %s', DATA_DIR);
end
fprintf('\n[run_sweep] Found %d data file(s) in %s\n', numel(listing), DATA_DIR);

%% ---- Main loop --------------------------------------------------------
for fi = 1:numel(listing)

    dataFile = fullfile(DATA_DIR, listing(fi).name);
    [~, dataName] = fileparts(listing(fi).name);

    fprintf('\n========================================================\n');
    fprintf('[run_sweep] Data file %d/%d: %s\n', fi, numel(listing), listing(fi).name);
    fprintf('========================================================\n');

    % ---- 1. Compute kappa from data ------------------------------------
    kappa = computeKappa(dataFile, GRID_AREA, WINDOW_STEPS);
    fprintf('  kappa (mean arrival rate) = %.6e\n', kappa);

    % Guard: kappa must be positive and finite
    if ~isfinite(kappa) || kappa <= 0
        warning('run_sweep:invalidKappa', ...
            'kappa=%.4g is invalid for %s — skipping file.', kappa, dataName);
        continue;
    end

    % lambda0 is sampled in [kappa * LO_frac, kappa].
    % This is always a valid range and keeps the init proportional
    % to the actual data rate, regardless of dataset size.
    lambda0_lo = kappa * LAMBDA0_KAPPA_LO;
    lambda0_hi = kappa * LAMBDA0_KAPPA_HI;
    fprintf('  lambda0 sample range: [%.4e, %.4e]\n', lambda0_lo, lambda0_hi);

    % ---- 2. Output directory for this data file -----------------------
    outDir = fullfile(RESULTS_DIR, dataName);
    if ~isfolder(outDir)
        mkdir(outDir);
        fprintf('  Created results dir: %s\n', outDir);
    end

    % ---- 3. Summary table (append-friendly) ---------------------------
    summaryFile = fullfile(outDir, 'summary.csv');
    if ~isfile(summaryFile)
        fid = fopen(summaryFile, 'w');
        fprintf(fid, 'run,init_lambda0,init_K0,init_omega,init_rho,');
        fprintf(fid, 'final_lambda0,final_K0,final_omega,final_rho,');
        fprintf(fid, 'iters,output_file\n');
        fclose(fid);
    end

    % ---- 4. Random sweep -----------------------------------------------
    for r = 1:N_RUNS

        % Sample initial parameters
        l0_init  = unifrnd(lambda0_lo, lambda0_hi);  % kappa-relative, always valid
        K0_init  = 10^(unifrnd(log10(K0_LO), log10(K0_HI)));   % log-uniform
        w_init   = unifrnd(OMEGA_LO, OMEGA_HI);
        rho_init = unifrnd(RHO_LO, RHO_HI);

        fprintf('\n  --- Run %d/%d ---\n', r, N_RUNS);
        fprintf('  Init: lambda0=%.4e  K0=%.4e  omega=%.4f  rho=%.4f\n', ...
                l0_init, K0_init, w_init, rho_init);

        % Build output filename encoding data name + initial params
        matName   = buildResultName(dataName, l0_init, K0_init, w_init, rho_init);
        matPath   = fullfile(outDir, matName);

        % Skip if already computed (rerun safety)
        if isfile(matPath)
            fprintf('  [SKIP] Already exists: %s\n', matName);
            continue;
        end

        % Call estimation
        try
            [fl0, fK0, fw, frho, iters] = ExpKernelEstimation( ...
                dataFile, l0_init, K0_init, w_init, rho_init, matPath, WINDOW_STEPS);

            fprintf('  Final: lambda0=%.4e  K0=%.4e  omega=%.4f  rho=%.4f  (iters=%d)\n', ...
                    fl0, fK0, fw, frho, iters);

            % Append row to summary CSV
            fid = fopen(summaryFile, 'a');
            fprintf(fid, '%d,%.6e,%.6e,%.6f,%.6f,%.6e,%.6e,%.6f,%.6f,%d,%s\n', ...
                    r, l0_init, K0_init, w_init, rho_init, ...
                    fl0, fK0, fw, frho, iters, matName);
            fclose(fid);

        catch ME
            fprintf('  [ERROR] Run %d failed: %s\n', r, ME.message);
            fprintf('          %s\n', ME.identifier);
        end

    end  % runs loop

    fprintf('\n  Summary saved to: %s\n', summaryFile);

end  % files loop

fprintf('\n[run_sweep] All done.\n\n');

%% =========================================================================
%  LOCAL HELPER FUNCTIONS
%% =========================================================================

% ---- computeKappa -------------------------------------------------------
function kappa = computeKappa(dataFile, gridArea, windowSteps)
% COMPUTEKAPPA  Estimate the mean arrival rate kappa = n / (GridArea * T).
%
%   The full observation window includes 15-minute blocks with zero events.

fid = fopen(dataFile, 'r');
if fid < 0
    error('computeKappa:fileNotFound', 'Cannot open: %s', dataFile);
end
raw = textscan(fid, '%s %d %d %d %f %f %d %d');
fclose(fid);

n_total = numel(raw{1});
kappa = n_total / (gridArea * windowSteps);
end

% ---- buildResultName ----------------------------------------------------
function name = buildResultName(dataName, l0, K0, omega, rho)
% BUILDRESULTNAME  Construct a filename encoding data name + init params.
%
%   Format:  result_<dataName>_l<lambda0>_K<K0>_w<omega>_r<rho>.mat
%   Floating point values are formatted with 4 significant figures and
%   dots replaced by 'p' to keep filenames valid on all platforms.

    function s = fmt(v)
        % 4 sig-fig scientific-ish string, dot → p
        s = strrep(sprintf('%.4g', v), '.', 'p');
    end

name = sprintf('result_%s_l%s_K%s_w%s_r%s.mat', ...
               dataName, fmt(l0), fmt(K0), fmt(omega), fmt(rho));
end
