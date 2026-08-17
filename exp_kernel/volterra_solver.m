function [t_grid, Q] = volterra_solver(dataFilePath, lambda0, K0, omega, rho, T, varargin)
%VOLTERRA_SOLVER  Numerical solution of the Volterra integral equation for
%                 the spatially-integrated expected intensity Q(t).
%
%   [t_grid, Q] = volterra_solver(dataFilePath, lambda0, K0, omega, rho, T)
%   [t_grid, Q] = volterra_solver(..., 'Name', Value)
%
%   Solves the Volterra integral equation of the second kind (paper Eq. 20):
%
%     Q(t) = lambda0*|S| + K0*E[|a|/b]*H_bar_S
%             * integral_0^t omega*exp(-omega*(t-s))*gamma_bar(t-s)*Q(s) ds
%
%   where gamma_bar(tau) = P(L > tau) is the lifetime survival function,
%   estimated empirically from the data.
%
%   Inputs:
%     dataFilePath  – path to a Lifetime-*.txt file (for empirical marks)
%     lambda0       – background rate per unit space per unit time
%     K0            – productivity parameter
%     omega         – temporal decay rate
%     rho           – spatial decay rate (retained for model/API consistency;
%                     it drops out of Q when H_bar_S = 1)
%     T             – time window length (in time-step units)
%
%   Name-Value Options:
%     'GridSide'  – side length of the square spatial domain (default 91)
%     'dt'        – quadrature step in 15-minute block units (default 0.25,
%                   corresponding to 3.75 minutes)
%     'Verbose'   – print info (default true)
%
%   Outputs:
%     t_grid – (M x 1) time grid
%     Q      – (M x 1) spatially-integrated expected intensity at each time
%
%   The per-unit-area expected intensity is lambda*(t) = Q(t) / |S|.
%
%   Example (using estimated parameters from first 14 days):
%     [tg, Q] = volterra_solver('../data/Lifetime-first14days.txt', ...
%                                0.00049833, 2.1154, 0.77657, ...
%                                0.020512, 1344);
%     plot(tg, Q / 91^2);  % per-unit-area intensity
%
%   The quadrature follows standard direct methods for second-kind Volterra
%   equations; see Linz (1985) and Brunner (2004).
%
% See also: simulation, ExpKernelEstimation

%% ---- Parse inputs -------------------------------------------------------
p = inputParser;
addRequired(p, 'dataFilePath', ...
    @(v) ischar(v) || (isstring(v) && isscalar(v)));
addRequired(p, 'lambda0',     @(v) isfinite(v) && v > 0);
addRequired(p, 'K0',          @(v) isfinite(v) && v > 0);
addRequired(p, 'omega',       @(v) isfinite(v) && v > 0);
addRequired(p, 'rho',         @(v) isfinite(v) && v > 0);
addRequired(p, 'T',           @(v) isfinite(v) && v > 0);
addParameter(p, 'GridSide', 91,   @(v) isfinite(v) && v > 0);
addParameter(p, 'dt',       0.25, @(v) isfinite(v) && v > 0);
addParameter(p, 'Verbose',  true, @islogical);
parse(p, dataFilePath, lambda0, K0, omega, rho, T, varargin{:});

GridSide = p.Results.GridSide;
GridArea = GridSide^2;
dt       = p.Results.dt;
verbose  = p.Results.Verbose;

%% ---- Load empirical marks for E[|a|/b] and survival function -----------
fid = fopen(dataFilePath, 'r');
if fid < 0
    error('volterra_solver:dataFile', ...
          'Cannot open empirical catalog: %s', dataFilePath);
end
cleanup = onCleanup(@() fclose(fid));
raw = textscan(fid, '%s%f%f%f%f%f%f%f');
t_raw = raw{2};
aa_raw = raw{5};
bb_raw = raw{6};
clear cleanup;

% Same scaling as ExpKernelEstimation
time_scalar = 1;
aa_emp = aa_raw * (60/100) * 15;
bb_emp = bb_raw * 125;
lt_emp = double(t_raw) * time_scalar / 15;

valid = bb_emp > 0;
aa_emp = aa_emp(valid);
bb_emp = bb_emp(valid);
lt_emp = lt_emp(valid);

% E[|a|/b]
mean_strain = mean(abs(aa_emp) ./ bb_emp);

% Empirical survival function gamma_bar(tau) = P(L > tau)
% Build from sorted lifetimes using the ECDF complement
lt_sorted = sort(lt_emp);
nEmp = numel(lt_sorted);

if verbose
    fprintf('[volterra_solver] Empirical marks: %d events\n', nEmp);
    fprintf('  E[|a|/b] = %.6f\n', mean_strain);
    fprintf('  Lifetime range: [%.1f, %.1f] time steps\n', min(lt_emp), max(lt_emp));
end

%% ---- Whole-plane spatial normalization --------------------------------
% The scalar mean-field model adopts H_bar_S = 1. The normalized Gaussian
% spatial kernel then integrates out, so rho does not enter Q(t).
H_bar_S = 1.0;

if verbose
    fprintf('  H_bar_S = %.4f (large-domain approximation)\n', H_bar_S);
end

%% ---- Build time grid ---------------------------------------------------
M = ceil(T / dt) + 1;
t_grid = linspace(0, T, M)';
dt_actual = t_grid(2) - t_grid(1);

%% ---- Precompute survival function on the grid -------------------------
% gamma_bar(tau) = fraction of lifetimes > tau
gamma_bar_grid = zeros(M, 1);
for k = 1:M
    tau = t_grid(k);
    gamma_bar_grid(k) = sum(lt_emp > tau) / nEmp;
end

%% ---- Precompute the convolution kernel K(tau) --------------------------
% K(tau) = K0 * E[|a|/b] * H_bar_S * omega * exp(-omega*tau) * gamma_bar(tau)
K_grid = K0 * mean_strain * H_bar_S * omega * exp(-omega * t_grid) .* gamma_bar_grid;

%% ---- Solve by the implicit composite trapezoidal rule ------------------
Q = zeros(M, 1);
Q(1) = lambda0 * GridArea;  % Q(0) = lambda0 * |S| (no offspring yet)
denominator = 1 - 0.5 * dt_actual * K_grid(1);
if denominator <= 0
    error('volterra_solver:stepTooLarge', ...
          'dt is too large for the trapezoidal update; reduce dt.');
end

for n = 2:M
    history = 0.5 * K_grid(n) * Q(1);
    if n > 2
        history = history + sum(K_grid(n-1:-1:2) .* Q(2:n-1));
    end
    Q(n) = (lambda0 * GridArea + dt_actual * history) / denominator;
end

%% ---- Summary -----------------------------------------------------------
if verbose
    % Branching ratio
    n_branch = K0 * mean_strain * H_bar_S * mean(1 - exp(-omega * lt_emp));
    
    fprintf('\n=== VOLTERRA SOLVER RESULTS ===\n');
    fprintf('  Branching ratio n = %.4f\n', n_branch);
    fprintf('  Q(0) = %.2f,  Q(T) = %.2f\n', Q(1), Q(end));
    fprintf('  Steady state Q_inf = lambda0*|S|/(1-n) = %.2f\n', ...
            lambda0 * GridArea / max(1 - n_branch, 1e-9));
    fprintf('  lambda*(0) = %.6f,  lambda*(T) = %.6f\n', Q(1)/GridArea, Q(end)/GridArea);
    fprintf('==============================\n\n');
end

end
