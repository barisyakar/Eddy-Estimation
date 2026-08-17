function [final_lambda0,final_K0,final_omega,final_rho,iter]=ExpKernelEstimation(dataFilePath,init_lambda0,init_K0,init_omega,init_rho,output_workspace_file,window_steps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EM Algorithm with EXPONENTIAL triggering kernel:
%   g(i,j) = K0 * strain_j * omega * exp(-omega*dt) * (rho/pi) * exp(-rho*r^2)
%
% Inputs:
%   dataFilePath        – path to the Lifetime-*.txt data file
%   init_lambda0        – initial background rate
%   init_K0             – initial productivity parameter
%   init_omega          – initial temporal decay rate
%   init_rho            – initial spatial decay rate
%   output_workspace_file – (optional) path to save full workspace .mat
%
% Analytical M-step updates:
%   lambda0 = sum(Pr0) / (GridArea * T)
%   omega   = sum_lhat / sum(Pr .* dt)         (= 1 / mean weighted dt)
%   rho     = sum_lhat / sum(Pr .* r^2)         (= 1 / mean weighted r^2)
%   K0      = sum_lhat / G_sum
%           where G_sum = sum_i strain_i * (1 - exp(-omega*l_i))
%           (spatial integral = 1 for normalized Gaussian kernel on R^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6; output_workspace_file = ''; end
if nargin < 7; window_steps = []; end

% ---- Validate initial parameters ----------------------------------------
% Reject nonfinite or non-positive values before starting the EM iteration.
params_check = [init_lambda0, init_K0, init_omega, init_rho];
params_names = {'init_lambda0', 'init_K0', 'init_omega', 'init_rho'};
for pi = 1:4
    v = params_check(pi);
    if ~isfinite(v) || v <= 0
        error('ExpKernelEstimation:invalidInit', ...
              '%s must be positive and finite, got %.6g.\n%s', ...
              params_names{pi}, v, ...
              'Hint: check that kappa > LAMBDA0_KAPPA_LO in run_sweep.m.');
    end
end

if ~isfile(dataFilePath)
    error('ExpKernelEstimation:fileNotFound', 'Data file not found: %s', dataFilePath);
end
fprintf('Loading data from: %s\n', dataFilePath);
[s,t,firstx,firsty,aa,bb,~,~] = ...
    textread(dataFilePath,'%s%d%d%d%f%f%d%d');  %#ok<DTXTRD>

[eventTimeBlocks, numericTimes] = timestamp_blocks(s);
[~, sortIdx] = sort(numericTimes);
s = s(sortIdx); t = t(sortIdx);
firstx = firstx(sortIdx); firsty = firsty(sortIdx);
aa = aa(sortIdx); bb = bb(sortIdx);
eventTimeBlocks = eventTimeBlocks(sortIdx);
fprintf('Data sorted by timestamp (%d entries)\n', length(s));

time_scalar = 1;
N = length(s);
numberofvortices = zeros(length(unique(eventTimeBlocks)), 1);
aa = aa * (60/100) * 15;
bb = bb * 125;
t = t * time_scalar / 15;

k = 1; m = 1;
vortexarray(k,m) = getvortex(1, s, t, firstx, firsty, aa, bb, 1, eventTimeBlocks(1));
for i = 2:N
    if eventTimeBlocks(i) == eventTimeBlocks(i-1)
        m = m + 1;
        vortexarray(k,m) = getvortex(i, s, t, firstx, firsty, aa, bb, k, eventTimeBlocks(i));
    else
        numberofvortices(k) = m;
        k = k + 1;
        vortexarray(k,1) = getvortex(i, s, t, firstx, firsty, aa, bb, k, eventTimeBlocks(i));
        m = 1;
    end
end
numberofvortices(k) = m;
K = k;
GridArea = 91^2;
if isempty(window_steps)
    if contains(lower(dataFilePath), '14days')
        window_steps = 14 * 24 * 4;
    else
        window_steps = max(eventTimeBlocks) + 1;
    end
end
T = double(window_steps);
if max(eventTimeBlocks) >= T
    error('ExpKernelEstimation:eventOutsideWindow', ...
          'Latest event block %.0f lies outside T=%.0f.', max(eventTimeBlocks), T);
end

% Predecessor structure
for k = 2:K
    n_k = numberofvortices(k);
    j = 0;
    for kk = 1:k-1
        nn = numberofvortices(kk);
        index_num = 0;
        if kk > 1; index_num = sum(numberofvortices(1:kk-1)); end
        for r = 1:nn
            elapsed = vortexarray(k,1).time - vortexarray(kk,r).time;
            if vortexarray(kk,r).t >= elapsed
                j = j + 1;
                for m = 1:n_k
                    predvor(k,m,j) = vortexarray(kk,r);
                    predvor_index(k,m,j) = r + index_num;
                end
            end
        end
    end
end

n = sum(numberofvortices);
Pr0 = ones(n, 1);
Pr = zeros(n, n);
g = zeros(n, n);
lam = zeros(n, 1);
lamzero = zeros(n, 1);
ll = zeros(n, 1);

% Precompute dt_matrix and dist_sq_matrix
dt_matrix = zeros(n, n);
dist_sq_matrix = zeros(n, n);
i = numberofvortices(1);
for k = 2:K
    for m = 1:numberofvortices(k)
        i = i + 1;
        num_pred = length(predvor(k,m,:));
        for j = 1:num_pred
            idx = predvor_index(k,m,j);
            if idx > 0
                ti = vortexarray(k,m).time;
                tj = predvor(k,m,j).time;
                dt_matrix(i, idx) = abs(ti - tj);
                xi = vortexarray(k,m).firstx;
                yi = vortexarray(k,m).firsty;
                xj = predvor(k,m,j).firstx;
                yj = predvor(k,m,j).firsty;
                dist_sq_matrix(i, idx) = (xi - xj)^2 + (yi - yj)^2;
            end
        end
    end
end
fprintf('Precomputed dt_matrix and dist_sq_matrix (%d x %d)\n', n, n);

%% Initialize
newlambda0 = init_lambda0;
newK0 = init_K0;
newomega = init_omega;
newrho = init_rho;

initial_params = [newlambda0, newK0, newomega, newrho];

% Initial E-step
i = numberofvortices(1);
for k = 2:K
    for m = 1:numberofvortices(k)
        i = i + 1;
        ll(i) = length(predvor(k,m,:));
        for j = 1:ll(i)
            src = predvor(k,m,j);
            index_g = predvor_index(k,m,j);
            if index_g > 0
                g(i, index_g) = trigger_exp(vortexarray(k,m), src, initial_params, time_scalar);
            end
        end
        lamzero(i) = newlambda0;
        lam(i) = lamzero(i) + sum(g(i,:));
        Pr(i,:) = g(i,:) / lam(i);
        Pr0(i) = 1 - sum(Pr(i,:));
    end
end

g_nonzero = g(g > 0);
fprintf('\n=== INITIAL E-STEP (Exponential Kernel) ===\n');
fprintf('  Pr0 mean=%.4f, min=%.4f\n', mean(Pr0), min(Pr0));
if ~isempty(g_nonzero)
    fprintf('  g: mean=%.4e, max=%.4e, min=%.4e\n', mean(g_nonzero), max(g_nonzero), min(g_nonzero));
end
fprintf('============================================\n\n');

%% EM Loop
max_iter = 1000;
params_tolerance = 100;
iter = 0;
tolerance_threshold = 5e-4;

fprintf('Initial: lambda0=%.4e, K0=%.4e, omega=%.4f, rho=%.6f\n\n', ...
    newlambda0, newK0, newomega, newrho);

while (params_tolerance > tolerance_threshold) && (iter < max_iter)
    iter = iter + 1;
    
    lambda0 = newlambda0;
    K0 = newK0;
    omega = newomega;
    rho = newrho;
    
    parameters = [lambda0, K0, omega, rho];

    %% M-STEP

    % lambda0: MLE
    newlambda0 = sum(Pr0) / (GridArea * T);

    sum_lhat = sum(sum(Pr));

    % omega: 1 / (weighted mean dt)
    weighted_dt = sum(sum(Pr .* dt_matrix));
    if weighted_dt > 0
        newomega = sum_lhat / weighted_dt;
    else
        newomega = omega;
        fprintf('  WARNING: omega update failed\n');
    end
    newomega = max(0.01, min(newomega, 10.0));  % wider range for exponential

    % rho: 1 / (weighted mean r^2)
    weighted_r2 = sum(sum(Pr .* dist_sq_matrix));
    if weighted_r2 > 0
        newrho = sum_lhat / weighted_r2;
    else
        newrho = rho;
        fprintf('  WARNING: rho update failed\n');
    end
    newrho = max(1e-6, min(newrho, 1.0));

    % K0: sum_lhat / G_sum
    % Time integral: int_0^l omega*exp(-omega*t) dt = 1 - exp(-omega*l)
    % Spatial integral: int (rho/pi)*exp(-rho*r^2) dA = 1 (over R^2)
    G_sum = 0;
    for k = 1:K
        for m = 1:numberofvortices(k)
            v = vortexarray(k,m);
            li = v.t;
            if v.bb > 0
                strain_i = abs(v.aa) / v.bb;
                time_integral = 1 - exp(-newomega * li);
                G_sum = G_sum + strain_i * time_integral;  % spatial = 1
            end
        end
    end
    
    if G_sum > 0
        newK0 = sum_lhat / G_sum;
    else
        newK0 = K0;
    end
    newK0 = max(newK0, 1e-9);

    fprintf('Iter %d: lambda0=%.4e, K0=%.4e, omega=%.4f, rho=%.6f\n', ...
        iter, newlambda0, newK0, newomega, newrho);
    fprintf('  sum_lhat=%.1f, G_sum=%.1f, branching=%.4f\n', ...
        sum_lhat, G_sum, newK0*G_sum/n);

    % Convergence check
    newparams = [newlambda0, newK0, newomega, newrho];
    oldparams = [lambda0, K0, omega, rho];
    rel_change = abs(newparams - oldparams) ./ max(abs(oldparams), 1e-12);
    params_tolerance = max(rel_change);
    
    if all(rel_change < tolerance_threshold)
        params_tolerance = 0;
    end

    %% E-STEP
    newparameters = [newlambda0, newK0, newomega, newrho];
    
    i = numberofvortices(1);
    for k = 2:K
        for m = 1:numberofvortices(k)
            i = i + 1;
            ll(i) = length(predvor(k,m,:));
            for j = 1:ll(i)
                src = predvor(k,m,j);
                index_g = predvor_index(k,m,j);
                if index_g > 0
                    g(i, index_g) = trigger_exp(vortexarray(k,m), src, newparameters, time_scalar);
                end
            end
            lamzero(i) = newlambda0;
            lam(i) = lamzero(i) + sum(g(i,:));
            Pr(i,:) = g(i,:) / lam(i);
            Pr0(i) = 1 - sum(Pr(i,:));
        end
    end

    g_nonzero = g(g > 0);
    fprintf('  Pr0: mean=%.4f, min=%.4f, g_max=%.4e\n\n', ...
        mean(Pr0), min(Pr0), max(g_nonzero));
end

%% Output
final_lambda0 = newlambda0;
final_K0 = newK0;
final_omega = newomega;
final_rho = newrho;

fprintf('=== ExpKernel-EM CONVERGED (iter=%d) ===\n', iter);
fprintf('lambda0=%.6e, K0=%.6e, omega=%.6f, rho=%.6f\n', ...
    final_lambda0, final_K0, final_omega, final_rho);
fprintf('Branching ratio: %.4f (%.1f%% triggered)\n', ...
    final_K0 * G_sum / n, 100 * (1 - mean(Pr0)));

% Characteristic scales
fprintf('\nCharacteristic scales:\n');
fprintf('  Temporal: 1/omega = %.2f time steps (%.0f min)\n', 1/final_omega, 15/final_omega);
fprintf('  Spatial:  1/sqrt(rho) = %.2f grid units\n', 1/sqrt(final_rho));

if ~isempty(output_workspace_file)
    fprintf('Saving workspace to %s\n', output_workspace_file);
    save(output_workspace_file);
end

end
