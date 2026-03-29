%% Final Project: Dynamic Gompertz Model for Salmonella Enteritidis
%  Sinusoidal temperature profile
%  Primary: Gompertz differential form
%  Secondary: Modified Ratkowsky equation
%  OLS parameter estimation via nlinfit
clear; clc; close all;

%% ==================== 1. Read Data ====================
global tTemp

% Temperature data: [time_hr, T_degC]
temp_raw = readmatrix('Salmonella sin growth Temps.xlsx');
tTemp = [temp_raw(:,1)/60, temp_raw(:,2)];  % convert min -> hr

% Growth data
growth_raw = readmatrix('Salmonella sin growth.xlsx');
t_all = [];
logN_all = [];
for i = 1:size(growth_raw, 1)
    t_hr = growth_raw(i, 1);
    cfu1 = growth_raw(i, 2);
    cfu2 = growth_raw(i, 3);
    if ~isnan(cfu1) && cfu1 > 0
        t_all = [t_all; t_hr];
        logN_all = [logN_all; log10(cfu1)];
    end
    if ~isnan(cfu2) && cfu2 > 0
        t_all = [t_all; t_hr];
        logN_all = [logN_all; log10(cfu2)];
    end
end
Ndata = length(t_all);
fprintf('Growth data: %d points\n', Ndata);
fprintf('Temp data:   %d points, %.2f to %.2f hr\n', ...
    size(tTemp,1), min(tTemp(:,1)), max(tTemp(:,1)));

%% ==================== 2. Initial Guesses ====================
% Parameters: [A, C, M, a, b, Tmin, Tmax]
A0    = log10(400);   % ~2.60
C0    = 11;
M0    = 7.5;
a0    = 0.000338;
b0    = 0.275;
Tmin0 = 6;
Tmax0 = 46.3;

p_all = [A0, C0, M0, a0, b0, Tmin0, Tmax0];
pnames_all = {'A','C','M','a','b','Tmin','Tmax'};

%% ==================== 3. SSC Analysis ====================
fprintf('\n--- Computing SSCs ---\n');
logN_base = gompertzFOR(p_all, t_all);

Np_all = length(p_all);
SSC = zeros(Ndata, Np_all);
delta = 1e-4;
for j = 1:Np_all
    p_pert = p_all;
    dp = max(abs(p_all(j)) * delta, 1e-10);
    p_pert(j) = p_all(j) + dp;
    logN_pert = gompertzFOR(p_pert, t_all);
    SSC(:, j) = p_all(j) * (logN_pert - logN_base) / dp;
end

figure('Name','Scaled Sensitivity Coefficients','Position',[100 100 800 500]);
plot(t_all, SSC, 'o-', 'LineWidth', 1.2, 'MarkerSize', 5);
xlabel('Time (hr)'); ylabel('SSC');
title('Scaled Sensitivity Coefficients');
legend(pnames_all, 'Location', 'best');
grid on;
saveas(gcf, 'SSC_plot.png');

fprintf('  Max |SSC| for each parameter:\n');
for j = 1:Np_all
    fprintf('    %5s: %.4f\n', pnames_all{j}, max(abs(SSC(:,j))));
end
fprintf('  >> Examine SSC plot to decide which params to fix.\n');

%% ==================== 4. Parameter Estimation ====================
% ---- Adjust after examining SSC plot ----
% Default: fix Tmin and Tmax (typically not estimable from one profile)
est_idx = [1, 2, 3, 4, 5];   % estimate A, C, M, a, b
fix_idx = [6, 7];             % fix Tmin = 6, Tmax = 46.3
% ------------------------------------------

p_fixed = p_all;
beta0 = p_all(est_idx);
pnames_est = pnames_all(est_idx);
Np_est = length(est_idx);

fprintf('\n--- Parameter Estimation (nlinfit) ---\n');
fprintf('  Estimating: '); fprintf('%s ', pnames_est{:}); fprintf('\n');
fprintf('  Fixed:      ');
for j = fix_idx
    fprintf('%s=%.2f ', pnames_all{j}, p_fixed(j));
end
fprintf('\n');

% nlinfit
opts = statset('MaxIter', 500, 'TolFun', 1e-10, 'TolX', 1e-10);
[beta_opt, resids, J, CovB, MSE] = nlinfit(t_all, logN_all, ...
    @(beta,t) gompertzINV(beta, t, est_idx, fix_idx, p_fixed), beta0, opts);

% Reconstruct full parameter vector
p_opt = p_fixed;
p_opt(est_idx) = beta_opt;

% Results
se = sqrt(diag(CovB));
ci = nlparci(beta_opt, resids, 'jacobian', J);
fprintf('\n  Estimated parameters:\n');
for j = 1:Np_est
    fprintf('    %5s = %10.6f  (SE=%.6f)  95%%CI [%.6f, %.6f]\n', ...
        pnames_est{j}, beta_opt(j), se(j), ci(j,1), ci(j,2));
end
for j = fix_idx
    fprintf('    %5s = %10.6f  (fixed)\n', pnames_all{j}, p_opt(j));
end

%% ==================== 5. Goodness of Fit ====================
SSR = sum(resids.^2);
SST = sum((logN_all - mean(logN_all)).^2);
RMSE = sqrt(SSR / Ndata);
R2 = 1 - SSR / SST;

fprintf('\n--- Goodness of Fit ---\n');
fprintf('  RMSE      = %.4f log10 CFU/mL\n', RMSE);
fprintf('  Pseudo-R2 = %.4f\n', R2);
fprintf('  MSE       = %.6f\n', MSE);
fprintf('  N         = %d\n', Ndata);

%% ==================== 6. Final Plot ====================
t_fine = linspace(0, max(t_all), 500)';
logN_pred_fine = gompertzFOR(p_opt, t_fine);
T_fine = interp1(tTemp(:,1), tTemp(:,2), t_fine, 'linear', 'extrap');

figure('Name','Model Prediction','Position',[100 100 900 550]);

yyaxis left
plot(t_all, logN_all, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
hold on;
plot(t_fine, logN_pred_fine, 'b-', 'LineWidth', 2);
ylabel('log_{10}(N) (CFU/mL)');
set(gca, 'YColor', 'b');

yyaxis right
plot(t_fine, T_fine, 'r--', 'LineWidth', 1.5);
ylabel('Temperature (°C)');
set(gca, 'YColor', 'r');

xlabel('Time (hr)');
title('Dynamic Gompertz Model: SE Growth under Sinusoidal Temperature');
legend('Observed log_{10}(N)', 'Predicted log_{10}(N)', 'Temperature', ...
    'Location', 'northwest');
grid on;
saveas(gcf, 'final_plot.png');

%% ==================== 7. Residual Plot ====================
logN_pred_data = gompertzFOR(p_opt, t_all);

figure('Name','Residuals','Position',[100 100 700 400]);
plot(t_all, logN_all - logN_pred_data, 'ko', 'MarkerSize', 8);
hold on; yline(0, 'r--', 'LineWidth', 1);
xlabel('Time (hr)');
ylabel('Residual (log_{10} CFU/mL)');
title('Residuals: Observed - Predicted');
grid on;
saveas(gcf, 'residual_plot.png');

fprintf('\nDone. Figures saved: SSC_plot.png, final_plot.png, residual_plot.png\n');

%% ==================== INVERSE FUNCTION (for nlinfit) ====================
function y = gompertzINV(beta, t, est_idx, fix_idx, p_fixed)
    % Maps estimated params -> full param vector -> ODE solve
    % Handles duplicate times via unique + expansion
    global tTemp

    % Reconstruct full parameter vector
    p = p_fixed;
    p(est_idx) = beta;

    A    = p(1);
    C    = p(2);
    M    = p(3);
    a    = p(4);
    b    = p(5);
    Tmin = p(6);
    Tmax = p(7);

    % Get unique times (ode45 needs monotonically increasing tspan)
    [t_unique, ~, ic] = unique(t);

    % Initial condition from Gompertz analytical form at t=0
    T0 = interp1(tTemp(:,1), tTemp(:,2), 0, 'linear', 'extrap');
    mu0 = mu_ratkowsky(T0, a, b, Tmin, Tmax);
    if mu0 > 0
        y0 = A + C * exp(-exp(mu0 * exp(1) * M / C + 1));
    else
        y0 = A;
    end

    % Solve ODE on unique times
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [~, y_unique] = ode45(@ff, t_unique, y0, opts);

    % Expand back to match all data points (including replicates)
    y = y_unique(ic);

    % Clamp
    y = max(y, A);
    y = min(y, A + C);

    % Nested ODE function — accesses beta, tTemp from parent scope
    function dydt = ff(t, y)
        T = interp1(tTemp(:,1), tTemp(:,2), t, 'linear', 'extrap');
        mu = mu_ratkowsky(T, a, b, Tmin, Tmax);
        ratio = (y - A) / C;
        if ratio <= 1e-15 || ratio >= (1 - 1e-15)
            dydt = 0;
        else
            dydt = -mu * exp(1) * ratio * log(ratio);
        end
    end
end

%% ==================== FORWARD FUNCTION (for plotting) ====================
function y = gompertzFOR(p, t)
    % Forward model: solve ODE on given times, return log10(N)
    global tTemp

    A    = p(1);
    C    = p(2);
    M    = p(3);
    a    = p(4);
    b    = p(5);
    Tmin = p(6);
    Tmax = p(7);

    [t_unique, ~, ic] = unique(t);

    % Initial condition
    T0 = interp1(tTemp(:,1), tTemp(:,2), 0, 'linear', 'extrap');
    mu0 = mu_ratkowsky(T0, a, b, Tmin, Tmax);
    if mu0 > 0
        y0 = A + C * exp(-exp(mu0 * exp(1) * M / C + 1));
    else
        y0 = A;
    end

    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [~, y_unique] = ode45(@ff, t_unique, y0, opts);

    y = y_unique(ic);
    y = max(y, A);
    y = min(y, A + C);

    function dydt = ff(t, y)
        T = interp1(tTemp(:,1), tTemp(:,2), t, 'linear', 'extrap');
        mu = mu_ratkowsky(T, a, b, Tmin, Tmax);
        ratio = (y - A) / C;
        if ratio <= 1e-15 || ratio >= (1 - 1e-15)
            dydt = 0;
        else
            dydt = -mu * exp(1) * ratio * log(ratio);
        end
    end
end

%% ==================== Secondary Model ====================
function mu = mu_ratkowsky(T, a, b, Tmin, Tmax)
    if T <= Tmin || T >= Tmax
        mu = 0;
    else
        mu = a * (T - Tmin)^2 * (1 - exp(b * (T - Tmax)));
    end
end
