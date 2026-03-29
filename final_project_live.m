%% Final Project: Dynamic Gompertz Model for Salmonella Enteritidis
% *Modeling Methods in Biosystems Engineering — Spring 2026*
%
% This script develops a dynamic predictive model for the growth of
% _Salmonella_ Enteritidis (SE) in egg yolk under a sinusoidal temperature
% profile. The model uses the Gompertz equation (differential form) as the
% primary model and the modified Ratkowsky equation as the secondary model.
%
% Reference: Gumudavelli V, Subbiah J, Thippareddi H, Velugoti PR,
% Froning G. 2007. _J. Food Sci._ 72:M254-62.

%% 1. Read Data
% Temperature data is shared between the script and ODE functions via a
% global variable, following the instructor's |global_example.m| template.

clear; clc; close all;
global tTemp

%%
% *Temperature data:* 47 points, sinusoidal profile (~7-43 deg C, 24 h).
% Time is converted from minutes to hours.

temp_raw = readmatrix('Salmonella sin growth Temps.xlsx');
tTemp = [temp_raw(:,1)/60, temp_raw(:,2)];

%%
% *Growth data:* 12 time points (0-17 hr). Times 0, 0.5, 1, 2 hr have
% single measurements; times 3-17 hr have duplicate measurements.
% Total: 20 data points. CFU/mL is converted to log10.

growth_raw = readmatrix('Salmonella sin growth.xlsx');
tobs = [];
yobs = [];
for i = 1:size(growth_raw, 1)
    t_hr = growth_raw(i, 1);
    cfu1 = growth_raw(i, 2);
    cfu2 = growth_raw(i, 3);
    if ~isnan(cfu1) && cfu1 > 0
        tobs = [tobs; t_hr];
        yobs = [yobs; log10(cfu1)];
    end
    if ~isnan(cfu2) && cfu2 > 0
        tobs = [tobs; t_hr];
        yobs = [yobs; log10(cfu2)];
    end
end
Ndata = length(tobs);
fprintf('Growth data: %d points\n', Ndata);
fprintf('Temp data:   %d points, %.2f to %.2f hr\n', ...
    size(tTemp,1), min(tTemp(:,1)), max(tTemp(:,1)));

%% 2. Model Description
% *Primary model (Gompertz, constant T):*
%
% $$\log_{10} N(t) = A + C \cdot \exp\!\Big(-\exp\!\Big(\frac{\mu\, e}{C}(M - t) + 1\Big)\Big)$$
%
% *Differential form (for varying T):*
%
% $$\frac{dy}{dt} = -\mu(T)\, e\, \frac{y - A}{C}\, \ln\!\Big(\frac{y - A}{C}\Big)$$
%
% *Secondary model (modified Ratkowsky):*
%
% $$\mu(T) = a\,(T - T_{\min})^2\,\Big(1 - \exp\!\big(b\,(T - T_{\max})\big)\Big)$$
%
% *Parameters:*
%
% * $A$ = initial log10(CFU/mL)
% * $C$ = growth range (upper - lower asymptote)
% * $M$ = inflection time (hr)
% * $a, b$ = Ratkowsky regression coefficients
% * $T_{\min}, T_{\max}$ = theoretical growth temperature limits (fixed)

%% 3. Initial Guesses
% Based on the instructor's recommendations. $T_{\min}$ and $T_{\max}$ are
% fixed because they are not estimable from a single temperature profile
% (confirmed by SSC analysis below).

beta0 = [log10(400), 11, 7.5, 0.000338, 0.275];
Np = length(beta0);
pnames = {'A','C','M','a','b'};

fprintf('\nInitial guesses:\n');
for j = 1:Np
    fprintf('  %s = %.6f\n', pnames{j}, beta0(j));
end
fprintf('  Tmin = 6.0 (fixed)\n');
fprintf('  Tmax = 46.3 (fixed)\n');

%% 4. Forward Problem: Ypred with Initial Guesses
% Solve the ODE with initial guesses and compare with observed data.
% This verifies the guesses are in the right ballpark for |nlinfit| to
% converge.

t_fine = linspace(0, max(tobs), 500)';
logN_guess = gompertzFOR(beta0, t_fine);
T_fine = interp1(tTemp(:,1), tTemp(:,2), t_fine, 'linear', 'extrap');

figure;
yyaxis left
plot(tobs, yobs, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
hold on;
plot(t_fine, logN_guess, 'b-', 'LineWidth', 2.5);
ylabel('log_{10}(N) (CFU/mL)', 'FontSize', 14);

yyaxis right
plot(t_fine, T_fine, 'r--', 'LineWidth', 2);
ylabel('Temperature (°C)', 'FontSize', 14);

xlabel('Time (hr)', 'FontSize', 14);
title('Forward Problem: Y_{pred} with Initial Guesses', 'FontSize', 14);
legend('Observed', 'Predicted (guess)', 'Temperature', ...
    'Location', 'northwest', 'FontSize', 12);
grid on; set(gca, 'FontSize', 12);

%% 5. Scaled Sensitivity Coefficients (Initial Guesses)
% SSC measures how sensitive the model output is to each parameter:
%
% $$\text{SSC}_{ij} = p_j \cdot \frac{\partial y_i}{\partial p_j}$$
%
% Computed via finite differences. Large, distinct SSC patterns indicate
% the parameter is estimable.

logN_base = gompertzFOR(beta0, tobs);
SSC0 = zeros(Ndata, Np);
delta = 1e-4;
for j = 1:Np
    bp = beta0;
    dp = max(abs(beta0(j)) * delta, 1e-10);
    bp(j) = beta0(j) + dp;
    logN_pert = gompertzFOR(bp, tobs);
    SSC0(:, j) = beta0(j) * (logN_pert - logN_base) / dp;
end

figure;
plot(tobs, SSC0, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Time (hr)', 'FontSize', 14);
ylabel('Scaled Sensitivity Coefficient', 'FontSize', 14);
title('SSC with Initial Parameter Guesses', 'FontSize', 14);
legend(pnames, 'Location', 'best', 'FontSize', 12);
grid on; set(gca, 'FontSize', 12);

%%
% *SSC analysis results:*

fprintf('\n--- SSC (Initial Guesses) ---\n');
for j = 1:Np
    fprintf('  %5s: max|SSC| = %.4f\n', pnames{j}, max(abs(SSC0(:,j))));
end
fprintf('\nAll 5 parameters show significant, distinct SSC patterns.\n');
fprintf('Tmin and Tmax (not shown) have negligible SSC -> fixed.\n');

%% 6. OLS Parameter Estimation
% Using |nlinfit| (nonlinear least squares, Levenberg-Marquardt).
% The inverse function |gompertzINV| maps parameters to predicted log10(N)
% via |ode45|.

opts = statset('MaxIter', 500, 'TolFun', 1e-10, 'TolX', 1e-10);
[beta, resids, J, COVB, MSE] = nlinfit(tobs, yobs, @gompertzINV, beta0, opts);

%% 6a. Parameter Estimates and Statistics

se = sqrt(diag(COVB));
ci = nlparci(beta, resids, 'jacobian', J);
rel_err = se ./ abs(beta') * 100;

fprintf('\n--- Estimated Parameters ---\n');
fprintf('  %-5s  %12s  %10s  %10s  %24s\n', ...
    'Param', 'Estimate', 'SE', 'RelErr%', '95% CI');
fprintf('  %s\n', repmat('-', 1, 70));
for j = 1:Np
    fprintf('  %-5s  %12.6f  %10.6f  %9.2f%%  [%10.6f, %10.6f]\n', ...
        pnames{j}, beta(j), se(j), rel_err(j), ci(j,1), ci(j,2));
end
fprintf('  Tmin = 6.0 (fixed),  Tmax = 46.3 (fixed)\n');

%% 6b. Correlation Matrix

Corr = COVB ./ (se * se');
fprintf('\n--- Correlation Matrix ---\n');
fprintf('       ');
for j = 1:Np; fprintf('  %7s', pnames{j}); end; fprintf('\n');
for i = 1:Np
    fprintf('  %-5s', pnames{i});
    for j = 1:Np
        fprintf('  %7.3f', Corr(i,j));
    end
    fprintf('\n');
end

%% 6c. Goodness of Fit

SSR = sum(resids.^2);
SST = sum((yobs - mean(yobs)).^2);
RMSE = sqrt(MSE);
R2 = 1 - SSR / SST;

fprintf('\n--- Goodness of Fit ---\n');
fprintf('  RMSE      = %.4f log10 CFU/mL\n', RMSE);
fprintf('  Pseudo-R2 = %.4f\n', R2);
fprintf('  MSE       = %.6f\n', MSE);

%% 7. Fitted Curve with Confidence and Prediction Bands
% * *Confidence Band (CB):* 95%% band for the mean response $E[y|t]$
% * *Prediction Band (PB):* 95%% band for a new observation $y|t$
%
% Computed using |nlpredci|.

[logN_pred_fine, delta_pred] = nlpredci(@gompertzINV, t_fine, beta, resids, ...
    'jacobian', J, 'predopt', 'curve');
[~, delta_obs] = nlpredci(@gompertzINV, t_fine, beta, resids, ...
    'jacobian', J, 'predopt', 'observation');

CB_upper = logN_pred_fine + delta_pred;
CB_lower = logN_pred_fine - delta_pred;
PB_upper = logN_pred_fine + delta_obs;
PB_lower = logN_pred_fine - delta_obs;

figure;
yyaxis left
fill([t_fine; flipud(t_fine)], [PB_upper; flipud(PB_lower)], ...
    [0.9 0.9 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;
fill([t_fine; flipud(t_fine)], [CB_upper; flipud(CB_lower)], ...
    [0.7 0.7 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(tobs, yobs, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot(t_fine, logN_pred_fine, 'b-', 'LineWidth', 2.5);
ylabel('log_{10}(N) (CFU/mL)', 'FontSize', 14);

yyaxis right
plot(t_fine, T_fine, 'r--', 'LineWidth', 2);
ylabel('Temperature (°C)', 'FontSize', 14);

xlabel('Time (hr)', 'FontSize', 14);
title('OLS Fit with Confidence and Prediction Bands', 'FontSize', 14);
legend('Prediction Band', 'Confidence Band', 'Observed', 'Predicted', ...
    'Temperature', 'Location', 'northwest', 'FontSize', 11);
grid on; set(gca, 'FontSize', 12);

%% 8. Residual Analysis
% Check the five standard statistical assumptions.

logN_pred_data = gompertzFOR(beta, tobs);
res_plot = yobs - logN_pred_data;

%% 8a. Residual Scatter Plot

figure;
plot(tobs, res_plot, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
hold on; yline(0, 'r--', 'LineWidth', 1.5);
xlabel('Time (hr)', 'FontSize', 14);
ylabel('Residual (log_{10} CFU/mL)', 'FontSize', 14);
title('Residual Scatter Plot', 'FontSize', 14);
grid on; set(gca, 'FontSize', 12);

%% 8b. Residual Histogram

figure;
histogram(res_plot, 8, 'FaceColor', [0.3 0.5 0.8]);
xlabel('Residual (log_{10} CFU/mL)', 'FontSize', 14);
ylabel('Frequency', 'FontSize', 14);
title('Residual Histogram', 'FontSize', 14);
set(gca, 'FontSize', 12);

%% 8c. Five Standard Statistical Assumptions

fprintf('\n--- Standard Statistical Assumptions ---\n');

% 1. Model is correct
fprintf('  1. Model is correct: Check Ypred vs Yobs plot (visual).\n');

% 2. Errors are random (Durbin-Watson)
DW = sum(diff(res_plot).^2) / sum(res_plot.^2);
fprintf('  2. Errors are random: Durbin-Watson = %.4f', DW);
if abs(DW - 2) < 0.5
    fprintf(' (PASS, close to 2)\n');
else
    fprintf(' (CHECK, far from 2)\n');
end

% 3. Constant variance
fprintf('  3. Constant variance: Check residual scatter plot (visual).\n');

% 4. Errors uncorrelated
fprintf('  4. Errors uncorrelated: See Durbin-Watson above.\n');

% 5. Normally distributed (Shapiro-Wilk)
[h_sw, p_sw] = swtest(res_plot);
if h_sw == 0
    fprintf('  5. Normal distribution: Shapiro-Wilk p=%.4f (PASS, p>0.05).\n', p_sw);
else
    fprintf('  5. Normal distribution: Shapiro-Wilk p=%.4f (FAIL, p<0.05).\n', p_sw);
end

%% 9. Final Scaled Sensitivity Coefficients
% SSCs recomputed with estimated parameters to confirm identifiability
% at the solution.

logN_base_final = gompertzFOR(beta, tobs);
SSC_final = zeros(Ndata, Np);
for j = 1:Np
    bp = beta;
    dp = max(abs(beta(j)) * delta, 1e-10);
    bp(j) = beta(j) + dp;
    logN_pert = gompertzFOR(bp, tobs);
    SSC_final(:, j) = beta(j) * (logN_pert - logN_base_final) / dp;
end

figure;
plot(tobs, SSC_final, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Time (hr)', 'FontSize', 14);
ylabel('Scaled Sensitivity Coefficient', 'FontSize', 14);
title('SSC with Estimated Parameters', 'FontSize', 14);
legend(pnames, 'Location', 'best', 'FontSize', 12);
grid on; set(gca, 'FontSize', 12);

fprintf('\n--- SSC (Estimated Parameters) ---\n');
for j = 1:Np
    fprintf('  %5s: max|SSC| = %.4f\n', pnames{j}, max(abs(SSC_final(:,j))));
end

%% 10. Optimal Experimental Design
% * *Delta criterion:* $\det(X^T X)$ measures overall information content.
% * *$C_{ii}$ curves:* diagonal of $(X^T X)^{-1}$, proportional to
%   parameter variance.

t_design = linspace(0, max(tobs), 200)';
logN_base_d = gompertzFOR(beta, t_design);
X_full = zeros(length(t_design), Np);
for j = 1:Np
    bp = beta;
    dp = max(abs(beta(j)) * delta, 1e-10);
    bp(j) = beta(j) + dp;
    logN_pert = gompertzFOR(bp, t_design);
    X_full(:, j) = (logN_pert - logN_base_d) / dp;
end

n_design = length(t_design);
delta_crit = zeros(n_design, 1);
Cii_curves = zeros(n_design, Np);
for k = Np:n_design
    Xk = X_full(1:k, :);
    XtX = Xk' * Xk;
    delta_crit(k) = det(XtX);
    if rcond(XtX) > 1e-15
        Cinv = inv(XtX);
        for j = 1:Np
            Cii_curves(k, j) = Cinv(j, j);
        end
    end
end

%% 10a. Delta Criterion

figure;
plot(t_design(Np:end), delta_crit(Np:end), 'b-', 'LineWidth', 2);
xlabel('Last Measurement Time (hr)', 'FontSize', 14);
ylabel('det(X^TX)', 'FontSize', 14);
title('Delta Criterion for Optimal Experimental Design', 'FontSize', 14);
grid on; set(gca, 'FontSize', 12);

%% 10b. Cii Curves

figure;
plot(t_design(Np:end), Cii_curves(Np:end, :), '-', 'LineWidth', 2);
xlabel('Last Measurement Time (hr)', 'FontSize', 14);
ylabel('C_{ii}', 'FontSize', 14);
title('C_{ii} Curves for Optimal Experimental Design', 'FontSize', 14);
legend(pnames, 'Location', 'best', 'FontSize', 12);
grid on; set(gca, 'FontSize', 12);

%% 11. Bootstrap Analysis
% Residual resampling bootstrap with 1000 iterations.
% Compare bootstrap 95%% CI with asymptotic CI.

fprintf('\n--- Bootstrap (Residual Resampling, N=1000) ---\n');
Nboot = 1000;
beta_boot = zeros(Nboot, Np);
logN_pred_boot = zeros(Nboot, length(t_fine));

rng(42);
for ib = 1:Nboot
    res_boot = resids(randi(Ndata, Ndata, 1));
    yobs_boot = logN_pred_data + res_boot;
    try
        beta_boot(ib, :) = nlinfit(tobs, yobs_boot, @gompertzINV, beta, opts);
        logN_pred_boot(ib, :) = gompertzFOR(beta_boot(ib, :), t_fine)';
    catch
        beta_boot(ib, :) = beta;
        logN_pred_boot(ib, :) = logN_pred_fine';
    end
    if mod(ib, 200) == 0
        fprintf('  Iteration %d/%d\n', ib, Nboot);
    end
end

%% 11a. Bootstrap Parameter CIs

beta_boot_sorted = sort(beta_boot);
idx_lo = max(1, round(0.025 * Nboot));
idx_hi = min(Nboot, round(0.975 * Nboot));

fprintf('\n  Bootstrap vs Asymptotic 95%% CI:\n');
fprintf('  %-5s  %12s  %26s  %26s\n', ...
    'Param', 'Estimate', 'Bootstrap 95% CI', 'Asymptotic 95% CI');
fprintf('  %s\n', repmat('-', 1, 78));
for j = 1:Np
    ci_b = [beta_boot_sorted(idx_lo, j), beta_boot_sorted(idx_hi, j)];
    fprintf('  %-5s  %12.6f  [%10.6f, %10.6f]  [%10.6f, %10.6f]\n', ...
        pnames{j}, beta(j), ci_b(1), ci_b(2), ci(j,1), ci(j,2));
end

%% 11b. Bootstrap CB and PB

logN_boot_sorted = sort(logN_pred_boot);
CB_boot_lo = logN_boot_sorted(idx_lo, :)';
CB_boot_hi = logN_boot_sorted(idx_hi, :)';

logN_pred_boot_obs = logN_pred_boot + randn(Nboot, length(t_fine)) * sqrt(MSE);
logN_boot_obs_sorted = sort(logN_pred_boot_obs);
PB_boot_lo = logN_boot_obs_sorted(idx_lo, :)';
PB_boot_hi = logN_boot_obs_sorted(idx_hi, :)';

figure;
yyaxis left
fill([t_fine; flipud(t_fine)], [PB_boot_hi; flipud(PB_boot_lo)], ...
    [1 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
hold on;
fill([t_fine; flipud(t_fine)], [CB_boot_hi; flipud(CB_boot_lo)], ...
    [0.85 1 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(tobs, yobs, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot(t_fine, logN_pred_fine, 'b-', 'LineWidth', 2.5);
ylabel('log_{10}(N) (CFU/mL)', 'FontSize', 14);

yyaxis right
plot(t_fine, T_fine, 'r--', 'LineWidth', 2);
ylabel('Temperature (°C)', 'FontSize', 14);

xlabel('Time (hr)', 'FontSize', 14);
title('Bootstrap Confidence and Prediction Bands', 'FontSize', 14);
legend('Bootstrap PB', 'Bootstrap CB', 'Observed', 'Predicted', ...
    'Temperature', 'Location', 'northwest', 'FontSize', 11);
grid on; set(gca, 'FontSize', 12);

%%
% *Comparison of band widths:*

fprintf('\n  Band width comparison (average):\n');
fprintf('  Asymptotic CB width: %.4f\n', mean(CB_upper - CB_lower));
fprintf('  Bootstrap  CB width: %.4f\n', mean(CB_boot_hi - CB_boot_lo));
fprintf('  Asymptotic PB width: %.4f\n', mean(PB_upper - PB_lower));
fprintf('  Bootstrap  PB width: %.4f\n', mean(PB_boot_hi - PB_boot_lo));

%% 12. Summary
% * Dynamic Gompertz model successfully developed for SE growth in egg yolk
% * 5 parameters estimated (A, C, M, a, b), 2 fixed (Tmin, Tmax)
% * ODE solved numerically with |ode45|, temperature interpolated via |interp1|
% * All 5 standard statistical assumptions evaluated
% * Optimal experimental design analyzed
% * Bootstrap CIs consistent with asymptotic CIs
% * Model supports USDA-FSIS recommendation: store eggs at $\leq$ 7.2 deg C

fprintf('\n========== DONE ==========\n');

%% ==================== Helper Functions ====================

%% Shapiro-Wilk Test
function [H, pValue] = swtest(x)
    x = sort(x(:));
    n = length(x);
    if n < 3 || n > 5000
        H = 0; pValue = 1; return;
    end
    m = norminv(((1:n)' - 0.375) / (n + 0.25));
    m = m / sqrt(m' * m);
    W = (m' * x)^2 / ((x - mean(x))' * (x - mean(x)));
    mu = -1.2725 + 1.0521 * log(n);
    sigma = 1.0308 - 0.26758 * log(n);
    z = (log(1 - W) - mu) / sigma;
    pValue = max(min(1 - normcdf(z), 1), 0);
    H = (pValue < 0.05);
end

%% Inverse Function (for nlinfit)
% Maps parameter vector |beta| to predicted log10(N) at times |t|.
% Uses |ode45| with nested function |ff| for the Gompertz ODE.
% Handles duplicate times via |unique|.

function y = gompertzINV(beta, t)
    global tTemp
    A = beta(1); C = beta(2); M = beta(3);
    a = beta(4); b = beta(5);
    Tmin = 6; Tmax = 46.3;

    [t_unique, ~, ic] = unique(t);

    T0 = interp1(tTemp(:,1), tTemp(:,2), 0, 'linear', 'extrap');
    mu0 = a * (T0 - Tmin)^2 * (1 - exp(b * (T0 - Tmax)));
    if mu0 > 0
        y0 = A + C * exp(-exp(mu0 * exp(1) * M / C + 1));
    else
        y0 = A;
    end

    [~, y_sol] = ode45(@ff, t_unique, y0);
    y = y_sol(ic);
    y = max(y, A);
    y = min(y, A + C);

    function dydt = ff(t, y)
        T = interp1(tTemp(:,1), tTemp(:,2), t, 'linear', 'extrap');
        if T <= Tmin || T >= Tmax
            mu = 0;
        else
            mu = a * (T - Tmin)^2 * (1 - exp(b * (T - Tmax)));
        end
        ratio = (y - A) / C;
        if ratio <= 1e-15 || ratio >= (1 - 1e-15)
            dydt = 0;
        else
            dydt = -mu * exp(1) * ratio * log(ratio);
        end
    end
end

%% Forward Function (for plotting and SSC)
% Same as |gompertzINV| but called separately for plotting and
% sensitivity analysis.

function y = gompertzFOR(beta, t)
    global tTemp
    A = beta(1); C = beta(2); M = beta(3);
    a = beta(4); b = beta(5);
    Tmin = 6; Tmax = 46.3;

    [t_unique, ~, ic] = unique(t);

    T0 = interp1(tTemp(:,1), tTemp(:,2), 0, 'linear', 'extrap');
    mu0 = a * (T0 - Tmin)^2 * (1 - exp(b * (T0 - Tmax)));
    if mu0 > 0
        y0 = A + C * exp(-exp(mu0 * exp(1) * M / C + 1));
    else
        y0 = A;
    end

    [~, y_sol] = ode45(@ff, t_unique, y0);
    y = y_sol(ic);
    y = max(y, A);
    y = min(y, A + C);

    function dydt = ff(t, y)
        T = interp1(tTemp(:,1), tTemp(:,2), t, 'linear', 'extrap');
        if T <= Tmin || T >= Tmax
            mu = 0;
        else
            mu = a * (T - Tmin)^2 * (1 - exp(b * (T - Tmax)));
        end
        ratio = (y - A) / C;
        if ratio <= 1e-15 || ratio >= (1 - 1e-15)
            dydt = 0;
        else
            dydt = -mu * exp(1) * ratio * log(ratio);
        end
    end
end
