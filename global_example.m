%% Final Project: Dynamic Gompertz Model for Salmonella Enteritidis
%  Based on global_example.m template
%  Primary: Gompertz differential form
%  Secondary: Modified Ratkowsky equation
clear; clc; close all;

%% Read data
global tTemp

% Temperature data: first column = time (hr), second column = temp (deg C)
temp_raw = readmatrix('Salmonella sin growth Temps.xlsx');
tTemp = [temp_raw(:,1)/60, temp_raw(:,2)];  % convert min to hr

% Growth data
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

%% Initial guesses
% beta = [A, C, M, a, b]
% Fixed: Tmin = 6, Tmax = 46.3
beta0 = [log10(400), 11, 7.5, 0.000338, 0.275];

%% SSC (Scaled Sensitivity Coefficients)
logN_base = gompertzFOR(beta0, tobs);
Np = length(beta0);
pnames = {'A','C','M','a','b'};
SSC = zeros(length(tobs), Np);
delta = 1e-4;
for j = 1:Np
    bp = beta0;
    dp = max(abs(beta0(j)) * delta, 1e-10);
    bp(j) = beta0(j) + dp;
    logN_pert = gompertzFOR(bp, tobs);
    SSC(:, j) = beta0(j) * (logN_pert - logN_base) / dp;
end

figure;
plot(tobs, SSC, 'o-', 'LineWidth', 1.2, 'MarkerSize', 5);
xlabel('Time (hr)'); ylabel('SSC');
title('Scaled Sensitivity Coefficients');
legend(pnames, 'Location', 'best');
grid on;

%% nlinfit
opts = statset('MaxIter', 500, 'TolFun', 1e-10, 'TolX', 1e-10);
[beta, resids, J, COVB, mse] = nlinfit(tobs, yobs, @gompertzINV, beta0, opts);

% Results
se = sqrt(diag(COVB));
ci = nlparci(beta, resids, 'jacobian', J);
fprintf('\n--- Estimated Parameters ---\n');
for j = 1:Np
    fprintf('  %s = %.6f  (SE=%.6f)  95%%CI [%.6f, %.6f]\n', ...
        pnames{j}, beta(j), se(j), ci(j,1), ci(j,2));
end
fprintf('  Tmin = 6.000000  (fixed)\n');
fprintf('  Tmax = 46.300000  (fixed)\n');

% Goodness of fit
SSR = sum(resids.^2);
SST = sum((yobs - mean(yobs)).^2);
RMSE = sqrt(SSR / length(yobs));
R2 = 1 - SSR / SST;
fprintf('\n--- Goodness of Fit ---\n');
fprintf('  RMSE      = %.4f log10 CFU/mL\n', RMSE);
fprintf('  Pseudo-R2 = %.4f\n', R2);

%% Plot: logN(t) observed, logN(t) predicted, T(t)
t_fine = linspace(0, max(tobs), 500)';
logN_pred = gompertzFOR(beta, t_fine);
T_fine = interp1(tTemp(:,1), tTemp(:,2), t_fine, 'linear', 'extrap');

figure;
yyaxis left
plot(tobs, yobs, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
hold on;
plot(t_fine, logN_pred, 'b-', 'LineWidth', 2);
ylabel('log_{10}(N) (CFU/mL)');

yyaxis right
plot(t_fine, T_fine, 'r--', 'LineWidth', 1.5);
ylabel('Temperature (°C)');

xlabel('Time (hr)');
title('Gompertz Dynamic Model: SE Growth (Sinusoidal Temperature)');
legend('Observed', 'Predicted', 'Temperature', 'Location', 'northwest');
grid on;

%% Residual plot
logN_pred_data = gompertzFOR(beta, tobs);
figure;
plot(tobs, yobs - logN_pred_data, 'ko', 'MarkerSize', 8);
hold on; yline(0, 'r--');
xlabel('Time (hr)');
ylabel('Residual (log_{10} CFU/mL)');
title('Residuals');
grid on;

%% ==================== INVERSE function (for nlinfit) ====================
function y = gompertzINV(beta, t)
    global tTemp
    A    = beta(1);
    C    = beta(2);
    M    = beta(3);
    a    = beta(4);
    b    = beta(5);
    Tmin = 6;       % fixed
    Tmax = 46.3;    % fixed

    % unique times for ode45 (handles duplicate times from replicates)
    [t_unique, ~, ic] = unique(t);

    % initial condition from Gompertz at t=0
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

%% ==================== FORWARD function (for plotting/SSC) ====================
function y = gompertzFOR(beta, t)
    global tTemp
    A    = beta(1);
    C    = beta(2);
    M    = beta(3);
    a    = beta(4);
    b    = beta(5);
    Tmin = 6;
    Tmax = 46.3;

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
