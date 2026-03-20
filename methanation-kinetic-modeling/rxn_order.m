% Evaluated at T=250 C

x = [0.01642, 0.01448, 0.00875, 0.00676];  % Cao
y_exp = [0.827040195, 0.820441989, 0.801142857, 0.804733728]; % Ca/Cao

%% Define explicit models
model_1   = @(p, x) (1 + sqrt(1 - 4*p.*x.^-2)) ./ 2;
model0    = @(p, x) 1 - p.*x.^-1;
model05   = @(p, x) 1 - ((p.*x.^(-1/2)).^2)./2 .* (sqrt(1 + 4./(p.*x.^(-1/2)).^2)-1);
model1    = @(p, x) (1 + p).^(-1) .* ones(size(x));
model2    = @(p, x) (sqrt(1 + 4*p.*x) - 1) ./ (2*p.*x);

%% Fit explicit models using lsqcurvefit
p_init = 0.1;  % initial guess
options = optimoptions('lsqcurvefit','Display','off');

p_fit_1   = lsqcurvefit(@(p,x) model_1(p,x), p_init, x, y_exp, [], [], options);
p_fit_0   = lsqcurvefit(@(p,x) model0(p,x), p_init, x, y_exp, [], [], options);
p_fit_05  = lsqcurvefit(@(p,x) model05(p,x), p_init, x, y_exp, [], [], options);
p_fit_1st = lsqcurvefit(@(p,x) model1(p,x), p_init, x, y_exp, [], [], options);
p_fit_2   = lsqcurvefit(@(p,x) model2(p,x), p_init, x, y_exp, [], [], options);

%% 1.5-order (implicit in y) using fsolve + fminsearch
model1_5_fun = @(p,x) arrayfun(@(xi) fsolve(@(y) p*xi^0.5*y^(3/2) + y - 1, 0.8), x);
sse1_5 = @(p) sum((y_exp - model1_5_fun(p,x)).^2);
p_fit_1_5 = fminsearch(sse1_5, p_init);
yfit_1_5 = model1_5_fun(p_fit_1_5, x);

%% 3rd-order (implicit in y) using fsolve + fminsearch
model3_fun = @(p,x) arrayfun(@(xi) fsolve(@(y) (p*xi^2)*y^3 + y - 1, 0.8), x);
sse3 = @(p) sum((y_exp - model3_fun(p,x)).^2);
p_fit_3 = fminsearch(sse3, p_init);
yfit_3 = model3_fun(p_fit_3, x);

%% Compute fitted y for explicit models given calculated p
yfit_1   = model_1(p_fit_1, x);
yfit_0   = model0(p_fit_0, x);
yfit_05  = model05(p_fit_05, x);
yfit_1st = model1(p_fit_1st, x);
yfit_2   = model2(p_fit_2, x);

%% Function to calculate SSE
calcSSE = @(ydata, yfit) sum((ydata - yfit).^2);

%% Compute SSE for all models
SSE_1   = calcSSE(y_exp, yfit_1);
SSE_0   = calcSSE(y_exp, yfit_0);
SSE_05  = calcSSE(y_exp, yfit_05);
SSE_1st = calcSSE(y_exp, yfit_1st);
SSE_1_5 = calcSSE(y_exp, yfit_1_5);
SSE_2   = calcSSE(y_exp, yfit_2);
SSE_3   = calcSSE(y_exp, yfit_3);

%% Display results
fprintf('Fitted parameters and SSE:\n');
fprintf('-1 order: p = %.4f, SSE = %.6f\n', p_fit_1, SSE_1);
fprintf('0 order: p = %.4f, SSE = %.6f\n', p_fit_0, SSE_0);
fprintf('0.5 order: p = %.4f, SSE = %.6f\n', p_fit_05, SSE_05);
fprintf('1st order: p = %.4f, SSE = %.6f\n', p_fit_1st, SSE_1st);
fprintf('1.5 order: p = %.4f, SSE = %.6f\n', p_fit_1_5, SSE_1_5);
fprintf('2nd order: p = %.4f, SSE = %.6f\n', p_fit_2, SSE_2);
fprintf('3rd order: p = %.4f, SSE = %.6f\n', p_fit_3, SSE_3);

%% Function to calculate RMSE
calcRMSE = @(ydata, yfit) sqrt(mean((ydata - yfit).^2));

%% Compute RMSE for explicit models
RMSE_1   = calcRMSE(y_exp, yfit_1);
RMSE_0   = calcRMSE(y_exp, yfit_0);
RMSE_05  = calcRMSE(y_exp, yfit_05);
RMSE_1st = calcRMSE(y_exp, yfit_1st);
RMSE_2   = calcRMSE(y_exp, yfit_2);

%% Compute RMSE for implicit models (1.5- and 3rd-order)
RMSE_1_5 = calcRMSE(y_exp, yfit_1_5);
RMSE_3   = calcRMSE(y_exp, yfit_3);

%% Display RMSE
fprintf('RMSE for each model:\n');
fprintf('-1 order: RMSE = %.6f\n', RMSE_1);
fprintf('0 order: RMSE = %.6f\n', RMSE_0);
fprintf('0.5 order: RMSE = %.6f\n', RMSE_05);
fprintf('1st order: RMSE = %.6f\n', RMSE_1st);
fprintf('1.5 order: RMSE = %.6f\n', RMSE_1_5);
fprintf('2nd order: RMSE = %.6f\n', RMSE_2);
fprintf('3rd order: RMSE = %.6f\n', RMSE_3);

%% Plot experimental data and fitted models
figure; hold on; grid on;
scatter(x, y_exp, 80, 'k', 'filled'); % experimental points

% Plot fitted models
plot(x, yfit_1, '-o', 'DisplayName','-1 order');
plot(x, yfit_0, '-s', 'DisplayName','0 order');
plot(x, yfit_05, '-d', 'DisplayName','0.5 order');
plot(x, yfit_1st, '-^', 'DisplayName','1st order');
plot(x, yfit_1_5, '-v', 'DisplayName','1.5 order');
plot(x, yfit_2, '-p', 'DisplayName','2nd order');
plot(x, yfit_3, '-h', 'DisplayName','3rd order');
set(gcf,'color','white')
xlabel('C_{A0}'); 
ylabel('Ca/Cao, y'); 
title('Reaction Order Model Fits at 250C'); 
legend('Location','best');

% Calculate rate constant
Q = 400/60; %convert sccm to mL/s
tau = 305/Q; %reactor volume 305 mL
k1   = p_fit_1   / tau;
k0   = p_fit_0   / tau;
k05  = p_fit_05  / tau;
k1st = p_fit_1st / tau;
k15  = p_fit_1_5 / tau;
k2   = p_fit_2   / tau;
k3   = p_fit_3   / tau;

fprintf('Fitted parameters and rate constants:\n');
fprintf('-1 order: p = %.6f, k = %.6f\n', p_fit_1, k1);
fprintf('0 order: p = %.6f, k = %.6f\n', p_fit_0, k0);
fprintf('0.5 order: p = %.6f, k = %.6f\n', p_fit_05, k05);
fprintf('1st order: p = %.6f, k = %.6f\n', p_fit_1st, k1st);
fprintf('1.5 order: p = %.6f, k = %.6f\n', p_fit_1_5, k15);
fprintf('2nd order: p = %.6f, k = %.6f\n', p_fit_2, k2);
fprintf('3rd order: p = %.6f, k = %.6f\n', p_fit_3, k3);

