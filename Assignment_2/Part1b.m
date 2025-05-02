clc;
clear;
close all;

% === System Parameters ===
m = 1.315;
b = 0.225;
k = 0.725;

% === Matrix A and B construction ===
a11 = 0;
a12 = 1;
a21 = -k/m;
a22 = -b/m;
A = [a11 a12; a21 a22];
B = [0; 1/m];

% === Input and adaptation parameters ===
c = 2.5;
d = 1;
f = 0;

gamma1 = 0.0475;
gamma2 = 0.0454;
gamma3 = 0.0167;

% === Time Setup ===
tspan = 0:0.01:100;
x0_par = [0; 0; 0; 0; -0.01; -0.01; 0.01];  % Initial conditions (parallel)

% === Pack parameters into struct ===
params = struct('m', m, 'b', b, 'k', k, 'A', A, 'B', B, ...
                'c', c, 'd', d, 'f', f, ...
                'gamma1', gamma1, 'gamma2', gamma2, 'gamma3', gamma3);

% === Solve ODE (Parallel Configuration) ===
[t, x] = ode45(@(t, x) lyap_parallel(t, x, params), tspan, x0_par);

% === Extract Results ===
x1     = x(:, 1); x2     = x(:, 2);
x1_hat = x(:, 3); x2_hat = x(:, 4);
a21_hat = x(:, 5); a22_hat = x(:, 6); b2_hat = x(:, 7);

% === Error ===
e1 = x1 - x1_hat;

% === Parameter Estimates ===
m_hat = 1 ./ b2_hat;
b_hat = -a22_hat .* m_hat;
k_hat = -a21_hat .* m_hat;

% === Plotting Results (Parallel) ===
figure;
plot(t, x1, 'r', t, x1_hat, 'b'); grid on;
title('x_1 and x̂_1 (Parallel Configuration)');
xlabel('time in sec'); ylabel('Displacement');
legend('x_1', 'x̂_1');

figure;
plot(t, e1, 'k'); grid on;
title('Estimation Error e_1 (Parallel)');
xlabel('time in sec'); ylabel('Error');

% Apply mask for plot smoothing
mask = false(size(t));  

for i = 1:length(t)
    if t(i) >= 0.5
        mask(i) = true;
    end
end
t_plot = t(mask);
m_hat_plot = m_hat(mask);
b_hat_plot = b_hat(mask);
k_hat_plot = k_hat(mask);

figure;
subplot(3,1,1);
plot(t_plot, m_hat_plot, 'LineWidth', 1.5); hold on;
yline(m, '--r'); ylabel('m̂(t)'); grid on;
title('Estimated Parameters – Parallel');

subplot(3,1,2);
plot(t_plot, b_hat_plot, 'LineWidth', 1.5); hold on;
yline(b, '--r'); ylabel('b̂(t)'); grid on;

subplot(3,1,3);
plot(t_plot, k_hat_plot, 'LineWidth', 1.5); hold on;
yline(k, '--r'); ylabel('k̂(t)'); xlabel('time in sec'); grid on;

fprintf('\n=== Final Parameter Estimates (Parallel) ===\n');
fprintf('Mass:      true = %.4f,  estimated = %.4f\n', m, m_hat(end));
fprintf('Damping:   true = %.4f,  estimated = %.4f\n', b, b_hat(end));
fprintf('Stiffness: true = %.4f,  estimated = %.4f\n', k, k_hat(end));


% === Series-Parallel Configuration ===

thetam = [0.2 0.2 0.2 0.2];
params_sp = params;
params_sp.gamma1 = 0.062;
params_sp.gamma2 = 0.0597;
params_sp.gamma3 = 1.503;
params_sp.thetam = thetam;

x0_sp = [0; 0; 0; 0; 0.01; 0.01; 0.01];

[t, x] = ode45(@(t, x) lyap_series_parallel(t, x, params_sp), tspan, x0_sp);

x1     = x(:, 1); x2     = x(:, 2);
x1_hat = x(:, 3); x2_hat = x(:, 4);
a21_hat = x(:, 5); a22_hat = x(:, 6); b2_hat = x(:, 7);
e1 = x1 - x1_hat;

m_hat = 1 ./ b2_hat;
b_hat = -a22_hat .* m_hat;
k_hat = -a21_hat .* m_hat;

mask = false(size(t));  

for i = 1:length(t)
    if t(i) >= 0.5
        mask(i) = true;
    end
end

t_plot = t(mask);
m_hat_plot = movmean(m_hat(mask), 150);
b_hat_plot = movmean(b_hat(mask), 150);
k_hat_plot = movmean(k_hat(mask), 150);

figure;
plot(t, x1, 'r', t, x1_hat, 'b'); grid on;
title('x_1 and x̂_1 (Series-Parallel)');
xlabel('time in sec'); ylabel('Displacement');
legend('x_1', 'x̂_1');

figure;
plot(t, e1, 'k'); grid on;
title('Estimation Error e_1 (Series-Parallel)');
xlabel('time in sec'); ylabel('Error');

figure;
subplot(3,1,1);
plot(t_plot, m_hat_plot, 'LineWidth', 1.5); hold on;
yline(m, '--r'); ylabel('m̂(t)'); grid on;
title('Estimated Parameters – Series-Parallel');

subplot(3,1,2);
plot(t_plot, b_hat_plot, 'LineWidth', 1.5); hold on;
yline(b, '--r'); ylabel('b̂(t)'); grid on;

subplot(3,1,3);
plot(t_plot, k_hat_plot, 'LineWidth', 1.5); hold on;
yline(k, '--r'); ylabel('k̂(t)'); xlabel('time in sec'); grid on;

fprintf('\n=== Final Parameter Estimates (Series-Parallel) ===\n');
fprintf('Mass:      true = %.4f,  estimated = %.4f\n', m, m_hat(end));
fprintf('Damping:   true = %.4f,  estimated = %.4f\n', b, b_hat(end));
fprintf('Stiffness: true = %.4f,  estimated = %.4f\n', k, k_hat(end));


% === Estimators ===

function dx = lyap_parallel(t, x, P)
    u = P.c * sin(P.d * t) + P.f;

    e1 = x(1) - x(3);
    e2 = x(2) - x(4);

    dx = zeros(7,1);
    dx(1) = P.A(1,:)*[x(1); x(2)] + P.B(1)*u;
    dx(2) = P.A(2,:)*[x(1); x(2)] + P.B(2)*u;
    dx(3) = x(4);
    dx(4) = x(5)*x(3) + x(6)*x(4) + x(7)*u;
    dx(5) = P.gamma1*e2*x(3);
    dx(6) = P.gamma2*e2*x(4);
    dx(7) = P.gamma3*e2*u;
end

function dx = lyap_series_parallel(t, x, P)
    u = P.c * sin(P.d * t) + P.f;

    e1 = x(1) - x(3);
    e2 = x(2) - x(4);

    dx = zeros(7,1);
    dx(1) = x(2);
    dx(2) = P.A(2,:)*[x(1); x(2)] + P.B(2)*u;
    dx(3) = x(2) + P.thetam(1)*e1 + P.thetam(2)*e2;
    dx(4) = x(5)*x(1) + x(6)*x(2) + x(7)*u + P.thetam(3)*e1 + P.thetam(4)*e2;
    dx(5) = P.gamma1*e2*x(1);
    dx(6) = P.gamma2*e2*x(2);
    dx(7) = P.gamma3*e2*u;
end
