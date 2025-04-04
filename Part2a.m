clc;
clear;
close all;

%%%%%%%%%%%%% 1) SYSTEM PARAMETERS %%%%%%%%%%%%%
m = 0.75;         % kg
L = 1.25;         % m
c = 0.15;         % N·m·s
g = 9.81;         % m/s^2
A0 = 4;           % N·m (input)
w0 = 2;           % rad/s (input frequency)

%%%%%%%%%%%%% 2) SIMULATION SETTINGS %%%%%%%%%%%%%
T_sim = 20;
dt = 1e-4;
tspan = 0:dt:T_sim;
x0 = [0.1; 0];

%%%%%%%%%%%%% 3) SIMULATE TRUE SYSTEM %%%%%%%%%%%%%
[t, x] = simulate_system(tspan, x0, m, L, c, A0, w0, g);
q = x(:,1);
q_dot = x(:,2);
u = A0 * sin(w0 * t);  % known input

%%%%%%%%%%%%% 4) SAMPLE DATA %%%%%%%%%%%%%
Ts = 0.1;
sample_indices = 1:round(Ts/dt):length(t);
t_sampled = t(sample_indices);
q_sampled = q(sample_indices);
q_dot_sampled = q_dot(sample_indices);
u_sampled = u(sample_indices);

%%%%%%%%%%%%% 5) FILTER & LEAST SQUARES ESTIMATION %%%%%%%%%%%%%
lambda = 0.1;

phi1 = first_order_filter(-q_dot_sampled, Ts, lambda);
phi2 = first_order_filter(-q_sampled,     Ts, lambda);
phi3 = first_order_filter( u_sampled,     Ts, lambda);

Phi = [phi1(:), phi2(:), phi3(:)];
Y = q_dot_sampled(:);

theta_hat = (Phi' * Phi) \ (Phi' * Y);
a1_hat = theta_hat(1);
a2_hat = theta_hat(2);
b0_hat = theta_hat(3);

L_hat = g / a2_hat;
m_hat = 1 / (b0_hat * L_hat^2);
c_hat = (a1_hat + lambda) * m_hat * L_hat^2;

%%%%%%%%%%%%% 6) DISPLAY RESULTS %%%%%%%%%%%%%
fprintf('--- True Parameters ---\n');
fprintf('True L: %.4f m\n', L);
fprintf('True m: %.4f kg\n', m);
fprintf('True c: %.4f N·m·s\n', c);
fprintf('\n--- Estimated Parameters ---\n');
fprintf('Estimated L: %.4f m\n', L_hat);
fprintf('Estimated m: %.4f kg\n', m_hat);
fprintf('Estimated c: %.4f N·m·s\n', c_hat);

%%%%%%%%%%%%% 7) RE-SIMULATE WITH ESTIMATED PARAMETERS %%%%%%%%%%%%%
[t_hat, x_hat] = simulate_system(tspan, x0, m_hat, L_hat, c_hat, A0, w0, g);
q_hat = x_hat(:,1);
q_dot_hat = x_hat(:,2);

%%%%%%%%%%%%% 8) ERROR ANALYSIS %%%%%%%%%%%%%
error_q = q - q_hat;
error_qdot = q_dot - q_dot_hat;

mse_q = mean(error_q.^2);
mse_qdot = mean(error_qdot.^2);

fprintf('\nMSE value for q(t):     %.6f\n', mse_q);
fprintf('MSE value for q_dot(t):    %.6f\n', mse_qdot);

%%%%%%%%%%%%% 9) PLOTS (ALL IN ONE FIGURE) %%%%%%%%%%%%%
figure('Name','All Results - Θέμα 2(α)','NumberTitle','off');

subplot(4,1,1);
plot(t, q, 'b', t_hat, q_hat, 'r--', 'LineWidth', 1.5);
legend('True q(t)', 'Estimated q̂(t)');
xlabel('Time [s]');
ylabel('q [rad]');
title('Angular Displacement');
grid on;

subplot(4,1,2);
plot(t, q_dot, 'b', t_hat, q_dot_hat, 'r--', 'LineWidth', 1.5);
legend('True q̇(t)', 'Estimated q̇̂(t)');
xlabel('Time [s]');
ylabel('q̇ [rad/s]');
title('Angular Velocity');
grid on;

subplot(4,1,3);
plot(t, error_q, 'k', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Error');
title('$e_q(t) = q(t) - \hat{q}(t)$', 'Interpreter', 'latex');
grid on;

subplot(4,1,4);
plot(t, error_qdot, 'k', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Error');
title('$e_qdot(t) = q_dot(t) - \hat{q_dot}(t)$', 'Interpreter', 'latex');
grid on;

%%%%%%%%%%%%%%% 10) FUNCTIONS %%%%%%%%%%%%%%%

function [t, x] = simulate_system(tspan, x0, m, L, c, A0, w0, g)
    dyn = @(t, x) [x(2); ...
                  -g / L * x(1) ...
                  - c / (m * L^2) * x(2) ...
                  + 1 / (m * L^2) * A0 * sin(w0 * t)];
    [t, x] = ode45(dyn, tspan, x0);
end

function y_filt = first_order_filter(x, Ts, lambda)
    N = length(x);
    y_filt = zeros(size(x));
    alpha = (2 - Ts * lambda) / (2 + Ts * lambda);
    beta  = Ts / (2 + Ts * lambda);
    for k = 2:N
        y_filt(k) = alpha * y_filt(k-1) + beta * (x(k) + x(k-1));
    end
end
