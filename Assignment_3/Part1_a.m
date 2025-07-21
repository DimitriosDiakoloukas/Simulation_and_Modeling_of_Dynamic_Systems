clc;
clear;
close all;

% Simulation parameters
T   = 50.0;           
dt  = 0.05;          
N   = round(T/dt);
time = linspace(0, T, N).';

A_true = [-2.15,  0.25;
          -0.75, -2.00];
B_true = [0.0; 1.5];

alpha   = 3.0;   
gamma_A = 80.0;     
gamma_B = 50.0;     

% Initial conditions
x = [0.01; 0.01];  
x_hat= [0.0; 0.0];
A_hat= zeros(2,2);
B_hat= zeros(2,1);

% Allocate storage for plotting
x_hist = zeros(N, 2);
x_hat_hist = zeros(N, 2);
e_hist = zeros(N, 2);
A_hat_hist = zeros(N, 2, 2);
B_hat_hist = zeros(N, 2);

% Simulation loop
for k = 1 : N
    t = (k-1)*dt;

    u =  1.0*sin(2*t) ...
       + 0.5*sin(3*t) ...
       + 1*cos(5*t) ...
       + 1*sin(7*t) ...
       + 1*sin(9*t);

    x_dot = A_true * x + B_true * u;

    e = x - x_hat;  
    x_hat_dot = A_hat * x + B_hat * u + alpha * e;

    denomA = 1 + (x.' * x);    % = 1 + ||x||^2
    denomB = 1 + (u^2);

    A_hat_dot = gamma_A * (e * x.') / denomA;   % 2×2 update
    B_hat_dot = gamma_B * (e * u)/ denomB;   % 2×1 update

    %—Euler integration—
    x= x + x_dot* dt;
    x_hat= x_hat + x_hat_dot * dt;
    A_hat= A_hat+ A_hat_dot * dt;
    B_hat = B_hat + B_hat_dot * dt;
    
    %- Apply known constraints-
    A_hat(1,1) = max(min(A_hat(1,1), -1), -3);
    B_hat(2) = max(B_hat(2), 1);
    
    %—Store histories—
    x_hist(k,:) = x.';
    x_hat_hist(k, :) = x_hat.';
    e_hist(k, :)= e.';
    A_hat_hist(k, :, :) = A_hat;
    B_hat_hist(k, :) = B_hat.';
end

A_final = A_hat;
B_final = B_hat;

fprintf('\n=== Final parameter estimates after t = %.2f s ===\n', T);
fprintf('A_hat = [ %.3f   %.3f ;\n          %.3f   %.3f ]\n', ...
        A_final(1,1), A_final(1,2), A_final(2,1), A_final(2,2));
fprintf('A_true = [ %.3f   %.3f ;\n           %.3f   %.3f ]\n', ...
        A_true(1,1),  A_true(1,2),  A_true(2,1),  A_true(2,2));
fprintf('||A_hat − A_true||_F = %.6f\n', norm(A_final - A_true, 'fro'));

fprintf('\nB_hat = [ %.3f ; %.3f ]\n', B_final(1), B_final(2));
fprintf('B_true = [ %.3f ; %.3f ]\n', B_true(1),  B_true(2));
fprintf('||B_hat − B_true||_2 = %.6f\n', norm(B_final - B_true, 2));
fprintf('===============================================\n\n');

% Plotting
figure('Units','normalized','Position',[0.1 0.1 0.8 0.7]);

% True states vs. estimated states
subplot(3,1,1);
plot(time, x_hist(:,1),'b-','LineWidth',1.2); hold on;
plot(time, x_hat_hist(:,1),'b--','LineWidth',1.2);
plot(time, x_hist(:,2),'r-','LineWidth',1.2);
plot(time, x_hat_hist(:,2),'r--','LineWidth',1.2);
title('True States vs. Estimated States','Interpreter','none');
xlabel('Time (s)','Interpreter','none');
ylabel('State Values','Interpreter','none');
legend({ '$x_{1}$ (true)', '$\hat{x}_{1}$ (est)', ...
         '$x_{2}$ (true)', '$\hat{x}_{2}$ (est)' }, ...
       'Interpreter','latex','Location','best');
grid on;

% State estimation error
subplot(3,1,2);
plot(time, e_hist(:,1),'m-','LineWidth',1.2); hold on;
plot(time, e_hist(:,2),'c-','LineWidth',1.2);
title('State Estimation Error','Interpreter','none');
xlabel('Time (s)','Interpreter','none');
ylabel('Error','Interpreter','none');
legend({ '$e_{1} = x_{1} - \hat{x}_{1}$', '$e_{2} = x_{2} - \hat{x}_{2}$' }, ...
       'Interpreter','latex','Location','best');
grid on;

% Parameter estimates vs. true
subplot(3,1,3);
A11_hat = squeeze(A_hat_hist(:,1,1));
A12_hat = squeeze(A_hat_hist(:,1,2));
A21_hat = squeeze(A_hat_hist(:,2,1));
A22_hat = squeeze(A_hat_hist(:,2,2));
B1_hat  = B_hat_hist(:, 1);
B2_hat  = B_hat_hist(:, 2);

plot(time, A11_hat,'b-','LineWidth',1.0); hold on;
plot(time, A12_hat,'r-','LineWidth',1.0);
plot(time, A21_hat,'g-','LineWidth',1.0);
plot(time, A22_hat,'k-','LineWidth',1.0);
plot(time, B1_hat, 'b--','LineWidth',1.0);
plot(time, B2_hat, 'r--','LineWidth',1.0);

% True‐value horizontal lines
plot([0,T],[A_true(1,1),A_true(1,1)], 'b:','LineWidth',1.2);
plot([0,T],[A_true(1,2),A_true(1,2)], 'r:','LineWidth',1.2);
plot([0,T],[A_true(2,1),A_true(2,1)], 'g:','LineWidth',1.2);
plot([0,T],[A_true(2,2),A_true(2,2)], 'k:','LineWidth',1.2);
plot([0,T],[B_true(1),B_true(1)], 'b:','LineWidth',1.2);
plot([0,T],[B_true(2),B_true(2)], 'r:','LineWidth',1.2);

title('Parameter Estimates $\hat{A}_{ij}(t),\;\hat{B}_i(t)$ vs.\ True Values', ...
      'Interpreter','latex');
xlabel('Time (s)','Interpreter','latex');
ylabel('Parameter Value','Interpreter','latex');
legend({ '$\hat{A}_{11}$', '$\hat{A}_{12}$', '$\hat{A}_{21}$', '$\hat{A}_{22}$', ...
         '$\hat{B}_{1}$', '$\hat{B}_{2}$', ...
         '$A_{11}^{\rm true}$', '$A_{12}^{\rm true}$', '$A_{21}^{\rm true}$', '$A_{22}^{\rm true}$', ...
         '$B_{1}^{\rm true}$', '$B_{2}^{\rm true}$' }, ...
       'Interpreter','latex','Location','northeastoutside');
grid on;
