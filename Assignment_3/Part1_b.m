clc;
clear;
close all;

% Plant & adaptation settings
A_true = [-2.15,  0.25; -0.75, -2.00];
B_true = [0; 1.5];

% input signal parameters
inputParams.c = 5.0;
inputParams.d = 0.5;
inputParams.f = 4.0;
inputParams.h = 1.3;
inputParams.i = 3.0;
inputParams.j = 2.7;

% adaptation gains [gamma_A, gamma_B]
adaptGains = [0.0286, 0.0239];

% sigma–modification / bias bounds
biasBounds.A = 0.01;
biasBounds.B = 0.01;
sigmaMod.sigma = 0.05;
sigmaMod.M = 5.0;

% Simulation time
T     = 150;       
dt    = 0.01;      
tspan = 0:dt:T;    

% Initial state [x; x_hat; vec(A_hat); B_hat]
X0 = [ ...
    0; 0;          % x1(0), x2(0)
    0; 0;          % xhat1(0), xhat2(0)
   -2; 0.2;        % Ahat11(0), Ahat12(0)
   -0.6; -1.5;     % Ahat21(0), Ahat22(0)
   -0.1; 1.0       % Bhat1(0), Bhat2(0)
];

% Integrate Lyapunov‐parallel ODE
odefun = @(t, X) lyapunovParallelDynamics( ...
    t, X, A_true, B_true, inputParams, adaptGains, biasBounds, sigmaMod );
[~, X] = ode45(odefun, tspan, X0);

% Pull out histories
x= X(:,1:2);
x_hat = X(:,3:4);
ahat11 = X(:,5);
ahat12 = X(:,6);
ahat21= X(:,7);
ahat22= X(:,8);
bhat1= X(:,9);
bhat2= X(:,10);
e = x - x_hat;

figure('Units','normalized','Position',[.1 .1 .8 .7]);

% true vs estimated
subplot(3,1,1)
plot(tspan, x(:,1),'b-',tspan, x_hat(:,1),'b--','LineWidth',1.2); hold on
plot(tspan, x(:,2),'r-',tspan, x_hat(:,2),'r--','LineWidth',1.2)
title('True States vs. Estimated States','Interpreter','none')
xlabel('Time (s)'); ylabel('States')
legend({'$x_1$ (true)','$\hat x_1$ (est)','$x_2$ (true)','$\hat x_2$ (est)'},...
       'Interpreter','latex','Location','best')
grid on

% estimation errors
subplot(3,1,2)
plot(tspan, e(:,1),'m-','LineWidth',1.2); hold on
plot(tspan, e(:,2),'c-','LineWidth',1.2)
title('State Estimation Error','Interpreter','none')
xlabel('Time (s)'); ylabel('Error')
legend({'$e_1$','$e_2$'},'Interpreter','latex','Location','best')
grid on

% parameter trajectories vs true
subplot(3,1,3)
plot(tspan, ahat11,'b-',tspan, ahat12,'r-',tspan, ahat21,'g-',...
     tspan, ahat22,'k-',tspan, bhat1,'b--',tspan, bhat2,'r--','LineWidth',1)
hold on
plot([0 T],[A_true(1,1) A_true(1,1)],'b:');
plot([0 T],[A_true(1,2) A_true(1,2)],'r:');
plot([0 T],[A_true(2,1) A_true(2,1)],'g:');
plot([0 T],[A_true(2,2) A_true(2,2)],'k:');
plot([0 T],[B_true(1) B_true(1)], 'b:');
plot([0 T],[B_true(2) B_true(2)], 'r:');
title('Parameter Estimates vs True','Interpreter','latex')
xlabel('Time (s)'); ylabel('Parameters')
legend({ ...
  '$\hat A_{11}$','$\hat A_{12}$','$\hat A_{21}$','$\hat A_{22}$', ...
  '$\hat B_{1}$','$\hat B_{2}$', ...
  '$A_{11}^{\rm true}$','$A_{12}^{\rm true}$','$A_{21}^{\rm true}$','$A_{22}^{\rm true}$', ...
  '$B_{1}^{\rm true}$','$B_{2}^{\rm true}$' },...
  'Interpreter','latex','Location','northeastoutside')
grid on

A_final = [ahat11(end), ahat12(end); ahat21(end), ahat22(end)];
B_final = [bhat1(end); bhat2(end)];
fprintf('\n=== Final parameter estimates at T=%.2f ===\n', T);
fprintf('A_hat = [ %.3f   %.3f ;\n          %.3f   %.3f ]\n', ...
        A_final(1,1), A_final(1,2), A_final(2,1), A_final(2,2));
fprintf('A_true = [ %.3f   %.3f ;\n           %.3f   %.3f ]\n', ...
        A_true(1,1),  A_true(1,2),  A_true(2,1),  A_true(2,2));
fprintf('||A_hat − A_true||_F = %.6f\n', norm(A_final - A_true, 'fro'));
fprintf('\nB_hat = [ %.3f ; %.3f ]\n', B_final(1), B_final(2));
fprintf('B_true = [ %.3f ; %.3f ]\n', B_true(1),  B_true(2));
fprintf('||B_hat − B_true||_2 = %.6f\n', norm(B_final - B_true, 2));
fprintf('===============================================\n\n');


function delta = sigmaDiscontinuity(thetaHat, s)
    absTh = abs(thetaHat);
    if absTh < s.M
        delta = 0;
    elseif  absTh <= 2*s.M
        delta = s.sigma*(absTh/s.M - 1);
    else
        delta = s.sigma;
    end
end

% --- ODE function ---
function dX = lyapunovParallelDynamics( ...
    t, X, A, B, inp, gam, bias, sig)

    % unpack state and estimates
    x1   = X(1);  
    x2   = X(2);
    x1h  = X(3);
    x2h  = X(4);
    a11h = X(5);
    a12h = X(6);
    a21h = X(7);
    a22h = X(8);
    b1h  = X(9);
    b2h  = X(10);

    % input + bias
    u = inp.c*sin(inp.d*t) + inp.f*sin(inp.h*t) + inp.i*sin(inp.j*t);
    omega = [ bias.A; bias.B*sin(inp.h*t) ];

    % estimation errors
    e1 = x1 - x1h;
    e2 = x2 - x2h;

    % allocate
    dX = zeros(10,1);

    % true plant
    dX(1) = A(1,1)*x1 + A(1,2)*x2 + B(1)*u + omega(1);
    dX(2) = A(2,1)*x1 + A(2,2)*x2 + B(2)*u + omega(2);

    % observer
    dX(3) = a11h*x1h + a12h*x2h + b1h*u;
    dX(4) = a21h*x1h + a22h*x2h + b2h*u;

    % sigma–modification for each estimate
    d5  = sigmaDiscontinuity(a11h, sig);
    d6  = sigmaDiscontinuity(a12h, sig);
    d7  = sigmaDiscontinuity(a21h, sig);
    d8  = sigmaDiscontinuity(a22h, sig);
    d9  = sigmaDiscontinuity(b1h, sig);
    d10 = sigmaDiscontinuity(b2h, sig);

    % parallel adaptation
    dX(5)  = gam(1)*e1*x1h - gam(1)*d5 * a11h;
    dX(6)  = gam(1)*e1*x2h - gam(1)*d6 * a12h;
    dX(7)  = gam(1)*e2*x1h - gam(1)*d7 * a21h;
    dX(8)  = gam(1)*e2*x2h - gam(1)*d8 * a22h;
    dX(9)  = gam(2)*e1*u - gam(2)*d9 * b1h;
    dX(10) = gam(2)*e2*u - gam(2)*d10* b2h;

    % projection
    if a11h < -3 || a11h > -1
        dX(5)  = 0;
    end
    if b2h  < 1
        dX(10) = 0;
    end
end

