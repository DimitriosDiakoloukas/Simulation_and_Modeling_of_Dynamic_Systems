clc;
clear;
close all;

%%%%%%%%%%%%% 1) TRUE SYSTEM PARAMETERS  %%%%%%%%%%%%%
m_true = 0.75;
L_true = 1.25;
c_true = 0.15;
g      = 9.81;

a1_true = c_true / (m_true * L_true^2);
a2_true = g / L_true;
b0_true = 1 / (m_true * L_true^2);

%%%%%%%%%%%%% 2) SIMULATE "TRUE" SYSTEM  %%%%%%%%%%%%%
Tend = 10;
x0   = [0.1; 0];  
[tCont, xCont] = ode45(@(t,x) odeSystem(t, x, a1_true, a2_true, b0_true), [0 Tend], x0);

qCont    = xCont(:,1);
qdotCont = xCont(:,2);
uCont    = input_u(tCont);

%%%%%%%%%%%%% 3) SAMPLE AT Ts=0.1 (BASELINE)  %%%%%%%%%%%%%
Ts          = 0.1;
tSamples    = 0:Ts:Tend; 
qSamples    = interp1(tCont, qCont,    tSamples);
uSamples    = interp1(tCont, uCont,    tSamples);
qdotSamples = interp1(tCont, qdotCont, tSamples);

%%%%%%%%%%%%% 4) APPLY THE FILTERED REGRESSION  %%%%%%%%%%%%%
[tFilt_q, y0_q, y1_q, y2_q] = filter_signals(tSamples, qSamples);
[tFilt_u, u0_u, ~,    ~   ] = filter_signals(tSamples, uSamples);

tFilter = tFilt_q;  
y0S = interp1(tFilt_q, y0_q, tFilter);
y1S = interp1(tFilt_q, y1_q, tFilter);
y2S = interp1(tFilt_q, y2_q, tFilter);
u0S = interp1(tFilt_u, u0_u, tFilter);

Yvec = y2S(:);
Zeta = [-y1S(:), -y0S(:), u0S(:)]; 
theta_est = (Zeta'*Zeta)\(Zeta'*Yvec);

a1_est = theta_est(1);
a2_est = theta_est(2);
b0_est = theta_est(3);

L_est = g / a2_est;
m_est = 1 / (b0_est * L_est^2);
c_est = a1_est * m_est * L_est^2;

%%%%%%%%%%%%% 5) DISPLAY AND COMMENT  %%%%%%%%%%%%%
fprintf('[Part2b Core Used] Base Ts=0.1: True vs Estimated:\n');
fprintf(' a1 = %.5f vs %.5f\n', a1_true, a1_est);
fprintf(' a2 = %.5f vs %.5f\n', a2_true, a2_est);
fprintf(' b0 = %.5f vs %.5f\n', b0_true, b0_est);
fprintf(' m  = %.5f vs %.5f\n', m_true, m_est);
fprintf(' L  = %.5f vs %.5f\n', L_true, L_est);
fprintf(' c  = %.5f vs %.5f\n', c_true, c_est);

%%%%%%%%%%%%% 6) RESIMULATE WITH ESTIMATED PARAMETERS  %%%%%%%%%%%%%
[tEst, xEst] = ode45(@(t,x) odeSystem(t,x,a1_est,a2_est,b0_est), [0 Tend], x0);
qEst    = xEst(:,1);
qdotEst = xEst(:,2);

%%%%%%%%%%%%% 7) PLOT THE RESULTS  %%%%%%%%%%%%%
figure('Name','Part2b Core Results','NumberTitle','off');
subplot(3,1,1);
plot(tCont, qCont,'b', tEst, qEst, 'r--','LineWidth',1.5);
xlabel('t (sec)');
ylabel('q(t)');
grid on;
legend('True','Estimated','Location','best');
title('Pendulum angle q(t) - True vs Identified');

subplot(3,1,2);
plot(tCont, qdotCont,'b', tEst, qdotEst,'r--','LineWidth',1.5);
xlabel('t (sec)');
ylabel('qdot(t)');
grid on;
legend('True','Estimated','Location','best');
title('Pendulum angular velocity');

eq = interp1(tCont, qCont, tEst) - qEst;
subplot(3,1,3);
plot(tEst, eq, 'k','LineWidth',1.5);
xlabel('t (sec)');
ylabel('e_q(t)');
title('$e_q(t) = q(t) - \hat{q}(t)$','Interpreter','latex');
grid on;

%%%%%%%%%%%%% [PART 3(b)] Vary Ts in a loop %%%%%%%%%%%%%
disp(' ');
disp('=== Part 3(b): Vary Ts and see the param estimation error ===');

Ts_list = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0];
err_a1 = zeros(size(Ts_list));
err_a2 = zeros(size(Ts_list));
err_b0 = zeros(size(Ts_list));
err_m  = zeros(size(Ts_list));
err_L  = zeros(size(Ts_list));
err_c  = zeros(size(Ts_list));

for i = 1:length(Ts_list)
    Ts_curr = Ts_list(i);
    tS2 = 0:Ts_curr:Tend;
    qS2 = interp1(tCont, qCont, tS2);
    uS2 = interp1(tCont, uCont, tS2);
    
    [tFq2, y0q2, y1q2, y2q2] = filter_signals(tS2, qS2);
    [tFu2, u0u2, ~,   ~   ]  = filter_signals(tS2, uS2);
    
    tCommon = tFq2;
    y0_2 = interp1(tFq2, y0q2, tCommon);
    y1_2 = interp1(tFq2, y1q2, tCommon);
    y2_2 = interp1(tFq2, y2q2, tCommon);
    u0_2 = interp1(tFu2, u0u2, tCommon);
    
    Yv2 = y2_2(:);
    Z2  = [-y1_2(:), -y0_2(:), u0_2(:)];
    th2 = (Z2'*Z2)\(Z2'*Yv2);
    
    a1_2 = th2(1);
    a2_2 = th2(2);
    b0_2 = th2(3);
    
    % Compute estimated m, L, c
    L_est_i = g / a2_2;
    m_est_i = 1 / (b0_2 * L_est_i^2);
    c_est_i = a1_2 * m_est_i * L_est_i^2;
    
    err_a1(i) = a1_2 - a1_true;
    err_a2(i) = a2_2 - a2_true;
    err_b0(i) = b0_2 - b0_true;
    
    err_m(i)  = m_est_i - m_true;
    err_L(i)  = L_est_i - L_true;
    err_c(i)  = c_est_i - c_true;
end

figure('Name','Part3b: Vary Ts','NumberTitle','off');
subplot(3,1,1);
plot(Ts_list, abs(err_a1),'o-','LineWidth',1.5);
grid on;
ylabel('Error in a1');
title('3(b): Param. Error vs Ts');

subplot(3,1,2);
plot(Ts_list, abs(err_a2),'o-','LineWidth',1.5);
grid on;
ylabel('Error in a2');

subplot(3,1,3);
plot(Ts_list, abs(err_b0),'o-','LineWidth',1.5);
grid on;
ylabel('Error in b0');
xlabel('Sampling Period Ts');

figure('Name','Errors in m, L, c vs Ts','NumberTitle','off');
subplot(3,1,1);
plot(Ts_list, abs(err_m),'o-','LineWidth',1.5);
grid on;
ylabel('Error in m');
title('Error in m, L, c vs Ts');

subplot(3,1,2);
plot(Ts_list, abs(err_L),'o-','LineWidth',1.5);
grid on;
ylabel('Error in L');

subplot(3,1,3);
plot(Ts_list, abs(err_c),'o-','LineWidth',1.5);
grid on;
xlabel('Sampling Period Ts');
ylabel('Error in c');


%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%

function dx = odeSystem(t, x, a1, a2, b0)
    q    = x(1);
    qdot = x(2);
    dx = [qdot; -a1*qdot - a2*q + b0*input_u(t)];
end

function u = input_u(t)
    A0 = 4.0;
    w = 2.0;
    u  = A0.*sin(w.*t);
end

function [tFilt,z0,z1,z2] = filter_signals(tIn, yIn)
    z_init = [0;0];
    [tFilt,zSol] = ode45(@(t,z) filter_ode(t,z,tIn,yIn), [tIn(1),tIn(end)], z_init);
    z0 = zSol(:,1);
    z1 = zSol(:,2);
    z2 = zeros(size(z1));
    for k=1:length(tFilt)
        val = interp1(tIn,yIn,tFilt(k));
        z2(k) = val - 3*z1(k) - 2*z0(k);
    end
end

function dz = filter_ode(t,z,tIn,yIn)
    val = interp1(tIn,yIn,t);
    dz = [ z(2); val - 3*z(2) - 2*z(1) ];
end
