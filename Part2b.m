clc;
clear;
close all;

%%%%%%%%%%%%% 1) TRUE SYSTEM PARAMETERS  %%%%%%%%%%%%%
m_true = 0.75;
L_true = 1.25;
c_true = 0.15;
g      = 9.81;

% Convert to second-order form: q'' + a1 q' + a2 q = b0 u
a1_true = c_true / (m_true * L_true^2);
a2_true = g / L_true;
b0_true = 1 / (m_true * L_true^2);

%%%%%%%%%%%%% 2) SIMULATE "TRUE" SYSTEM  %%%%%%%%%%%%%
% ODE: x1=q, x2=qdot,  => q''= -a1 qdot - a2 q + b0 u(t)
Tend = 10;
x0   = [0.1; 0];  
[tCont, xCont] = ode45(@(t,x) odeSystem(t, x, a1_true, a2_true, b0_true), [0 Tend], x0);

% True outputs
qCont    = xCont(:,1);
qdotCont = xCont(:,2);
uCont    = input_u(tCont);

%%%%%%%%%%%%% 3) SAMPLE AT Ts=0.1  %%%%%%%%%%%%%
Ts          = 0.1;
tSamples    = 0:Ts:Tend; 
qSamples    = interp1(tCont, qCont,    tSamples);
uSamples    = interp1(tCont, uCont,    tSamples);
qdotSamples = interp1(tCont, qdotCont, tSamples);

%%%%%%%%%%%%% 4) APPLY THE FILTERED REGRESSION  %%%%%%%%%%%%%
% We'll feed qSamples and uSamples through the same 2nd-order filter
[tFilt_q, y0_q, y1_q, y2_q] = filter_signals(tSamples, qSamples);
[tFilt_u, u0_u, ~,    ~   ] = filter_signals(tSamples, uSamples);

% Interpolate both sets onto a single time vector (use tFilt_q)
tFilter = tFilt_q;  
y0S = interp1(tFilt_q, y0_q, tFilter);
y1S = interp1(tFilt_q, y1_q, tFilter);
y2S = interp1(tFilt_q, y2_q, tFilter);
u0S = interp1(tFilt_u, u0_u, tFilter);

% We want:  y2(t) = a1*(-y1) + a2*(-y0) + b0*(+u0)
% So we define:
Yvec = y2S(:);
Zeta = [-y1S(:), -y0S(:), u0S(:)]; 
theta_est = (Zeta'*Zeta)\(Zeta'*Yvec);

a1_est = theta_est(1);
a2_est = theta_est(2);
b0_est = theta_est(3);

% Recover m, L, c from a1,a2,b0
L_est = g / a2_est;
m_est = 1 / (b0_est * L_est^2);
c_est = a1_est * m_est * L_est^2;

%%%%%%%%%%%%% 5) DISPLAY AND COMMENT  %%%%%%%%%%%%%
fprintf('\nTrue vs Estimated parameters (Part 2b):\n');
fprintf(' a1 = c/(mL^2):  %.5f vs %.5f\n', a1_true, a1_est);
fprintf(' a2 = g/L:       %.5f vs %.5f\n', a2_true, a2_est);
fprintf(' b0 = 1/(mL^2):  %.5f vs %.5f\n', b0_true, b0_est);
fprintf(' m :             %.5f vs %.5f\n', m_true, m_est);
fprintf(' L :             %.5f vs %.5f\n', L_true, L_est);
fprintf(' c :             %.5f vs %.5f\n', c_true, c_est);

%%%%%%%%%%%%% 6) RESIMULATE WITH ESTIMATED PARAMETERS  %%%%%%%%%%%%%
[tEst, xEst] = ode45(@(t,x) odeSystem(t,x,a1_est,a2_est,b0_est), [0 Tend], x0);
qEst    = xEst(:,1);
qdotEst = xEst(:,2);

%%%%%%%%%%%%% 7) PLOT THE RESULTS  %%%%%%%%%%%%%
figure('Name','Part2b Results','NumberTitle','off');

% (A) q(t) vs q_hat(t)
subplot(3,1,1);
plot(tCont, qCont,'b', tEst, qEst, 'r--','LineWidth',1.5);
xlabel('t (sec)');
ylabel('q(t)');
grid on;
legend('True','Estimated','Location','best');
title('Pendulum angle q(t) - True vs Identified');

% (B) qdot(t) vs qdot_hat(t) for reference
subplot(3,1,2);
plot(tCont, qdotCont,'b', tEst, qdotEst,'r--','LineWidth',1.5);
xlabel('t (sec)');
ylabel('qdot(t)');
grid on;
legend('True','Estimated','Location','best');
title('Pendulum angular velocity');

% (C) Error eq(t) = q(t) - qhat(t)
eq = interp1(tCont, qCont, tEst) - qEst;

MSE_q = mean(eq.^2);
fprintf('\nMean Squared Error of q: %.6f\n', MSE_q);

subplot(3,1,3);
plot(tEst, eq, 'k','LineWidth',1.5);
xlabel('t (sec)');
ylabel('e_q(t)');
title('$e_q(t) = q(t) - \hat{q}(t)$', 'Interpreter', 'latex');
grid on;


%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%

function dx = odeSystem(t, x, a1, a2, b0)
    % 2nd-order system: q'' = -a1 qdot - a2 q + b0 u(t)
    q    = x(1);
    qdot = x(2);
    dx = [ qdot; -a1*qdot - a2*q + b0*input_u(t) ];
end

function u = input_u(t)
    A0 = 4.0;  w = 2.0;
    u  = A0.*sin(w.*t);
end

function [tFilt,z0,z1,z2] = filter_signals(tIn, yIn)
    % Realize 1/(s^2+3s+2) by ODE: z'' + 3z' + 2z = y(t)
    % Then z0 = z(t), z1 = z'(t), and z2(t) = y(t)-3z'(t)-2z(t) = s^2/Λ(s)*y(t).
    z_init = [0;0];
    [tFilt, zSol] = ode45(@(t,z) ...
        filter_ode(t,z,tIn,yIn), [tIn(1), tIn(end)], z_init);
    z0 = zSol(:,1); 
    z1 = zSol(:,2);
    z2 = zeros(size(z1));
    for k=1:length(tFilt)
        yVal = interp1(tIn,yIn,tFilt(k));
        z2(k) = yVal - 3*z1(k) - 2*z0(k);
    end
end

function dz = filter_ode(t,z,tIn,yIn)
    % z(1)'=z(2),  z(2)'= y(t)-3z(2) - 2z(1).
    val = interp1(tIn,yIn,t);
    dz = [z(2); val - 3*z(2) - 2*z(1)];
end
