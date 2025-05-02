clc;
clear;
close all;

%     Plant: r¨ = –a1 ṙ – a2 sin(r) + a3 ṙ² sin(2r) + b u
%     Goal : r(t) : 0  →  π/10  →  0   within 20 s   (d(t)=0)

P.a1 = 1.315;           % [ 1/s ]
P.a2 = 0.725;           % [ 1/s² ]
P.a3 = 0.225;           % [ 1/s³ ]
P.b  = 1.175;           % control effectiveness

r_peak  =  pi/10;       % desired maximum roll angle  [rad]
Tfinal  = 20;           % reach 0 → r_peak → 0 in 20 s
r_d     = @(t) r_peak * sin(pi*t/Tfinal).^2;     % C² reference (zero ends)

C.phi0   = 0.50;        % initial funnel width
C.phiInf = 0.01;        % final funnel width
C.lambda = 0.5;         % funnel convergence rate
C.rho    = 1;           % gain in 2-nd step
C.k1     = 2;           % log-barrier gain (step-1)
C.k2     = 2;           % log-barrier gain (step-2)

dt     = 1e-4;                          % integration step  [s]
tEnd   = 20;                            % simulate 20 s
tSpan  = 0:dt:tEnd;
x0     = [0 ; 0];                       % r(0)=0 , ṙ(0)=0

[t,x] = ode45(@(t,x) plantODE(t,x,P,C,r_d), tSpan, x0);
r     = x(:,1);
r_dot = x(:,2);

u = arrayfun(@(ti,ri,rdi) controlLaw(ti,ri,rdi,C,r_d), t, r, r_dot);

figure;
plot(t, r          , 'b' , 'LineWidth',1.5); hold on
plot(t, r_d(t)     , 'r--', 'LineWidth',1.5);
grid on
title('r(t) and r_d(t)  –  closed-loop response');
xlabel('time  [s]');  ylabel('roll angle  [rad]');
legend('r(t)','r_d(t)','Location','Best');

figure;
plot(t, u, 'LineWidth',1.4); grid on
title('Control input  u(t)');
xlabel('time  [s]');  ylabel('u(t)');

fprintf('Reach r_peak = %.3f rad and return to 0 in %.0f s.\n', ...
        r_peak, Tfinal);

function dx = plantODE(t,x,P,C,r_d_fun)
    r     = x(1);
    r_dot = x(2);

    u     = controlLaw(t, r, r_dot, C, r_d_fun);

    r_ddot = -P.a1*r_dot ...
             -P.a2*sin(r) ...
             +P.a3*r_dot^2 * sin(2*r) ...
             +P.b * u;         % disturbance d(t)=0

    dx = [ r_dot ;
           r_ddot ];
end

function u = controlLaw(t, r, r_dot, C, r_d_fun)
    rd  = r_d_fun(t);

    phi = (C.phi0 - C.phiInf)*exp(-C.lambda*t) + C.phiInf;
    z1  = (r - rd) / phi;
    a   = -C.k1 * log((1+z1)/(1-z1));

    z2  = (r_dot - a) / C.rho;
    u   = -C.k2 * log((1+z2)/(1-z2));
end
