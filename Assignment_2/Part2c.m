clc;
clear;
close all;

% ============  Simulation & reference ===============================
dt    = 1e-4;                 % integration step [s]
tSpan = 0:dt:20;              % simulation horizon

r_d_peak = pi/10;  Tref = 20;                     % C² reference
r_d      = @(t) r_d_peak * sin(pi*t/Tref).^2;     % 0→π/10→0 in 20 s

d_fun = @(t) 0.15 * sin(0.5*t);                   % disturbance

a1 = 1.315;   a2 = 0.725;   a3 = 0.225;   b  = 1.175;

phi0 = 0.50;   phiInf = 0.01;   lambda = 0.5;
rho  = 1;      k1      = 2;     k2     = 2;

x0 = [0 0];     %  r(0)=0 , ṙ(0)=0
[t,x] = ode45(@(t,x) plant_dyn(t,x, ...
                    a1,a2,a3,b,phi0,phiInf,lambda,rho,k1,k2,r_d,d_fun), ...
              tSpan,x0);

r      = x(:,1);     r_dot = x(:,2);
r_d_vec = arrayfun(r_d,t);

u = zeros(size(t));
for k = 1:numel(t)
    phi_k = (phi0 - phiInf)*exp(-lambda*t(k)) + phiInf;
    z1    = (r(k) - r_d_vec(k))/phi_k;
    a_k   = -k1*log((1+z1)/(1-z1));
    z2    = (r_dot(k) - a_k)/rho;
    u(k)  = -k2*log((1+z2)/(1-z2));
end

figure; plot(t,u,'LineWidth',1.5); grid on;
xlabel('time in sec'); ylabel('u(t)'); title('Control input  u(t)  (Part 2-c)');

% ============  Series–Parallel estimator (with disturbance) ==========
gamma1 = 6;   gamma2 = 0.01;   gamma3 = 0.4;   gamma4 = 0.03;
thetaM1 = 10;  thetaM2 = 10;

N = numel(t);
a1_hat = zeros(1,N); a1_hat(1) = 1.0;
a2_hat = zeros(1,N); a2_hat(1) = 1.0;
a3_hat = zeros(1,N); a3_hat(1) = 0.5;
b_hat  = zeros(1,N); b_hat (1) = 1.0;

r_hat     = zeros(1,N); r_hat(1)     = r(1);
r_dot_hat = zeros(1,N); r_dot_hat(1) = r_dot(1);

for k = 1:N-1
    f  = sin( r(k) );
    g  = r_dot(k)^2 * sin(2*r(k));
    ed = r_dot(k) - r_dot_hat(k);                 % error in ṙ

    a1_hat(k+1) = a1_hat(k) - gamma1*r_dot(k)*ed;
    a2_hat(k+1) = a2_hat(k) - gamma2*f*ed;
    a3_hat(k+1) = a3_hat(k) + gamma3*g*ed;
    b_hat (k+1) = b_hat (k) + gamma4*u(k)*ed;

    r_ddot_hat     = - a1_hat(k)*r_dot(k) - a2_hat(k)*f + a3_hat(k)*g ...
                     + b_hat(k)*u(k) + thetaM2*ed + d_fun(t(k));
    r_dot_hat(k+1) = r_dot_hat(k) + dt*r_ddot_hat;
    r_hat(k+1)     = r_hat(k) + dt*( r_dot_hat(k) ...
                        + thetaM1*( r(k) - r_hat(k) ) );
end

e_r = r - r_hat(:); 
figure; plot(t,r,'b',t,r_hat,'r--','LineWidth',1.5); grid on;
xlabel('time in sec'); ylabel('roll angle  [rad]');
legend('r (true)','r̂ (est)'); title('True vs estimated roll angle  (d≠0)');

figure; plot(t,e_r,'k','LineWidth',1.5); grid on;
xlabel('time in sec'); ylabel('error  r - r̂');
title('Estimation error  (Series–Parallel, d≠0)');

figure;
subplot(4,1,1); plot(t,a1_hat,'LineWidth',1.5); hold on; yline(a1,'--r');
ylabel('a1̂'); grid on; title('Parameter estimates  (d≠0)')
subplot(4,1,2); plot(t,a2_hat,'LineWidth',1.5); hold on; yline(a2,'--r');
ylabel('a2̂'); grid on;
subplot(4,1,3); plot(t,a3_hat,'LineWidth',1.5); hold on; yline(a3,'--r');
ylabel('a3̂'); grid on;
subplot(4,1,4); plot(t,b_hat ,'LineWidth',1.5); hold on; yline(b ,'--r');
ylabel('b̂'); xlabel('time in sec'); grid on;

fprintf('\n=== Final parameter estimates (Series–Parallel, d≠0) ===\n');
fprintf('a1 : true %8.4f   est %8.4f\n', a1, a1_hat(end));
fprintf('a2 : true %8.4f   est %8.4f\n', a2, a2_hat(end));
fprintf('a3 : true %8.4f   est %8.4f\n', a3, a3_hat(end));
fprintf('b  : true %8.4f   est %8.4f\n',  b,  b_hat(end));
fprintf('---------------------------------------------------------\n');

% ====================== Helper function ===============================
function dx = plant_dyn(t,x,a1,a2,a3,b,phi0,phiInf,lambda,rho,k1,k2,r_d,d_fun)
    r     = x(1);               % roll angle
    r_dot = x(2);               % roll-rate
    rd    = r_d(t);             % reference
    d     = d_fun(t);           % disturbance

    phi = (phi0 - phiInf)*exp(-lambda*t) + phiInf;
    z1  = (r - rd)/phi;
    a   = -k1*log((1+z1)/(1-z1));

    z2  = (r_dot - a)/rho;
    u   = -k2*log((1+z2)/(1-z2));

    dx = zeros(2,1);
    dx(1) = r_dot;
    dx(2) = -a1*r_dot - a2*sin(r) + a3*r_dot^2*sin(2*r) ...
            + b*u + d;                        % disturbance added
end
