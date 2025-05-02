clc;
clear;
close all;

dt     = 1e-4;                          % integration step  [s]
tSpan  = 0:dt:20;                       % simulation horizon

r_d_val = pi/10;  Tref = 20;            % smooth reference r_d(t)
r_d     = @(t) r_d_val * sin(pi*t/Tref).^2;

a1 = 1.315;   a2 = 0.725;   a3 = 0.225;   b  = 1.175;

phi0 = 0.50;  phiInf = 0.01;  lambda = 0.5;
rho  = 1;     k1      = 2;     k2     = 2;

x0 = [0 0];   % [r(0)  ṙ(0)]
[t,x] = ode45(@(t,x) plant_dynamics(t,x, ...
                    a1,a2,a3,b,phi0,phiInf,lambda,rho,k1,k2,r_d), ...
              tSpan,x0);

r      = x(:,1);    r_dot = x(:,2);    r_d_vec = arrayfun(r_d,t);

u = zeros(size(t));
for k = 1:numel(t)
    phi_k = (phi0 - phiInf)*exp(-lambda*t(k)) + phiInf;
    z1    = (r(k) - r_d_vec(k))/phi_k;
    a_k   = -k1*log((1+z1)/(1-z1));
    z2    = (r_dot(k) - a_k)/rho;
    u(k)  = -k2*log((1+z2)/(1-z2));
end

figure; plot(t,u,'LineWidth',1.5); grid on;
xlabel('time in sec'); ylabel('u(t)');
title('Control input  u(t)  (Part 2-a)');

% === Series–Parallel estimator (part 2-b) =============================
gamma1 = 8.52; gamma2 = 0.0267; gamma3 = 0.01; gamma4 = 0.03; % adaptation gains
thetaM1 = 50;  thetaM2 = 90;                                  % observer gains

a1_hat=zeros(size(t)); a1_hat(1)=1.0;
a2_hat=zeros(size(t)); a2_hat(1)=1.0;
a3_hat=zeros(size(t)); a3_hat(1)=0.5;
b_hat =zeros(size(t)); b_hat(1) =1.0;

r_hat     = zeros(size(t)); r_hat(1)     = r(1);
r_dot_hat = zeros(size(t)); r_dot_hat(1) = r_dot(1);

for k = 1:numel(t)-1
    f = sin( r(k) );
    g = r_dot(k)^2 * sin(2*r(k));
    err_dot = r_dot(k) - r_dot_hat(k);

    a1_hat(k+1) = a1_hat(k) - gamma1*r_dot(k)*err_dot;
    a2_hat(k+1) = a2_hat(k) - gamma2*f*err_dot;
    a3_hat(k+1) = a3_hat(k) + gamma3*g*err_dot;
    b_hat (k+1) = b_hat (k) + gamma4*u(k)*err_dot;

    r_ddot_hat     = -a1_hat(k)*r_dot(k) - a2_hat(k)*f + a3_hat(k)*g ...
                     + b_hat(k)*u(k) + thetaM2*err_dot;
    r_dot_hat(k+1) = r_dot_hat(k) + dt*r_ddot_hat;
    r_hat(k+1)     = r_hat(k) + dt*( r_dot_hat(k) ...
                          + thetaM1*( r(k) - r_hat(k) ) );
end

% === Plots (part 2-b) =================================================
e_r = r - r_hat;
figure; plot(t,r,'b',t,r_hat,'r--','LineWidth',1.5); grid on;
xlabel('time in sec'); ylabel('roll angle  [rad]');
legend('r (true)','r_hat (est)'); title('True vs estimated roll angle');

figure; plot(t,e_r,'k','LineWidth',1.5); grid on;
xlabel('time in sec'); ylabel('error  r - r_hat');
title('Estimation error  (Series–Parallel)');

figure;
subplot(4,1,1); plot(t,a1_hat,'LineWidth',1.5); hold on; yline(a1,'--k');
ylabel('a1̂'); grid on; title('Parameter estimates (Series-Parallel)');
subplot(4,1,2); plot(t,a2_hat,'LineWidth',1.5); hold on; yline(a2,'--k');
ylabel('a2̂'); grid on;
subplot(4,1,3); plot(t,a3_hat,'LineWidth',1.5); hold on; yline(a3,'--k');
ylabel('a3̂'); grid on;
subplot(4,1,4); plot(t,b_hat ,'LineWidth',1.5); hold on; yline(b ,'--k');
ylabel('b̂'); xlabel('time in sec'); grid on;

fprintf('\n=== Final parameter estimates (Series–Parallel) ===\n');
fprintf('a1 : true %8.4f   est %8.4f\n', a1, a1_hat(end));
fprintf('a2 : true %8.4f   est %8.4f\n', a2, a2_hat(end));
fprintf('a3 : true %8.4f   est %8.4f\n', a3, a3_hat(end));
fprintf('b  : true %8.4f   est %8.4f\n',  b,  b_hat(end));
fprintf('---------------------------------------------------\n');

% ======================= Helper function =================================
function dx = plant_dynamics(t,x,a1,a2,a3,b,phi0,phiInf,lambda,rho,k1,k2,r_d)
    r     = x(1);            % roll angle
    r_dot = x(2);            % roll-rate
    rd    = r_d(t);          % reference

    phi   = (phi0 - phiInf)*exp(-lambda*t) + phiInf;
    z1    = (r - rd)/phi;
    a     = -k1*log((1+z1)/(1-z1));

    z2    = (r_dot - a)/rho;
    u     = -k2*log((1+z2)/(1-z2));

    dx = zeros(2,1);
    dx(1) = r_dot;
    dx(2) = -a1*r_dot - a2*sin(r) + a3*r_dot^2*sin(2*r) + b*u;   % d(t)=0
end
