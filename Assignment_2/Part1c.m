clc;
clear;
close all;

%  === True system and simulation data ===
m = 1.315;                % mass        [kg]
b = 0.225;                % damping     [N·s/m]
k = 0.725;                % stiffness   [N/m]

A = [ 0   1 ;
     -k/m  -b/m ];        % state matrix
B = [0 ; 1/m];            % input matrix

% Input signal u(t) = c*sin(d*t)+f
c = 2.5;  d = 1;  f = 0;

% Measurement-noise model  h(t) = h0*sin(2*pi*f0*t)
h0 = 0.25;      
f0 = 20;       

tStep = 0.01;
tSpan = 0 : tStep : 100;

% === Parallel configuration with noise ===
gamma1 = 0.059;    % gain for a21_hat
gamma2 = 0.064;    % gain for a22_hat
gamma3 = 0.0191;    % gain for b2_hat

x0_par = [0;0;0;0;  -0.01;-0.01;0.01];  % initial state & estimates

P = struct('m',m,'b',b,'k',k,'A',A,'B',B, ...
           'c',c,'d',d,'f',f, ...
           'h0',h0,'f0',f0, ...
           'gamma1',gamma1,'gamma2',gamma2,'gamma3',gamma3);

[t,X] = ode45(@(t,x) lyap_parallel_noise(t,x,P), tSpan, x0_par);

x1 = X(:,1);   x2 = X(:,2);
x1_hat = X(:,3);   x2_hat = X(:,4);
a21_hat = X(:,5);  a22_hat = X(:,6);  b2_hat = X(:,7);

m_hat = 1 ./ b2_hat;
b_hat = -a22_hat .* m_hat;
k_hat = -a21_hat .* m_hat;

e1 = x1 - x1_hat;

% ----- Plots (parallel) -----
figure; plot(t,x1,'r',t,x1_hat,'b'); grid on;
title('x1 and x1_hat  (Parallel - noise)');
xlabel('time in sec'); ylabel('Displacement');
legend('x1','x1_hat');

figure; plot(t,e1,'k'); grid on;
title('Estimation error e1  (Parallel - noise)');
xlabel('time in sec'); ylabel('e1');

mask = false(size(t));  

for i = 1:length(t)
    if t(i) >= 0.5
        mask(i) = true;
    end
end
t_plot = t(mask);
figure;
subplot(3,1,1); plot(t_plot,m_hat(mask),'LineWidth',1.5); hold on; yline(m,'--r');
ylabel('m_hat(t)'); grid on; title('Estimated parameters – Parallel (noise)');
subplot(3,1,2); plot(t_plot,b_hat(mask),'LineWidth',1.5); hold on; yline(b,'--r');
ylabel('b_hat(t)'); grid on;
subplot(3,1,3); plot(t_plot,k_hat(mask),'LineWidth',1.5); hold on; yline(k,'--r');
ylabel('k_hat(t)'); xlabel('time in sec'); grid on;

fprintf('\n=== Final parameter estimates  (Parallel – noise) ===\n');
fprintf('Mass      : true %.4f,  estimated %.4f\n', m, m_hat(end));
fprintf('Damping   : true %.4f,  estimated %.4f\n', b, b_hat(end));
fprintf('Stiffness : true %.4f,  estimated %.4f\n', k, k_hat(end));

%  ==== Series-Parallel configuration with noise === 
thetam = [0.2 0.2 0.2 0.2];

Psp            = P;
Psp.gamma1     = 0.06;
Psp.gamma2     = 0.057;
Psp.gamma3     = 1.5;
Psp.thetam     = thetam;

x0_sp = [0;0;0;0;  0.01;0.01;0.01];

[t,X] = ode45(@(t,x) lyap_series_parallel_noise(t,x,Psp), tSpan, x0_sp);

x1 = X(:,1);   x2 = X(:,2);
x1_hat = X(:,3);   x2_hat = X(:,4);
a21_hat = X(:,5);  a22_hat = X(:,6);  b2_hat = X(:,7);

m_hat = 1 ./ b2_hat;
b_hat = -a22_hat .* m_hat;
k_hat = -a21_hat .* m_hat;

e1 = x1 - x1_hat;

mask = false(size(t));  

for i = 1:length(t)
    if t(i) >= 0.5
        mask(i) = true;
    end
end
t_plot = t(mask);
smooth = @(v) movmean(v,500);

figure; plot(t,x1,'r',t,x1_hat,'b'); grid on;
title('x1 and x1_hat  (Series-Parallel - noise)');
xlabel('time in sec'); ylabel('Displacement');
legend('x1','x1_hat');

figure; plot(t,e1,'k'); grid on;
title('Estimation error e1  (Series-Parallel - noise)');
xlabel('time in sec'); ylabel('e1');

figure;
subplot(3,1,1); plot(t_plot, smooth(m_hat(mask)),'LineWidth',1.5); hold on; yline(m,'--r');
ylabel('m_hat(t)'); grid on; title('Estimated parameters – Series-Parallel (noise)');
subplot(3,1,2); plot(t_plot, smooth(b_hat(mask)),'LineWidth',1.5); hold on; yline(b,'--r');
ylabel('b_hat(t)'); grid on;
subplot(3,1,3); plot(t_plot, smooth(k_hat(mask)),'LineWidth',1.5); hold on; yline(k,'--r');
ylabel('k_hat(t)'); xlabel('time in sec'); grid on;

fprintf('\n=== Final parameter estimates  (Series-Parallel – noise) ===\n');
fprintf('Mass      : true %.4f  –  est %.4f\n', m, m_hat(end));
fprintf('Damping   : true %.4f  –  est %.4f\n', b, b_hat(end));
fprintf('Stiffness : true %.4f  –  est %.4f\n', k, k_hat(end));

% ======================  Dynamics ===============================
function dx = lyap_parallel_noise(t,x,P)
    u = P.c * sin(P.d*t) + P.f;
    h = P.h0 * sin(2*pi*P.f0*t);
    e1 = (x(1)+h) - x(3);
    e2 =  x(2)    - x(4);

    dx = zeros(7,1);
    dx(1:2) = P.A*[x(1);x(2)] + P.B*u;
    dx(3)   = x(4);
    dx(4)   = x(5)*x(3) + x(6)*x(4) + x(7)*u;
    dx(5)   = P.gamma1*e2*x(3);
    dx(6)   = P.gamma2*e2*x(4);
    dx(7)   = P.gamma3*e2*u;
end

function dx = lyap_series_parallel_noise(t,x,P)
    u = P.c * sin(P.d*t) + P.f;
    h = P.h0 * sin(2*pi*P.f0*t);
    e1 = (x(1)+h) - x(3);
    e2 =  x(2)    - x(4);

    dx = zeros(7,1);
    dx(1) = x(2);
    dx(2) = P.A(2,:)*[x(1);x(2)] + P.B(2)*u;
    dx(3) = x(2) + P.thetam(1)*e1 + P.thetam(2)*e2;
    dx(4) = x(5)*(x(1)+h) + x(6)*x(2) + x(7)*u + ...
            P.thetam(3)*e1 + P.thetam(4)*e2;
    dx(5) = P.gamma1*e2*(x(1)+h);
    dx(6) = P.gamma2*e2*x(2);
    dx(7) = P.gamma3*e2*u;
end
