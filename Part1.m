clc;
clear;
close all;

m = 0.7;   
L = 1.25;        
c = 0.15;        
g = 9.81;      
A0 = 4.0;       
w = 2.0;    

tFinal = 20;      
tspan = [0 tFinal];

x0 = [0.0; 0.0];
options = odeset('MaxStep',1e-3);

[tSol, xSol] = ode45(@(t,x) pendulum_ode(t, x, m, L, c, g, A0, w), ...
                     tspan, x0, options);

figure;
subplot(2,1,1);
plot(tSol, xSol(:,1), 'LineWidth',2); 
xlabel('Time [s]'); 
ylabel('Angle q(t) [rad]');
title('Pendulum Angle vs. Time'); 
grid on;

subplot(2,1,2);
plot(tSol, xSol(:,2), 'LineWidth',2);
xlabel('Time [s]'); 
ylabel('Angular Velocity dq(t)/dt [rad/s]');
title('Pendulum Angular Velocity vs. Time'); 
grid on;

function dx = pendulum_ode(t, x, m, L, c, g, A0, w)
    u = A0 * sin(w * t);
    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = (1/(m*L^2)) * (u - c*x(2) - m*g*L*x(1));
end
