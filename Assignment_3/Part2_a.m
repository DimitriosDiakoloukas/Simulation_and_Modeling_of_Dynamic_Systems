clc;
clear;
close all;

% True system and simulation settings
theta_true = [1.5; 1.0];      % [θ1; θ2] in [0.5, 2]
Ts = 0.01;            % sampling interval
N = 10000;       % number of samples
t = (0:N-1)' * Ts;   % time vector

% input signal u(t): sinusoid with white noise
u = 0.5*sin(2*pi*0.5*t) + 0.1*randn(N,1);

% simulate x(t) via Euler integration of x_dot = f_true(x,u)
x = zeros(N,1);
for k = 1:N-1
    f_val = - x(k)^3 ...
            + theta_true(1)*tanh( x(k) ) ...
            + theta_true(2)/(1 + x(k)^2 ) ...
            + u(k);
    x(k+1) = x(k) + Ts * f_val;
end

% approximate derivative ẋ ≈ Δx/Ts
dx = diff(x) / Ts;
x_mid = x(1:end-1);
u_mid = u(1:end-1);

% Train Test split
train_frac = 0.5;
N_train = floor(length(dx)*train_frac);

x_tr = x_mid(1:N_train);
u_tr = u_mid(1:N_train);
y_tr = dx(1:N_train) - u_mid(1:N_train);   % target = x_dot - u

x_te = x_mid(N_train+1:end);
u_te = u_mid(N_train+1:end);
y_te = dx(N_train+1:end) - u_mid(N_train+1:end);

basis{1}.name = 'Polynomial (x,x^2,x^3)';
basis{1}.Phi = @(z)[ z, z.^2, z.^3 ];

basis{2}.name = 'Mixed (tanh,1/(1+x^2))';
basis{2}.Phi = @(z)[ tanh(z), 1./(1+z.^2) ];

basis{3}.name = 'Full (poly + mixed)';
basis{3}.Phi = @(z)[ z, z.^2, z.^3, tanh(z), 1./(1+z.^2) ];

% RLS parameters
lambda = 1.0;       % forgetting factor
P0     = 1e3;       % initial P matrix diagonal

results = cell(length(basis),1);

for i = 1:length(basis)
    % Build of my regression matrix for train and test
    Phi_tr = basis{i}.Phi(x_tr);
    Phi_te = basis{i}.Phi(x_te);
    
    npar = size(Phi_tr,2);
    theta = zeros(npar,1);
    P = eye(npar)*P0;
    
    theta_hist = zeros(N_train, npar);
    
    % Online RLS loop
    for k = 1:N_train
        phi_k = Phi_tr(k,:)';      % column vector
        yk = y_tr(k);
        
        % RLS gain
        K = (P * phi_k) / ( lambda + phi_k'*P*phi_k );
        
        % update estimate
        theta = theta + K * ( yk - phi_k'*theta );
        
        % update covariance
        P = (P - K*phi_k'*P)/lambda;
        
        theta_hist(k,:) = theta';
    end
    
    % Test‐set MSE
    y_pred = Phi_te * theta;
    MSE= mean( (y_te - y_pred).^2 );
    
    results{i}.name = basis{i}.name;
    results{i}.MSE  = MSE;
    results{i}.theta_hist = theta_hist;
    results{i}.theta_est  = theta;
    
    % Plot convergence
    figure('Name',basis{i}.name,'NumberTitle','off');
    plot(1:N_train, theta_hist, 'LineWidth',1.2);
    xlabel('Iteration k');
    ylabel('Parameter value');
    title(['Parameter Convergence – ' basis{i}.name]);
    legend(arrayfun(@(n) sprintf('\\theta_{%d}',n), 1:npar,'uni',false), ...
           'Location','best');
    grid on;
end

fprintf('\nModel Structure\t\t\t Test MSE\n');
fprintf('---------------------------------------------\n');
for i=1:length(results)
    fprintf('%-30s %12.4e\n', results{i}.name, results{i}.MSE);
end
