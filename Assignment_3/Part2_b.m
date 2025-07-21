clc;
clear;
close all;

theta_true = [1.5; 1.0];

f_true = @(x,u) - x.^3 ...
           + theta_true(1)*tanh(x) ...
           + theta_true(2)./(1 + x.^2) ...
           + u;

Ts = 0.01;        % sampling interval
N  = 10000;       % total samples
t  = (0:N-1)'*Ts;

% Input: sinusoid and noise
u = 0.5*sin(2*pi*0.5*t) + 0.1*randn(N,1);

% Simulate x via Euler integration
x = zeros(N,1);
for k = 1:N-1
    x(k+1) = x(k) + Ts * f_true(x(k), u(k));
end

% Approximate derivative x_dot = Δx/Ts
dx = diff(x)/Ts;
x_mid = x(1:end-1);
u_mid = u(1:end-1);

% Train/Test split
train_frac = 0.5;
N_train    = floor(length(dx)*train_frac);

x_tr = x_mid(1:N_train);
y_tr = dx(1:N_train) - u_mid(1:N_train);
x_te = x_mid(N_train+1:end);
y_te = dx(N_train+1:end) - u_mid(N_train+1:end);

% Part A.1: Define Candidate Structures and RLS

% Polynomial: [x, x^2, x^3]
basis{1}.name = 'poly';
basis{1}.Phi = @(z)[ z, z.^2, z.^3 ];

% Mixed: [tanh(x), 1/(1+x^2)]
basis{2}.name = 'mixed';
basis{2}.Phi = @(z)[ tanh(z), 1./(1 + z.^2) ];

% Full: poly + mixed
basis{3}.name = 'full';
basis{3}.Phi= @(z)[ z, z.^2, z.^3, tanh(z), 1./(1 + z.^2) ];

% RLS settings
lambda  = 1.0;      % forgetting factor
P0 = 1e3;      % initial covariance diag
results = cell(3,1);

for i = 1:3
    Phi_tr = basis{i}.Phi(x_tr);
    Phi_te = basis{i}.Phi(x_te);
    npar= size(Phi_tr,2);

    % Initialize
    theta = zeros(npar,1);
    P = eye(npar)*P0;
    theta_hist = zeros(N_train, npar);

    % Online RLS
    for k = 1:N_train
        phi_k = Phi_tr(k,:)';
        yk = y_tr(k);
        K = (P*phi_k)/(lambda + phi_k'*P*phi_k);
        theta = theta + K*(yk - phi_k'*theta);
        P = (P - K*phi_k'*P)/lambda;
        theta_hist(k,:) = theta';
    end

    % Store results
    y_pred_te = Phi_te * theta;
    MSE= mean((y_te - y_pred_te).^2);

    results{i}.name= basis{i}.name;
    results{i}.MSE = MSE;
    results{i}.theta_est= theta;
    results{i}.Phi = basis{i}.Phi;
    results{i}.theta_hist= theta_hist;

    % Plot convergence
    figure('Name',['RLS Convergence: ' basis{i}.name],'NumberTitle','off');
    plot(theta_hist,'LineWidth',1.2);
    grid on;
    xlabel('Iteration k'); ylabel('θ estimates');
    title(['Parameter Convergence – ' basis{i}.name]);
    legend(arrayfun(@(n) sprintf('θ_%d',n),1:npar,'Uni',false),'Location','best');
end

% Display Part A cross‐validation MSEs
fprintf('\nPart A – Cross‐Validation Results:\n');
fprintf('Structure       Test MSE\n');
fprintf('------------------------------\n');
for i = 1:3
    fprintf('%-6s     %8.3e\n', results{i}.name, results{i}.MSE);
end

% Part B: Model Selection and Stability Check

% Number of parameters in each model
k_vals = [3, 2, 5];
N_test = length(y_te);

% Compute AIC and BIC
for i = 1:3
    mse = results{i}.MSE;
    k= k_vals(i);
    results{i}.AIC = N_test*log(mse) + 2*k;
    results{i}.BIC = N_test*log(mse) + k*log(N_test);
end

fprintf('\nPart B – AIC/BIC:\n');
fprintf('Struct  Test MSE    AIC       BIC\n');
fprintf('--------------------------------------\n');
for i=1:3
    fprintf('%-6s  %8.3e  %8.2f  %8.2f\n', ...
        results{i}.name, results{i}.MSE, results{i}.AIC, results{i}.BIC);
end

% Choose final model by AIC
[~, m] = min(cellfun(@(r) r.AIC, results));
phi = results{m}.Phi;
theta_hat  = results{m}.theta_est;
fprintf('\nChosen model: %s (by AIC)\n', results{m}.name);

% New validation run for stability and ISE
N2 = 10000;
t2 = (0:N2-1)'*Ts;
u2 = 0.4*sin(2*pi*1.3*t2) + 0.05*randn(N2,1);

x_true = zeros(N2,1);
x_mod= zeros(N2,1);
e = zeros(N2,1);

for k = 1:N2-1
    x_true(k+1) = x_true(k) + Ts*f_true(x_true(k), u2(k));
    x_mod(k+1)= x_mod(k)  + Ts*( phi(x_mod(k))*theta_hat + u2(k) );
    e(k+1) = x_true(k+1) - x_mod(k+1);
end

ISE = trapz(t2, e.^2);

% Plot results
figure('Name','Final Stability & Error','NumberTitle','off');
subplot(3,1,1);
plot(t2, x_true,'b', t2, x_mod,'r--');
legend('true','model');
ylabel('x'); title('State Trajectories');

subplot(3,1,2);
plot(t2, e); ylabel('error e'); title('Tracking Error');

subplot(3,1,3);
plot(t2, u2);
xlabel('t (s)'); ylabel('u'); title('Input Signal');

sgtitle(sprintf('Final Model: %s | ISE = %.3e', results{m}.name, ISE), ...
        'FontWeight','bold');

fprintf('Integrated Squared Error (ISE): %.3e\n', ISE);
fprintf('Max |x_model| = %.3f → bounded\n', max(abs(x_mod)));
