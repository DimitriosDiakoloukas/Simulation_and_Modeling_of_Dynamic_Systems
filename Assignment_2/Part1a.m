clc;
clear;
close all;

% === System parameters (Παράμετροι Συστήματος) ===
m = 1.315;   % mass (μάζα) [kg]
b = 0.225;   % damping coefficient (συντελεστής απόσβεσης) [N·s/m]
k = 0.725;   % stiffness (δυσκαμψία) [N/m]

% === Filter poles (Πόλοι φίλτρων) ===
poly_coeffs = [1 -1 1];    % s^2 - s + 1
poles = roots(poly_coeffs);

p1 = poles(1);   % 0.5 + 0.8660i
p2 = poles(2);   % 0.5 - 0.8660i


% === Simulation settings (Ρυθμίσεις προσομοίωσης) ===
tspan = 0:0.01:20;   % simulation time (χρονικό διάστημα προσομοίωσης)
x0 = [  % initial conditions (αρχικές συνθήκες)
    0;    % x1(0): initial displacement (αρχική μετατόπιση)
    0;    % x2(0): initial velocity (αρχική ταχύτητα)
    0;    % j1(0): initial filter output for x_dot (αρχικό φιλτραρισμένο ẋ)
    0;    % j2(0): initial filter output for x (αρχικό φιλτραρισμένο x)
    0;    % j3(0): initial filter output for u (αρχικό φιλτραρισμένο u)
    0;    % dj1(0): derivative of j1 (παράγωγος του αρχικού φιλτραρισμένου ẋ)
    0;    % dj2(0): derivative of j2 (παράγωγος του αρχικού φιλτραρισμένου x)
    0;    % dj3(0): derivative of j3 (παράγωγος του αρχικού φιλτραρισμένου u)
    0.01; % theta1(0): initial estimate for theta1 (αρχική εκτίμηση για θ1)
    0.01; % theta2(0): initial estimate for theta2 (αρχική εκτίμηση για θ2)
    0.01  % theta3(0): initial estimate for theta3 (αρχική εκτίμηση για θ3)
];

for i = 1:2
    % === Learning rates (Ρυθμοί μάθησης για κάθε παράμετρο) ===
    if i == 1
        % Constant input case
        gamma1 = 0.114;   % learning rate for theta1 (σχετίζεται με απόσβεση b)
        gamma2 = 0.058;   % learning rate for theta2 (σχετίζεται με δυσκαμψία k)
        gamma3 = 0.062;   % learning rate for theta3 (σχετίζεται με μάζα m)
    else
        % Sinusoidal input case
        gamma1 = 0.05;
        gamma2 = 0.05;
        gamma3 = 0.0495;
    end

    if i == 1
        c = 0; d = 0; f = 2.5;  % Constant input: U(t) = 2.5
        curr_case = "Constant Input: U(t) = 2.5";
    else
        c = 2.5; d = 1; f = 0;  % Sinusoidal input: U(t) = 2.5·sin(t)
        curr_case = "Sinusoidal Input: U(t) = 2.5·sin(t)";
    end

    % Parameters struct (όλα τα δεδομένα σε μία δομή για το ODE solver)
    params = struct('m',m,'b',b,'k',k,...
                    'gamma1',gamma1,'gamma2',gamma2,'gamma3',gamma3,...
                    'c',c,'d',d,'f',f,'p1',p1,'p2',p2);

    % === Solve the system (Επίλυση του συστήματος) ===
    [t, x] = ode45(@(t, x) grad_estimator(t, x, params), tspan, x0);

    % === Extract states (Ανάλυση των μεταβλητών) ===
    x_true = x(:,1);    % true displacement (πραγματική μετατόπιση)
    j1 = x(:,3); j2 = x(:,4); j3 = x(:,5);  % filtered signals (φιλτραρισμένα σήματα)
    th1 = x(:,9); th2 = x(:,10); th3 = x(:,11); % parameter estimates (εκτιμήσεις παραμέτρων)

    % === Reconstructed estimation ===
    x_hat = th1 .* j1 + th2 .* j2 + th3 .* j3;  % estimated displacement (εκτιμώμενη μετατόπιση)
    e_x = x_true - x_hat;                      % estimation error (σφάλμα εκτίμησης)

    % === Parameter recovery ===
    m_hat = 1 ./ th3;                 % estimated mass (εκτιμώμενη μάζα)
    b_hat = (th1 + (p1+p2)) .* m_hat;  % estimated damping (εκτιμώμενη απόσβεση)
    k_hat = (th2 + (p1*p2)) .* m_hat;  % estimated stiffness (εκτιμώμενη δυσκαμψία)

    % === Plots (Γραφήματα) ===
    figure;
    plot(t, x_true, 'b', 'LineWidth', 1.5); hold on;
    plot(t, x_hat, 'r--', 'LineWidth', 1.5);
    legend('x(t)', 'x̂(t)', 'Location', 'best'); grid on;
    xlabel('time in sec'); ylabel('Displacement');
    title(['x(t) and x̂(t) – ', curr_case]);

    figure;
    plot(t, e_x, 'k', 'LineWidth', 1.5);
    xlabel('time in sec'); ylabel('e_x(t) = x(t) - x̂(t)');
    title(['Estimation Error – ', curr_case]); grid on;

    figure;
    subplot(3,1,1);
    sgtitle(['Estimated Parameters – ', curr_case], 'FontWeight', 'bold');
    plot(t, m_hat, 'b', 'LineWidth', 1.5); hold on;
    yline(m, '--r', 'LineWidth', 1.5);
    ylabel('m̂(t)');
    xlabel('time in sec'); grid on;

    subplot(3,1,2);
    plot(t, b_hat, 'b', 'LineWidth', 1.5); hold on;
    yline(b, '--r', 'LineWidth', 1.5);
    ylabel('b̂(t)');
    xlabel('time in sec'); grid on;

    subplot(3,1,3);
    plot(t, k_hat, 'b', 'LineWidth', 1.5); hold on;
    yline(k, '--r', 'LineWidth', 1.5);
    ylabel('k̂(t)');
    xlabel('time in sec'); grid on;

    % === Print final results (Εκτύπωση τελικών εκτιμήσεων) ===
    fprintf('\n--- Input Case: %s ---\n', curr_case);
    fprintf('Final Estimates at t = %.2fs\n', t(end));
    fprintf('Mass:      true = %.4f,  estimated = %.4f\n', m, m_hat(end));
    fprintf('Damping:   true = %.4f,  estimated = %.4f\n', b, b_hat(end));
    fprintf('Stiffness: true = %.4f,  estimated = %.4f\n', k, k_hat(end));
end

% === Gradient estimator dynamics (Δυναμικά Εκτιμητή) ===
function dx = grad_estimator(t, x, P)
    % Extract states
    x1 = x(1); x2 = x(2);
    j1 = x(3); j2 = x(4); j3 = x(5);
    dj1 = x(6); dj2 = x(7); dj3 = x(8);
    th1 = x(9); th2 = x(10); th3 = x(11);

    % Input signal (είσοδος u(t))
    u = P.c * sin(P.d * t) + P.f;

    % True system parameters (αληθινοί συντελεστές συστήματος)
    theta_true = [-P.b/P.m, -P.k/P.m, 1/P.m];
    
    % True system dynamics (δυναμική του συστήματος)
    x2_dot = theta_true(1) * x2 + theta_true(2) * x1 + theta_true(3) * u;

    % Estimated output (εκτίμηση εξόδου)
    x_hat = th1*j1 + th2*j2 + th3*j3;

    % Estimation error (σφάλμα εκτίμησης)
    e = x1 - x_hat;

    % Dynamics (Διαφορικές εξισώσεις)
    dx = zeros(11,1);
    dx(1) = x2;
    dx(2) = x2_dot;
    dx(3) = dj1;
    dx(4) = dj2;
    dx(5) = dj3;
    dx(6) = - (P.p1 + P.p2) * dj1 - (P.p1 * P.p2) * j1 - x2;
    dx(7) = - (P.p1 + P.p2) * dj2 - (P.p1 * P.p2) * j2 - x1;
    dx(8) = - (P.p1 + P.p2) * dj3 - (P.p1 * P.p2) * j3 + u;
    dx(9) = P.gamma1 * e * j1;
    dx(10)= P.gamma2 * e * j2;
    dx(11)= P.gamma3 * e * j3;
end
