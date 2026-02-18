%% ============================================================
%  Kalman Filter (LQG) for Nonlinear Spacecraft Model
%  Based on: "PID Design for Nonlinear Spacecraft Model with
%             Reaction Wheels" (Osman et al.)
%
%  Combines the LQR controller with a Kalman Filter observer
%  to create a full LQG (Linear Quadratic Gaussian) regulator.
%  Noisy sensor measurements of theta are used to reconstruct
%  all three system states.
%% ============================================================
clc; clear; close all;

%% ── 1. SYSTEM PARAMETERS ─────────────────────────────────────
I_sc  = 1;    % Spacecraft moment of inertia
A_sc  = 200;  % Cross-inertia term
J     = 3;    % Reaction wheel inertia
K     = 0.01; % Motor torque constant
L_arm = 0.5;  % Armature inductance
R_arm = 1.0;  % Armature resistance
I_eq  = I_sc + A_sc + 2*J;   % 207

%% ── 2. STATE-SPACE MATRICES (same as LQR script) ─────────────
A_mat = [0,         1,           0;
         0,         0,     K/I_eq;
         0,   -K/L_arm,  -R_arm/L_arm];

B_mat = [0; 0; 1/L_arm];

% Full-state output (all states measured — Kalman will clean up noise)
C_full = eye(3);     % Used for Kalman design
C_out  = [1, 0, 0];  % Only theta is the "true" plant output

D_mat = zeros(3,1);  % No direct feed-through

%% ── 3. LQR GAINS (mirror from lqr_spacecraft.m) ──────────────
Q_lqr   = diag([500, 50, 1]);
R_lqr   = 0.01;
[K_lqr, ~, ~] = lqr(A_mat, B_mat, Q_lqr, R_lqr);

% Pre-compensator
A_cl   = A_mat - B_mat * K_lqr;
N_bar  = -1 / (C_out * (A_cl \ B_mat));

%% ── 4. KALMAN FILTER DESIGN ───────────────────────────────────
% Process noise covariance Q_kf: uncertainty in the dynamics
% Measurement noise covariance R_kf: sensor noise power
% Higher R_kf  → filter trusts model more
% Higher Q_kf  → filter trusts measurements more

Q_kf = diag([1e-6, 1e-5, 1e-4]);   % process noise
R_kf = diag([0.01, 100, 100]);      % measurement noise
%   theta measured (low noise), omega & current not measured directly
%   (high noise forces filter to estimate them from the model)

% Kalman gain via discrete DLQE / continuous Kalman
sys_kf = ss(A_mat, eye(3), C_full, zeros(3,3));
[~, L_kf, ~, ~] = kalman(sys_kf, Q_kf, R_kf);

fprintf('=== Kalman Observer Gain L (3×3) ===\n');
disp(L_kf);

% Observer A matrix: A_obs = A - L*C
A_obs = A_mat - L_kf * C_full;
fprintf('Observer eigenvalues:\n');
disp(eig(A_obs));

%% ── 5. COMBINED LQG SIMULATION ───────────────────────────────
dt   = 0.005;           % time step (s)
T    = 10;              % simulation horizon (s)
t    = (0:dt:T)';
N    = length(t);
ref  = 1;               % 1-rad step reference

% State and estimate initialisation
x_true = zeros(3,1);    % true state
x_hat  = zeros(3,1);    % Kalman estimate

% Storage
x_true_hist = zeros(N,3);
x_hat_hist  = zeros(N,3);
y_noisy_hist= zeros(N,3);
u_hist      = zeros(N,1);
theta_hist  = zeros(N,1);

% Noise standard deviations
sigma_theta = 0.05;   % rad  (sensor noise on angle)
sigma_omega = 0.2;    % rad/s – not directly measured
sigma_i     = 0.5;    % A    – not directly measured

rng(42);  % reproducible

for k = 1:N
    % ── A. Compute control from ESTIMATED states ──────────────
    u = N_bar * ref - K_lqr * x_hat;
    u = max(min(u, 50), -50);    % saturation ±50 V

    % ── B. True plant (Euler integration + process noise) ─────
    w_proc = [1e-3*randn; 1e-3*randn; 1e-3*randn];
    x_dot  = A_mat * x_true + B_mat * u + w_proc;
    x_true = x_true + dt * x_dot;

    % ── C. Noisy measurement (only theta directly sensed) ─────
    v_meas   = [sigma_theta * randn; sigma_omega * randn; sigma_i * randn];
    y_noisy  = C_full * x_true + v_meas;  % all 3 "channels" (2&3 very noisy)

    % ── D. Kalman Predict ─────────────────────────────────────
    x_hat_pred = x_hat + dt * (A_mat * x_hat + B_mat * u);

    % ── E. Kalman Update ──────────────────────────────────────
    innovation = y_noisy - C_full * x_hat_pred;
    x_hat      = x_hat_pred + L_kf * (dt * innovation);

    % ── Store ─────────────────────────────────────────────────
    x_true_hist(k,:) = x_true';
    x_hat_hist(k,:)  = x_hat';
    y_noisy_hist(k,:)= y_noisy';
    u_hist(k)        = u;
    theta_hist(k)    = x_true(1);
end

%% ── 6. PERFORMANCE METRICS ────────────────────────────────────
theta_true = x_true_hist(:,1);
theta_est  = x_hat_hist(:,1);

% Settling time: last time |theta - ref| > 5% of ref
band = 0.05 * ref;
settled_idx = find(abs(theta_true - ref) > band, 1, 'last');
Ts_sim = t(settled_idx);

% Overshoot
OS_sim = max(0, (max(theta_true) - ref) / ref * 100);

fprintf('\n=== LQG Simulation Performance ===\n');
fprintf('  Settling Time : %.3f s\n', Ts_sim);
fprintf('  Overshoot     : %.3f %%\n', OS_sim);
fprintf('  Final theta   : %.5f rad (ref = 1 rad)\n', theta_true(end));
fprintf('  Est. RMSE     : %.5f rad\n', rms(theta_true - theta_est));

%% ── 7. PLOT: Angle Tracking (true vs estimated vs reference) ──
fig1 = figure('Name','Kalman State Estimation','Position',[50 50 1000 620]);

subplot(2,2,1);
plot(t, theta_true,   'b-',  'LineWidth', 2); hold on;
plot(t, theta_est,    'r--', 'LineWidth', 1.5);
plot(t, y_noisy_hist(:,1), 'Color',[0.7 0.7 0.7], 'LineWidth', 0.8);
yline(1.0,'k--','Reference','LineWidth',1.2);
xlabel('Time (s)'); ylabel('\theta (rad)');
title('Angle \theta: True vs Estimated vs Noisy Measurement');
legend('True','Kalman Estimate','Noisy Sensor','Location','southeast');
grid on; xlim([0 T]);

subplot(2,2,2);
omega_true = x_true_hist(:,2);
omega_est  = x_hat_hist(:,2);
plot(t, omega_true, 'b-',  'LineWidth', 2); hold on;
plot(t, omega_est,  'r--', 'LineWidth', 1.5);
plot(t, y_noisy_hist(:,2), 'Color',[0.7 0.7 0.7], 'LineWidth', 0.8);
xlabel('Time (s)'); ylabel('\omega (rad/s)');
title('Angular Rate \omega: True vs Estimated');
legend('True','Kalman Estimate','Noisy Measurement','Location','northeast');
grid on; xlim([0 T]);

subplot(2,2,3);
i_true = x_true_hist(:,3);
i_est  = x_hat_hist(:,3);
plot(t, i_true, 'b-',  'LineWidth', 2); hold on;
plot(t, i_est,  'r--', 'LineWidth', 1.5);
plot(t, y_noisy_hist(:,3), 'Color',[0.7 0.7 0.7], 'LineWidth', 0.8);
xlabel('Time (s)'); ylabel('i_{arm} (A)');
title('Armature Current: True vs Estimated');
legend('True','Kalman Estimate','Noisy Measurement','Location','northeast');
grid on; xlim([0 T]);

subplot(2,2,4);
est_err = theta_true - theta_est;
plot(t, est_err, 'm-', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Error (rad)');
title('Kalman Estimation Error (\theta_{true} - \hat{\theta})');
yline(0,'k--'); grid on; xlim([0 T]);

sgtitle('Kalman Filter State Estimation — LQG Controller', ...
        'FontSize', 14, 'FontWeight', 'bold');
saveas(fig1, '/home/claude/kalman_state_estimation.png');
fprintf('Saved: kalman_state_estimation.png\n');

%% ── 8. PLOT: LQG vs LQR (clean) comparison ───────────────────
% Run clean LQR (no noise) for comparison
sys_lqr_cl = ss(A_cl, B_mat*N_bar, C_out, 0);
[y_lqr, t_lqr] = step(sys_lqr_cl, t);

fig2 = figure('Name','LQG vs LQR Comparison','Position',[100 200 900 480]);
plot(t_lqr, y_lqr,     'b-',  'LineWidth', 2.5); hold on;
plot(t,     theta_true, 'r--', 'LineWidth', 2.0);
plot(t,     theta_est,  'g:',  'LineWidth', 2.0);
yline(1.0,'k--','Reference','LineWidth',1.2);
xlabel('Time (s)', 'FontSize', 13);
ylabel('Angular Position \theta (rad)', 'FontSize', 13);
title('LQR (ideal) vs LQG (noisy sensor + Kalman filter)', ...
      'FontSize', 14, 'FontWeight', 'bold');
legend('LQR – no noise','LQG – true state (noisy plant)','LQG – Kalman estimate', ...
       'Location','southeast');
grid on; xlim([0 T]); ylim([-0.05 1.3]);
saveas(fig2, '/home/claude/lqg_vs_lqr.png');
fprintf('Saved: lqg_vs_lqr.png\n');

%% ── 9. PLOT: Control Input ────────────────────────────────────
fig3 = figure('Name','LQG Control Input','Position',[100 200 900 400]);
plot(t, u_hist, 'r-', 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 13);
ylabel('Control Voltage V (V)', 'FontSize', 13);
title('LQG Control Effort (with Kalman state estimates)', ...
      'FontSize', 14, 'FontWeight', 'bold');
grid on; xlim([0 T]);
saveas(fig3, '/home/claude/lqg_control_effort.png');
fprintf('Saved: lqg_control_effort.png\n');

%% ── 10. PLOT: Innovation / Residuals ─────────────────────────
residuals = y_noisy_hist(:,1) - theta_est;
fig4 = figure('Name','Kalman Residuals','Position',[100 200 900 350]);
plot(t, residuals, 'Color',[0.2 0.6 0.2], 'LineWidth', 1.2);
xlabel('Time (s)', 'FontSize', 13);
ylabel('Innovation (rad)', 'FontSize', 13);
title('Kalman Filter Innovation (Measurement Residuals)', ...
      'FontSize', 14, 'FontWeight', 'bold');
yline(0,'k--'); grid on; xlim([0 T]);
saveas(fig4, '/home/claude/kalman_residuals.png');
fprintf('Saved: kalman_residuals.png\n');

fprintf('\nKalman + LQG script complete.\n');
