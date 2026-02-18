%% ============================================================
%  LQR CONTROLLER DESIGN — Nonlinear Spacecraft with Reaction Wheels
%  Reference: "PID Design for Nonlinear Spacecraft Model with
%              Reaction Wheels" — Osman et al., Elektrika Journal (UTM)
%
%  Plant Transfer Function (from paper):
%    G(s) = -0.01 / [207 × (0.05s² + 0.06s + 0.1001)]
%
%  Performance Requirements:  Ts ≤ 4 s   |   OS ≤ 10%
%% ============================================================
clc; clear; close all;

fprintf('==============================================\n');
fprintf('   LQR SPACECRAFT ATTITUDE CONTROLLER\n');
fprintf('==============================================\n\n');

%% ── 1. SYSTEM PARAMETERS (from paper) ────────────────────────
I_sc  = 1;      % Spacecraft moment of inertia (kg.m²)
A_sc  = 200;    % Reaction wheel assembly inertia
J_rw  = 3;      % Reaction wheel inertia
I_eq  = I_sc + A_sc + 2*J_rw;  % Effective coupling denominator = 207

K     = 0.01;   % Motor torque constant (N.m/A)
L_arm = 0.5;    % Armature inductance  (H)     [JL=0.05 → L=0.5]
R_arm = 1.0;    % Armature resistance  (Ω)     [bR=0.1  → R=1.0]
K_e   = K;      % Back-EMF constant (equal to torque constant)

fprintf('System Parameters:\n');
fprintf('  I_eq  = %d   (I + A + 2J)\n', I_eq);
fprintf('  K     = %.3f N.m/A\n', K);
fprintf('  L_arm = %.1f H,  R_arm = %.1f Ω\n\n', L_arm, R_arm);

%% ── 2. STATE-SPACE MODEL ──────────────────────────────────────
% States:  x = [θ (rad),  ω (rad/s),  i_a (A)]
%
% Equations:
%   dθ/dt    =  ω
%   dω/dt    =  (K / I_eq) × i_a          [linearised spacecraft + wheel]
%   di_a/dt  = -(K_e/L)*ω  - (R/L)*i_a + (1/L)*V
%
A_mat = [0,       1,          0;
         0,       0,      K/I_eq;
         0,  -K_e/L_arm,  -R_arm/L_arm];

B_mat = [0; 0; 1/L_arm];

C_theta = [1, 0, 0];   % Observe attitude angle θ
C_full  = eye(3);       % Full-state observation (for analysis)

D_mat = zeros(1,1);

sys_ol = ss(A_mat, B_mat, C_theta, 0);
sys_ol.StateName  = {'theta','omega','i_a'};
sys_ol.InputName  = {'V (input voltage)'};
sys_ol.OutputName = {'theta (rad)'};

fprintf('Open-Loop State-Space:\n');
disp(sys_ol);

%% ── 3. STABILITY & CONTROLLABILITY CHECK ─────────────────────
fprintf('Controllability rank : %d / %d  → ', ...
    rank(ctrb(A_mat, B_mat)), size(A_mat,1));
if rank(ctrb(A_mat, B_mat)) == size(A_mat,1)
    fprintf('CONTROLLABLE ✓\n');
else; fprintf('NOT CONTROLLABLE ✗\n'); end

fprintf('Observability rank   : %d / %d  → ', ...
    rank(obsv(A_mat, C_theta)), size(A_mat,1));
if rank(obsv(A_mat, C_theta)) == size(A_mat,1)
    fprintf('OBSERVABLE ✓\n\n');
else; fprintf('NOT OBSERVABLE ✗\n\n'); end

fprintf('Open-Loop Eigenvalues:\n');
disp(eig(A_mat));

%% ── 4. LQR WEIGHT SELECTION (Bryson's Rule) ──────────────────
% Max acceptable deviations:
%   θ   ≤ 0.1 rad   → Q(1,1) = 1/0.1² = 100  (most critical)
%   ω   ≤ 0.5 rad/s → Q(2,2) = 1/0.5² = 4    → rounded to 10 for speed
%   i_a ≤ 2.0 A     → Q(3,3) = 1/2.0² = 0.25 → 1
% Max voltage ≤ 50 V → R = 1/50² ≈ 4e-4       → 1e-6 for aggressiveness
Q_lqr = diag([100, 10, 1]);
R_lqr = 1e-6;

[K_lqr, S_lqr, P_cl] = lqr(A_mat, B_mat, Q_lqr, R_lqr);

fprintf('LQR Gain Matrix K = [%.3f,  %.3f,  %.3f]\n', K_lqr);
fprintf('Closed-Loop Eigenvalues (LQR):\n');
disp(P_cl);

%% ── 5. REFERENCE FEEDFORWARD GAIN ────────────────────────────
% Nbar ensures zero steady-state error for θ step reference
A_cl  = A_mat - B_mat * K_lqr;
Nbar  = -1 / (C_theta * (A_cl \ B_mat));
fprintf('Feedforward gain Nbar = %.5f\n\n', Nbar);

%% ── 6. SIMULATION (10 s, unit step reference) ────────────────
t   = 0:0.005:10;
ref = ones(size(t));   % 1 rad step

sys_cl = ss(A_cl, B_mat*Nbar, C_theta, 0);
[y_theta, t_out] = lsim(sys_cl, ref, t, [0;0;0]);

% Full-state for control effort reconstruction
sys_cl_full = ss(A_cl, B_mat*Nbar, eye(3), zeros(3,1));
[x_out, ~]  = lsim(sys_cl_full, ref, t, [0;0;0]);
u_lqr       = Nbar * ref' - (K_lqr * x_out')';

%% ── 7. PERFORMANCE METRICS ────────────────────────────────────
info = stepinfo(y_theta, t_out, 1, 'SettlingTimeThreshold', 0.02);
ess  = abs(1 - y_theta(end));

fprintf('=== LQR Performance Metrics ===\n');
fprintf('  Rise Time     : %.4f s\n', info.RiseTime);
fprintf('  Settling Time : %.4f s   (target ≤ 4 s)\n', info.SettlingTime);
fprintf('  Overshoot     : %.2f%%   (target ≤ 10%%)\n', info.Overshoot);
fprintf('  Steady-State Error : %.6f rad\n\n', ess);

spec_met = (info.Overshoot <= 10) && (info.SettlingTime <= 4);
if spec_met
    fprintf('✓ ALL DESIGN SPECIFICATIONS MET\n\n');
else
    fprintf('✗ Some specifications not met — adjust Q or R\n\n');
end

%% ── 8. PLOT: Step Response ────────────────────────────────────
fig1 = figure('Name','LQR Step Response','NumberTitle','off', ...
              'Position',[50 80 1000 420], 'Color','w');

plot(t_out, y_theta, 'b-', 'LineWidth', 2.5); hold on;
plot(t_out, ref,     'k--','LineWidth', 1.2);
% Spec bands
yline(1.10, 'r:', 'LineWidth', 1.5, 'Label', '+10% OS limit');
patch([0 max(t_out) max(t_out) 0], [0.98 0.98 1.02 1.02], ...
      [0.2 0.8 0.2], 'FaceAlpha', 0.1, 'EdgeColor','none');
xline(4, 'm:', 'LineWidth', 1.5, 'Label', 'Ts target=4s');
% Annotation box
text(6.5, 0.45, sprintf('OS = %.2f%%\nTs = %.2f s\ne_{SS} = %.4f rad', ...
    info.Overshoot, info.SettlingTime, ess), ...
    'FontSize', 11, 'BackgroundColor', [0.95 0.95 1], ...
    'EdgeColor', 'b', 'Margin', 5);

xlabel('Time (s)', 'FontSize', 13);
ylabel('\theta (rad)', 'FontSize', 13);
title(sprintf('LQR Closed-Loop Step Response\n(Q = diag[%g, %g], R = %g×10^{-6})', ...
      Q_lqr(1,1), Q_lqr(2,2), R_lqr*1e6), 'FontSize', 13, 'FontWeight','bold');
legend('LQR Response', 'Reference (1 rad)', 'Location','SouthEast');
grid on; xlim([0 max(t_out)]); ylim([-0.05 1.4]);
saveas(fig1, 'lqr_step_response.png');
fprintf('Saved → lqr_step_response.png\n');

%% ── 9. PLOT: Control Effort ───────────────────────────────────
fig2 = figure('Name','LQR Control Effort','NumberTitle','off', ...
              'Position',[50 550 1000 380], 'Color','w');
plot(t_out, u_lqr, 'r-', 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 13);
ylabel('Control Input V (volts)', 'FontSize', 13);
title('LQR Control Effort — Reaction Wheel Drive Voltage', ...
      'FontSize', 13, 'FontWeight','bold');
grid on; xlim([0 max(t_out)]);
saveas(fig2, 'lqr_control_effort.png');
fprintf('Saved → lqr_control_effort.png\n');

%% ── 10. PLOT: Pole Migration ──────────────────────────────────
ol_eig = eig(A_mat);
cl_eig = eig(A_cl);

fig3 = figure('Name','Pole Migration','NumberTitle','off', ...
              'Position',[1060 80 600 480], 'Color','w');
plot(real(ol_eig), imag(ol_eig), 'rx', 'MarkerSize',14,'LineWidth',3); hold on;
plot(real(cl_eig), imag(cl_eig), 'b^', 'MarkerSize',12,'LineWidth',2,'MarkerFaceColor','b');
xline(0,'k--'); yline(0,'k--');
% Label poles
for k = 1:length(ol_eig)
    text(real(ol_eig(k))+0.02, imag(ol_eig(k))+0.04, ...
         sprintf('%.3f%+.3fj', real(ol_eig(k)), imag(ol_eig(k))), ...
         'Color','r','FontSize',9);
end
for k = 1:length(cl_eig)
    text(real(cl_eig(k))+0.02, imag(cl_eig(k))+0.04, ...
         sprintf('%.3f%+.3fj', real(cl_eig(k)), imag(cl_eig(k))), ...
         'Color','b','FontSize',9);
end
xlabel('Real Axis (s^{-1})','FontSize',12);
ylabel('Imaginary Axis (s^{-1})','FontSize',12);
title('Pole Migration: Open-Loop → LQR Closed-Loop', ...
      'FontSize', 13, 'FontWeight','bold');
legend('Open-Loop Poles','LQR CL Poles','Location','Best');
grid on;
saveas(fig3, 'lqr_poles.png');
fprintf('Saved → lqr_poles.png\n\n');

fprintf('LQR design script complete.\n');
fprintf('Run  kalman_lqg_spacecraft.m  next for the Kalman filter + LQG.\n');
