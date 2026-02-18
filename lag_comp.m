%% MATLAB Implementation of Lag Compensator
clear; clc; close all;

% 1. Plant Definition (provided)
num_plant = -0.01/207;
den_plant = [0.05, 0.06, 0.1001];
G = tf(num_plant, den_plant);

% 2. Gain Adjustment (provided)
K = -39583.33; 
G_u = K * G; % Uncompensated system with gain K

% 3. Lag Compensator Parameters (from your calculations)
zc = 0.103;     % Zero location
pc = 0.0042;    % Pole location

% 4. Construct Lag Compensator Transfer Function
% Standard form: Gc(s) = (s/zc + 1) / (s/pc + 1)
num_gc = [1/zc, 1];
den_gc = [1/pc, 1];
Gc = tf(num_gc, den_gc);

% 5. Closed-Loop System
G_comp = series(Gc, G_u); % Total Open-Loop: Gc(s) * K * G(s)

%% 6. Verification and Plots
figure(1);
margin(G_u); % Blue line: Uncompensated
hold on;
margin(G_comp); % Yellow/Orange line: Compensated
grid on;
title('Bode Plot Comparison: Uncompensated vs Lag Compensated');
legend('Uncompensated (with K)', 'Lag Compensated');
figure(2)
step(feedback(G_comp,1))
% Display Performance Data
[gm_u, pm_u, wcp_u, wcg_u] = margin(G_u);
[gm_c, pm_c, wcp_c, wcg_c] = margin(G_comp);

fprintf('--- Design Results ---\n');
fprintf('Uncompensated PM: %.2f deg at %.2f rad/s\n', pm_u, wcg_u);
fprintf('Compensated PM:   %.2f deg at %.2f rad/s\n', pm_c, wcg_c);
fprintf('Target PM was:    52.06 deg\n');