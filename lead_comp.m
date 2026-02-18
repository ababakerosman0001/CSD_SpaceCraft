% Plant Definition
num = -0.01/207;
den = [0.05, 0.06, 0.1001];
G = tf(num,den);

%% 1. Steady-State Error Requirement (< 5%)
% Kp > 19 to achieve ess < 0.05
K = -39583.33; 
G_u = K * G; 

%% 2. Transient Specifications
OS = 10; 
zeta = -log(OS/100) / sqrt(pi^2 + (log(OS/100))^2);
% Target PM (58.59 deg)
PM_target = atan2d(2*zeta, sqrt(-2*zeta^2 + sqrt(1 + 4*zeta^4)));

%% 3. Iterative Lead Compensator Design
% This loop ensures we hit the required 58.59 degree target exactly
phi_corr = 0; 
tolerance = 0.01;
max_iterations = 15;

for i = 1:max_iterations
    [~, pm_curr, ~, ~] = margin(G_u);
    phi_m = PM_target - pm_curr + phi_corr;
    
    % Calculate beta (alpha)
    beta = (1 - sind(phi_m)) / (1 + sind(phi_m));
    
    % Find new crossover frequency
    mag_shift = -10*log10(1/beta);
    [mag, ~, w] = bode(G_u);
    w_new = interp1(20*log10(squeeze(mag)), w, mag_shift);
    
    % Calculate Zero and Pole
    zc = w_new * sqrt(beta);
    pc = w_new / sqrt(beta);
    D_lead = tf([1/zc 1], [1/pc 1])
    
    % Check achieved PM
    [~, achieved_pm, ~, ~] = margin(D_lead * G_u);
    error = PM_target - achieved_pm;
    if abs(error) < tolerance, break; end
    phi_corr = phi_corr + error; 
end

%% 4. Verification
G_open_loop = D_lead * G_u;
%[gm, pm_final, ~, wcp] = margin(G_open_loop);

figure(1);
margin(G_open_loop); 
grid on;

figure(2);
step(feedback(G_open_loop, 1));
title('Step Response: ess < 5% and PM = 58.59 deg');
grid on;
figure(3)
rlocus(G_open_loop)
fprintf('Achieved Phase Margin: %.2f deg at %.2f rad/s\n', pm_final, wcp);
fprintf('Required Gain K: %.2f\n', K);

