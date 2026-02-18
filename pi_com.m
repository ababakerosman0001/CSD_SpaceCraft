% Plant
G = tf(-0.01/207, [0.05 0.06 0.1001]);

% Pd zero
z = 1.74;
Cpd_unscaled = tf([1 z],1)   % (s+z)/s
Cpid = tf([-2000  -5014  -5912],[1 0])
% Desired pole
OS = 0.10; Ts = 4;
zeta = -log(OS)/sqrt(pi^2 + log(OS)^2);
sigma = 4/Ts;
wd = sigma * sqrt(1-zeta^2)/zeta;
sd = -sigma + 1i*wd;

k = -39583.33

% Root locus
figure(1);
rlocus( Cpid * G); hold on;

sgrid(zeta, []);
%plot(real(sd), imag(sd), 'g*', 'MarkerSize', 12);
figure(2)
step(feedback(Cpid*G,1))
figure(3)
margin(G*k)
% Step response
%figure;
%step(feedback(Cpi * G, 1));
grid on;

disp('Final PID controller:');
disp(zpk(Cpi));
