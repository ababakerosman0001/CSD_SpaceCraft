G = tf(-0.01/207, [0.05 0.06 0.1001]);
K = -800
% Desired pole
OS = 0.10; Ts = 4;
zeta = -log(OS)/sqrt(pi^2 + log(OS)^2);
sigma = 4/Ts;
wd = sigma * sqrt(1-zeta^2)/zeta;
sd = -sigma + 1i*wd;

%PI Controller

Cpd = tf([1 1],[1])
Cpi = tf([1 0.9],[1 0])
figure(1)
rlocus(Cpd*G*Cpi)
figure(2)
step(feedback(K*Cpi*G*Cpi,1))