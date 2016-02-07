syms V c_eq c_str k D Gamma m G omega k_omega

% V       =  2e-6;
% c_eq    =  0.02;
% k       =  0.1389;
% D       =  1e-9;
% Gamma   =  1.03e-6;
% m       =  246;
% G       =  -1*c_eq*(1-k)/D;
% c_str   =  c_eq*(k-1);
% k_omega =  V/(2.0*D) + sqrt(V*V/(4.0*D*D) + omega*omega);
% b       =  -1.0*(Gamma*omega*omega)/m;

V       =  0.1;
c_eq    =  1.0;
k       =  0.4;
D       =  1;
Gamma   =  1;
m       =  246;
G       =  -1*c_eq*(1-k)/(D*V);
c_str   =  c_eq*(k-1);
k_omega =  V/(2.0*D) + sqrt(V*V/(4.0*D*D) + omega*omega);
b       =  (Gamma*omega*omega)/c_str;

eqn     = -D*b/(c_str)*(G/c_eq + k_omega) + G*D*k_omega/c_str*(1-V/(k_omega*D)) == 0;
solx    =  solve(eqn,omega);
sol1x   =  diff(eqn, omega);
sol2x   =  solve(sol1x, omega);
x0      =  vpa(solx)
x1      =  vpa(sol2x)