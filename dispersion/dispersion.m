%================================Constants==============================
V       =  2*10e-6;
Gamma   =  1.0;
D       =  10e-9;
k       =  0.4;
mu_eq   =  1.0;
G       = -1*mu_eq*(1-k)/D
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
fname = "dispersion.dat";
fp    = fopen(fname, 'w');
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
omega_init  = 10e5;
omega_final = 100000.0;
omega_step  = 1000000.0;
indx        = nsteps/omega_step;
ampl_fac    = zeros(indx);
omega       = zeros(indx);
omega(0)    = omega_init;
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
for i = 1:indx
  omega(i)      =  omega(i-1) + omega_step;
  first_term    =  V/(2.0*D);
  first_term_sq =  first_term*first_term;
  k_omega       =  first_term + sqrt(first_term_sq + (omega(i)*omega(i)));
  omega_tilda   =  k_omega - (V*(1.0-k)/D);
  b             =  Gamma*omega(i)*omega(i)/(mu_eq*(1-k));
  ampl_fac(i)   =  omega_tilda*V*(-(b/G) + (1.0/omega_tilda)*(k_omega - (V/D)));
endfor
disp(omega);
plot(om, ampl_fac)
fclose(fp);
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
