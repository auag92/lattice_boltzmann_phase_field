
/*Phasefield module*/
double *phi_new, *phi_old;
double *mu_new, *mu_old;
double *lap_phi, *lap_mu;
double *conc;
double *dphi_now,  *dphi_next;
//----------------------------------------------------------
/*Fluid module*/
double *P; //Pressure
double *fn;
double *u_old, *u_now; // velocity in x direction
double *v_old, *v_now; // velocity in y direction
double *v_str, *u_str;
double *a_x, *a_y;
double *rhs_fn; // rhs function that is passed on to the multigrid module
double *Hx, *Hy; // saving the Hx and Hy stuff from the fluid based calculations
