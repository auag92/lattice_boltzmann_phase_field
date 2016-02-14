#include "constants.h"
#include "variables.h"
//*****************************************************************************//
void    phi_allocate_memory();
void    phi_initialize();
void    phi_boundary(double *c);
void    phi_update();
void    laplacian(double *f, double *lap);
void    phi_solver();
void    concentration();
void    isotropic_solverloop();
void    anisotropic_solverloop();
void    grad_phi(int i, double *d_phi);
double  dqdx( double phi_x, double phi_y);
double  div_phi(int i);
void    fnupdate();
void    phi_free_memory();
double  Mobility(int indx);
//*****************************************************************************//
double inv_deltax2 = (1.0/deltax2);
//*****************************************************************************//
void phi_allocate_memory(){
  phi_old   = (double *)malloc(Mx*My*sizeof(double));
  phi_new   = (double *)malloc(Mx*My*sizeof(double));
  mu_old    = (double *)malloc(Mx*My*sizeof(double));
  mu_new    = (double *)malloc(Mx*My*sizeof(double));
  conc      = (double *)malloc(Mx*My*sizeof(double));
  lap_phi   = (double *)malloc(Mx*My*sizeof(double));
  lap_mu    = (double *)malloc(Mx*My*sizeof(double));
  dphi_now  = (double *)malloc(Mx*4*sizeof(double));
  dphi_next = (double *)malloc(Mx*4*sizeof(double));
}
//*****************************************************************************//
void laplacian(double *f, double *lap) {
  long i, j, z;
  for (i = 1; i < Mx -1; i++)
  {
    for (j = 1; j < My -1; j++)
    {
      z = i*My + j;
      lap[z] = (f[z-1] + f[z+1] -2.0*f[z])*inv_deltax2 + (f[z+My] + f[z-My] -2.0*f[z])*inv_deltax2 ;
    }
  }
 }
//*****************************************************************************//
void phi_solver(){
  laplacian(phi_old,     lap_phi);
  laplacian(mu_old,     lap_mu);
  concentration();
  #ifdef ISO
    isotropic_solverloop();
  #endif
  #ifdef ANISO
    anisotropic_solverloop();
  #endif
  phi_boundary(phi_new);
  phi_boundary(mu_new);
  phi_update();
}
//*****************************************************************************//
void phi_update() {
  long i, j, indx;
  for (i = 0; i < Mx; i++) {
    for (j = 0; j < My; j++){
      indx          = i*My + j;
      phi_old[indx] = phi_new[indx];
      mu_old[indx]  = mu_new[indx];
    }
  }
}
//*****************************************************************************//
void concentration(){
  double p, mu,h;
  int i, j, z;
  for ( i = 1; i < Mx-1; i++)
  {
    for ( j = 1; j < My-1; j++){
      z       =  i*My + j;
      p       =  phi_new[z];
      mu      =  mu_new[z];
      h       =  p*p*(3-2*p);
      conc[z] = mu*(1-h) + (h)*K*mu;
    }
  }
}
//*****************************************************************************//
void phi_initialize() {
  long i,j,z;
  double r;
#ifdef Centre
  for ( i = 0; i < Mx; i++)
  {
    for ( j=0; j < My; j++)
    {
      r= (i-Mx*0.5)*(i-Mx*0.5) + (j-My*0.5)*(j-My*0.5);
      z= i*My + j;
      if(r < RADIUS2){
	      phi_old[z] = 1.0;
      }
      else{
	      phi_old[z] = 0.0;
      }
      mu_old[z] = Mu - deltaMu;
    }
  }
#endif
#ifdef Corner
  for ( i = 0; i < Mx; i++)
  {
    for ( j=0; j < My; j++)
    {
      r= (i)*(i) + (j)*(j);
      z= i*My + j;
      if(r < RADIUS2){
	      phi_old[z] = 1.0;
      }
      else{
	      phi_old[z] = 0.0;
      }
      mu_old[z] = Mu - deltaMu;
    }
  }
#endif
#ifdef Nothing
  for ( i = 0; i < MESHX; i++)
  {
    for ( j=0; j < MESHX; j++)
    {
      z= i*MESHX + j;
      phi_old[z] = 0.0;
      mu_old[z] = 0.0;
    }
  }
#endif
#ifdef PlaneFront
for ( i = 0; i < Mx; i++)
{
  for ( j=0; j < My; j++)
  {
    z= i*My + j;
    if(i < Init_Length){
      phi_old[z] = 1.0;
    }
    else{
      phi_old[z] = 0.0;
    }
    mu_old[z] = Mu - deltaMu;
  }
}
#endif
#ifdef SineFront
  double delta;
  for ( i = 0; i < Mx; i++)
  {
    for ( j=0; j < My; j++)
    {
      z= i*My + j;
      delta = Init_Length + Amp*sin(M_PI/My*j + M_PI/2.0);
      if(i <= delta){
        phi_old[z] = 1.0;
      }
      else{
        phi_old[z] = 0.0;
      }
      if( i == delta){
        mu_old[z] = Mu;
      }else{
        mu_old[z] = Mu - deltaMu;
      }
    }
  }
#endif
}
//*****************************************************************************//
void phi_boundary(double *c) {

  int i ,j ,y ,z;
  //Top - Bottom
  for (j = 0; j <= My -1; j++)
  {
    c[j]            = c[My + j];
    c[MESHX2-My+j]  = c[MESHX2-2*My+j];
  }
  //Left - Right
  for (i = 1; i < Mx -1; i++)
  {
    y     = i*My;
    z     = i*My + My-1;

    c[y]  = c[y+1];
    c[z]  = c[z-1];
  }
}
//*****************************************************************************//
void isotropic_solverloop() {
  int    i,j,z;
  double p,dp_dt,dmu_dt, kai;
  double dc_dx, dc_dy, V_gradC = 0.0;
  double mob;

  for (i = 1; i < Mx-1; i++) {
    for (j = 1; j < My-1; j++){
      z = i*My + j;

      p = phi_old[z];
      dp_dt = (G*E*lap_phi[z] - (G/E)*18.0*(p)*(1.0-p)*(1.0-2.0*p) + (mu_old[z] - Mu)*(K-1)*(mu_old[z])*6*p*(1-p))/(tau*E);
      phi_new[z] = p + deltat*(dp_dt);

      dc_dx = (conc[z+1]-conc[z-1])*0.5*inv_deltax;
      dc_dy = (conc[z+My]-conc[z-My])*0.5*inv_deltax;
      // V_gradC = u[z]*dc_dx + v[z]*dc_dy;

      mob   = Mobility(z);
      // dmu_dt = Mob*lap_mu[z] - V_gradC - (K-1)*mu_old[z]*6*p*(1-p)*dp_dt;
      dmu_dt = mob - V_gradC - (K-1)*mu_old[z]*6*p*(1-p)*dp_dt;
      kai    = 1+(K-1)*p*p*(3-2*p);

      mu_new[z] = mu_old[z]  + deltat*dmu_dt/kai;
    }
  }
}
/******************************************************************************/
// void anisotropic_solverloop(){
//
//   int       i, j, z;
//   double    p,dp_dt,dmu_dt;
//   double    drv_frce, alln_chn;
//   double    Gamma, kai;
//   double    dc_dx, dc_dy, V_gradC = 0.0;
//
//   grad_phi(1, dphi_now);
//
//   for (i = 1; i < (Mx-1); i++) {
//     grad_phi(i+1, dphi_next);
//     for (j=1; j < (My-1); j++){
//
//       z =   i*My + j;
//       p =   phi_old[z];
//
//       Gamma         =     div_phi(j);           // Gamma = 2*G*lap_phi[z];
//       drv_frce      =     (mu_old[z] - Mu)*(K-1)*(mu_old[z])*6*p*(1-p);
//       alln_chn      =     E*Gamma - (G/E)*18.0*(p)*(1.0-p)*(1.0-2.0*p);
//       dp_dt         =     (alln_chn + drv_frce)/(tau*E);
//
//       phi_new[z]    =     p + deltat*dp_dt;
//
//       dc_dx         =     (conc[z+1]-conc[z-1])*0.5*inv_deltax;
//       dc_dy         =     (conc[z+MESHX]-conc[z-MESHX])*0.5*inv_deltax;
//       // V_gradC       =     u[z]*dc_dx + v[z]*dc_dy;
//
//       dmu_dt        =     Mob*lap_mu[z] - V_gradC - (K-1)*mu_old[z]*6*p*(1-p)*dp_dt;
//       // dmu_dt        =     Mob*lap_mu[z] - (K-1)*mu_old[z]*6*p*(1-p)*dp_dt;
//       kai           =     1+(K-1)*p*p*(3-2*p);
//       mu_new[z]     =     mu_old[z]  + deltat*dmu_dt/kai;
//     }
//     fnupdate();
//   }
// }
// /******************************************************************************/
// void grad_phi(int i, double *d_phi){
//
//  	int j,z;
//
//   for(j=1; j<MESHX-1; j++){
//     z = i * MESHX + j;
//     if (  i == MESHX -1 ){
//       d_phi[2*MESHX+j] = (phi_old[z] - phi_old[z-MESHX])*inv_deltax;
//       d_phi[3*MESHX+j] = (phi_old[z+1] - phi_old[z-1] + phi_old[z+1-MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax;
//     }
//     else {
//       d_phi[j] = (phi_old[z] - phi_old[z-1])*inv_deltax;
//       d_phi[MESHX+j] = (phi_old[z+MESHX] - phi_old[z-MESHX] + phi_old[z-1+MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax;
//       d_phi[2*MESHX+j] = (phi_old[z] - phi_old[z-MESHX])*inv_deltax;
//       d_phi[3*MESHX+j] = (phi_old[z+1] - phi_old[z-1] + phi_old[z+1-MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax;
//      }
//   }
//   if (  i != MESHX-1  ) {
//     z = i * MESHX + j;
//     d_phi[j]           = (phi_old[z] - phi_old[z-1])*inv_deltax;
//     d_phi[MESHX+j]     = (phi_old[z+MESHX] - phi_old[z-MESHX] + phi_old[z-1+MESHX] - phi_old[z-1-MESHX])*0.25*inv_deltax;
//   }
// }
// /******************************************************************************/
// double dqdx( double phi_x, double phi_y) {
//
// 	double   a, phi_x2, phi_x4, phi_y2, phi_y4, inv_phi;
// 	int      z;
//   double   ans = 0;
//   double   part1, part2, part3, part4;
//   phi_x2    =   phi_x *phi_x;
//   phi_y2    =   phi_y *phi_y;
//   phi_y4    =   phi_y2 *phi_y2;
//   phi_x4    =   phi_x2 *phi_x2;
//
//   if ((phi_x2> 1e-15) && (phi_y2> 1e-15)){
//     inv_phi    =    1/(phi_x2+phi_y2);
//     part1      =    (1-Dab*(3-4*(phi_x4+phi_y4)*inv_phi*inv_phi));
//     part2      =    2*G*E*part1*part1*phi_x;
//     part3      =    32*G*E*Dab*(phi_x2+phi_y2)*(part1);
//     part4      =    phi_x2*phi_x*inv_phi*inv_phi - phi_x*(phi_x4+phi_y4)*inv_phi*inv_phi*inv_phi;
//
//     ans        =    part2 + part3*part4;
//   }
//
//   return ans;
// }
// /******************************************************************************/
// double div_phi(int i){
//
//   double    ans;
// 	double    x_next,    x_now;
// 	double    y_next,    y_now;
//
//   x_now     = dqdx(dphi_now[i], dphi_now[i+MESHX]);
//   x_next    = dqdx(dphi_now[i+1], dphi_now[i+1+MESHX]);
//   y_now     = dqdx(dphi_now[i+2*MESHX], dphi_now[i+3*MESHX]);
//   y_next    = dqdx(dphi_next[i+2*MESHX], dphi_next[i+3*MESHX]);
// 	ans       = ((x_next - x_now) + ( y_next - y_now))*inv_deltax;
//
//   return ans;
// }
// //*****************************************************************************//
// //*****************************************************************************//
// void fnupdate()
// {
//   int i;
//
//   for( i=0; i < MESHX; i++ ) {
//     dphi_now[i]           =   dphi_next[i];
//     dphi_now[MESHX+i]     =   dphi_next[MESHX+i];
//     dphi_now[2*MESHX+i]   =   dphi_next[2*MESHX+i];
//     dphi_now[3*MESHX+i]   =   dphi_next[3*MESHX+i];
//   }
// }

//*****************************************************************************//
double Mobility(int indx){
  double m1, m2;
  double dmu_x1, dmu_x2;
  double dmu_y1, dmu_y2;
  double ans = 0.0;

  m1 = 1.0 - 0.5*(phi_old[indx]+phi_old[indx-1]);
  m2 = 1.0 - 0.5*(phi_old[indx]+phi_old[indx+1]);

  dmu_x1 = (mu_old[indx]-mu_old[indx-1])*inv_deltax;
  dmu_x2 = (mu_old[indx+1]-mu_old[indx])*inv_deltax;

  ans = ans + (m2*dmu_x2-m1*dmu_x1)*inv_deltax;

  m1 = 1.0 - 0.5*(phi_old[indx]+phi_old[indx-My]);
  m2 = 1.0 - 0.5*(phi_old[indx]+phi_old[indx+My]);

  dmu_y1 = (mu_old[indx]-mu_old[indx-My])*inv_deltax;
  dmu_y2 = (mu_old[indx+My]-mu_old[indx])*inv_deltax;

  ans = ans + (m2*dmu_y2-m1*dmu_y1)*inv_deltax;

  return ans;
}
//*****************************************************************************//
void phi_free_memory() {
  free(phi_old);
  free(phi_new);
  free(mu_old);
  free(mu_new);
  free(conc);
  free(lap_phi);
  free(lap_mu);
  free(dphi_now);
  free(dphi_next);
}
//*****************************************************************************//
