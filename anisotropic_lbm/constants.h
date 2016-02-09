//--------------Parameters for the D2Q9 LBM Solver----------------------------//
#define dx       1
#define dy       dx
#define dt       1.0
#define Mx       200
#define My       Mx
#define M2       Mx*My
#define Q        9   // no. of nodes in the LBM model
#define nu       1.0 // vsicosity
#define Re       u_lid*Mx/nu
#define inv_tau  (1.*dt)/(3.*nu*dt + 0.5)//(deltax*Cs)/(3.*nu+0.5*Cs*deltax)
#define omega    0.9 //dt*inv_tau
#define Rho_init 1.0
//****************************************************************************//
//****************************************************************************//
// Lid driven cavity boundary condition
// #define ldc
#define u_lid   0.5
//****************************************************************************//
//****************************************************************************//
//Pipe flow boundary Condition
#define pipeflow
#define rho_in      1.01
#define rho_out     1.0
#define inv_rho_in  1./rho_in
#define inv_rho_out 1./rho_out
/*****************************************************************************/
/*****************************************************************************/
// Models for calculating Rho
// #define model_1  //compressible flow
#define model_2  //incompressible flow
// #define mdoel_3  //modified incompressible flow

//****************************************************************************//
/*********************Parameters for Phi Solver*******************************/

#define MESHX  Mx
#define MESHX2 (MESHX*MESHX)
#define deltat (0.04)
#define deltax (2.0)
#define inv_deltax (1./deltax)
#define deltax2 (deltax*deltax)

#define K       (0.2)             /*Partition Coefficient*/
#define G       (1.0)             /*Surface Energy*/
#define Mob     (1.0)             /*Mobility*/
#define E       (8.0)            /*epsilon - dimensions of length [m]*/
#define tau     (0.8)
#define radius2  200
#define deltaMu (0.4)
#define Mu      (1.0)
#define Dab     (0.04)      // Strength of Anisotropy
//****************************************************************************//
// #define Corner
// #define Nothing
#define Centre
//****************************************************************************//
int t;
/*************Activate Phi Solver**********************************************/
#define PHI
/*************Activate lbm fluid solver****************************************/
// #define LBM
/******************************************************************************/
#define tol     10e-6
#define tsteps  100000 // no. of iterations
#define savet   500  //file saving steps
#define ftag    106
#ifndef PHI
  #define SMOOTH  1
#endif

#ifdef PHI
  #define SMOOTH 1000
#endif

/******************************************************************************/
// #define ISO
#define ANISO
