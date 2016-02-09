
/*********************Parameters for Phi Solver*******************************/
#define Mx     600
#define My     100
#define MESHX  Mx
#define MESHX2 (Mx*My)
#define deltat (0.02)
#define deltax (1.0)
#define inv_deltax (1./deltax)
#define deltax2 (deltax*deltax)

#define K       (0.2)             /*Partition Coefficient*/
#define G       (1.0)             /*Surface Energy*/
#define Mob     (1.0)             /*Mobility*/
#define E       (4.0)            /*epsilon - dimensions of length [m]*/
#define tau     (1.0)
#define deltaMu (0.4)
#define Mu      (1.0)
#define Dab     (0.04)          // Strength of Anisotropy
//**********************Initian Configurati**********************************//
// #define Corner
// #define Nothing
// #define Centre
#ifdef Centre
  #define RADIUS2  200
#endif
// #define SINE
#define PlaneFront
//****************************************************************************//
int t;
/*************Activate Phi Solver**********************************************/
#define PHI
/******************************************************************************/
#define tsteps  100000 // no. of iterations
#define savet   500  //file saving steps
#define ftag    101
/*****************Toggle Isotropy-Anisotropy*************************************************************/
#define ISO
// #define ANISO
