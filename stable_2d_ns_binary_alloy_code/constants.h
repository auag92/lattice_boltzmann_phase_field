#define pmesh               600                     //Mesh dimension for the pressure grid
#define pmesh2              (pmesh*pmesh)
#define MESHX               (pmesh + 1)             //Mesh dimension for the phi, U and V grids
#define MESHX2              (MESHX*MESHX)
#define deltax              1.0
#define inv_deltax          (1.0/deltax)
#define inv_deltax2         (inv_deltax*inv_deltax)
#define deltat              0.01
//----------Constants related to solidification equations-----------------------
#define K         (0.5)       // Partition Coefficient
#define G         (1.0)       // Surface Energy
#define Mob       (1.0)       // Mobility
#define E         (4.0)       // relaxation length - dimensions of length [m]
#define tau       (1.0)       // relaxation time
#define deltaMu   (0.4)       // driving force in terms of deviation from mu_equilibriu,
#define Mu        (1.0)       // chemical potential
#define Dab       (0.04)      // Strength of Anisotropy
//-------Switches for toggling between the output configuration-----------------
#define           growth          //for activating phi growth
#define           FLUID           // for activating fluid solver
#define           Centre          // seed in centre
// #define        Corner          // seed in corner
// #define        Nothing         // no seed

#define           radius2   1000  // size of the seed
// #define ISO                       // For activating Isotropic Solver
#define ANISO                     // For activating Anisotropic Solver

//------------------------------------------------------------------------------
#define GS_MPI
// #define GS_UNI                    //Serial GS solver
//-------File output parameters-------------------------------------------------
#define save_phi           (100)
#define save_fluid         (100)
#define phi_timesteps     (300000)
#define ftag                 3
//---------Tolerences-----------------------------------------------------------
#ifdef growth
  #define SMOOTH             (10000)
#endif
#ifndef growth
#define SMOOTH              (1)
#endif
#define phi_tol             0.5
#define gs_tol             10e-6
#define dirichlet_pressure
// #define LDC
#define PipeFlow
//--------------Defnitions for fluid flow---------------------------------------
#ifdef PipeFlow
#define inv_Re   (1)
#define Ue      (0.0) // east wall velocity
#define Uw      (0.0) // west wall velocity
#define Un      (0.0) // north wall velocity
#define Us      (0.0) // south wall velocity
#define Ve      (0.0) // east wall velocity
#define Vw      (0.0) // west wall velocity
#define Vn      (0.0) // north wall velocity
#define Vs      (0.0) // south wall velocity
//------------------------------------------------------------------------------
#define p_up      (0)
#define p_down    (0)
#define p_left    (5)
#define p_right   (0)
#endif
//------------------------------------------------------------------------------
#ifdef LDC
#define inv_Re (1)
#define Ue (0.0) // east wall velocity
#define Uw (0.0) // west wall velocity
#define Un (3.0) // north wall velocity
#define Us (0.0) // south wall velocity
#define Ve (0.0) // east wall velocity
#define Vw (0.0) // west wall velocity
#define Vn (0.0) // north wall velocity
#define Vs (0.0) // south wall velocity
//----------------------------------------
#define p_up      (0)
#define p_down    (0)
#define p_left    (0)
#define p_right   (0)
#endif
//------------------------------------------------------------------------------
