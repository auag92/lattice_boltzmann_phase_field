#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "phi_definitions.c"
#include "lbm_definitions.c"
#include "file_io.c"
//****************************************************************************//
void main(){
  int i,j,k;
//*********************Init >>> PhiSolver Stuff ******************************//
  phi_allocate_memory();
  phi_initialize();
  phi_boundary(phi_old);
  phi_boundary(mu_old);
//*********************Init >>> LBM stuff *************************************//
  lbm_allocate_memory();
  lbm_init();
//*****************************************************************************//
  for (t = 1; t <= tsteps; t++) {

    if(t < SMOOTH){
      #ifdef PHI
        phi_solver();
      #endif
    }else{
      #ifdef LBM
        lbm_solver();
      #endif
      #ifdef PHI
        phi_solver();
      #endif
    }
//**********File output routine**************************************************//
    if(t%savet == 0){
      file_output();
    }
    printf("iteration no. %d \n",t);
  }
#ifdef LBM
  lbm_free_memory();
#endif
#ifdef PHI
  phi_free_memory();
#endif
}
