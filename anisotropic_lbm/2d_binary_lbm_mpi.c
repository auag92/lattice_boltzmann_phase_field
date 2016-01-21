#include "stdio.h"
#include "stdlib.h"
#include "math.h"
// #include "constants.h"
#include "phi_definitions.c"
#include "lbm_definitions.c"
//*****************************************************************************//
//*****************************************************************************//
void    file_output();
void    write2file1( int t, int m, double *c, double *d, char *fname);  // for vector array - velocity
void    write2file( int t, int m, double *c, char *fname);              // for scalar array - phi
//-------------------------------------------------------------//
void main(){
  int i,j,k;

  phi_allocate_memory();
  phi_initialize();
  phi_boundary(phi_old);
  phi_boundary(mu_old);

  lbm_allocate_memory();
  lbm_init();

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

void file_output(){
  char fname1[100], fname2[100], fname3[100];
  int  m = sprintf(fname1, "velocities");
  int  n = sprintf(fname2, "rho");
  int  o = sprintf(fname3, "phi");
#ifdef LBM
  write2file1(t,Mx,u,v,fname1);
#endif
#ifdef PHI
  write2file(t,Mx,phi_old,fname3);
#endif
}
void write2file1( int t, int m, double *c, double *d, char *fname) {
  int i,j,z;
  FILE *fp;
  char filename[1000];
  sprintf(filename,"./datafiles%d/%s_%d.dat",ftag, fname, t);
  fp = fopen(filename,"w");
  for ( i = 0; i < m; i++)
  {
    for ( j = 0; j < m; j++)
    {
      z= i*m + j;
      fprintf(fp,"%d %d %le %le\n",j,i,c[z],d[z]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}
void write2file( int t, int m, double *c, char *fname) {
  int i,j,z;
  FILE *fp;
  char filename[1000];
  sprintf(filename,"./datafiles%d/%s_%d.dat",ftag,fname, t);
  fp = fopen(filename,"w");
  for ( i = 0; i < m; i++)
  {
    for ( j=0; j < m; j++)
    {
      z= i*m + j;
      fprintf(fp,"%d %d %le\n",j,i,c[z]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}
