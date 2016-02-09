#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "phi_definitions.c"

//*****************************************************************************//
void    file_output();
void    write2file1( int t, double *c, double *d, char *fname);  // for vector array - velocity
void    write2file( int t, double *c, char *fname);              // for scalar array - phi
//*****************************************************************************//

void main(){
  int i,j,k;

  phi_allocate_memory();
  phi_initialize();
  phi_boundary(phi_old);
  phi_boundary(mu_old);

  for (t = 1; t <= tsteps; t++) {

    phi_solver();

    if(t%savet == 0){
      file_output();
    }

    printf("iteration no. %d \n",t);
  }
  phi_free_memory();
}

void file_output(){
  char fname1[100], fname2[100], fname3[100];
  int  m = sprintf(fname1, "phi");
  int  n = sprintf(fname2, "mu");
  int  o = sprintf(fname3, "conc");

  write2file(t,phi_old,fname1);
  write2file(t,mu_old,fname2);
  write2file(t,mu_old,fname3);
}
void write2file1( int t, double *c, double *d, char *fname) {
  int i,j,z;
  FILE *fp;
  char filename[1000];
  sprintf(filename,"./datafiles%d/%s_%d.dat",ftag, fname, t);
  fp = fopen(filename,"w");
  for ( i = 0; i < Mx; i++)
  {
    for ( j = 0; j < My; j++)
    {
      z= i*Mx + j;
      fprintf(fp,"%d %d %le %le\n",j,i,c[z],d[z]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}
void write2file( int t, double *c, char *fname) {
  int i,j,z;
  FILE *fp;
  char filename[1000];
  sprintf(filename,"./datafiles%d/%s_%d.dat",ftag,fname, t);
  fp = fopen(filename,"w");
  for ( i = 0; i < Mx; i++)
  {
    for ( j=0; j < My; j++)
    {
      z= i*My + j;
      fprintf(fp,"%d %d %le\n",j,i,c[z]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}
