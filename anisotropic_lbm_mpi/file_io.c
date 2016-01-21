//*****************************************************************************//
void    file_output();
void    write2file1( int t, int m, double *c, double *d, char *fname);  // for vector array - velocity
void    write2file( int t, int m, double *c, char *fname);              // for scalar array - phi
//*****************************************************************************//
void file_output(){
  char fname1[100], fname2[100], fname3[100];
  int  m = sprintf(fname1, "velocities");
  int  n = sprintf(fname2, "mu");
  int  o = sprintf(fname3, "phi");
#ifdef LBM
  write2file1(t,Mx,u,v,fname1);
#endif
#ifdef PHI
  write2file(t,Mx,mu_old, fname2);
  write2file(t,Mx,phi_old,fname3);
#endif
}
//****************************************************************************//
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
//*****************************************************************************//
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
//*****************************************************************************//
