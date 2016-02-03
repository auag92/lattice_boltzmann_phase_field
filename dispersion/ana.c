#include"stdio.h"
#include"stdlib.h"
#include"math.h"

#define NUMCOMPONENTS 2
#define gamma 0.1

double V;


#define omega_min 0.000001
#define omega_max 0.01
#define omega_step 0.000001


double d_mu_d_c_liq;
double D;
double c_l_eq, c_s_eq, c_s_min_c_l_eq;
double K, G;


//===========================================================================//
//function for computing k_omega
double K_omega(double omega)
{
 double k_omega, first_term, first_term_sq;

 first_term    =  V/(2.0*D);
 first_term_sq =  first_term*first_term;

 k_omega           =  first_term + sqrt(first_term_sq + (omega*omega));

 return(k_omega);
}
//===========================================================================//
//===========================================================================//
//routine for computing omega tilde
double Omega_tilda(double k_omega)
{
  double omega_tilda;
  omega_tilda = k_omega - (V*(1.0-K)/D);
  return(omega_tilda);
}
//===========================================================================//
//===========================================================================//

#include"mat_prop.c"//for setting the setting_d_c_d_mu in the liquid and other constants mentioned above


main(int argc, char *argv[])
{

  char fname[100], str[100];
  FILE *fp,*fp_read,*fp_write;

  long m;
  long i,j;

  double omega;
  double lambda_min, lambda_max, lambda_step;

  double k_omega, omega_tilde;
  double b, rel_b_fir_ter, rel_b_sec_ter;
  double ampl_fac;
  double lhs, rhs;
  //=====================================================================//
  //creating the different property tensors
  //====================================================================//
  // init_mat_prop();
  //=====================================================================//
  //opening the file for writing the amplification files
  //====================================================================//
   sprintf(fname,"ampl_fac_10_D_11_%.1lf_D_22_%.1lf_D_33_%.1lf_calc_3.dat",D);

   if((fp=fopen(fname,"w"))==NULL)
   {
    printf("amplfication profile can't be written\n");
    exit(1);
   }
  //=========================================================================//
  //=========================================================================//
  //opening the file for reading in the equilibrium profiles
   sprintf(fname,"parameters_D_11_%.1lf_D_22_%.1lf_D_33_%.1lf.dat",D_L[0][0],D_L[1][1],D_L[2][2]);

   if((fp_read=fopen(fname,"r"))==NULL)
   {
    printf("parameters can't be read\n");
    exit(1);
   }
   //reading the parameters
   fscanf(fp_read,"%s%lf",str,&V);
   fscanf(fp_read,"%s%lf",str,&(c_s_eq[0]));
   fscanf(fp_read,"%s%lf",str,&(c_l_eq[0]));
   fscanf(fp_read,"%s%lf",str,&(K[0]));
   fscanf(fp_read,"%s%lf",str,&(G[0]));
   fscanf(fp_read,"%s%lf",str,&(c_s_eq[1]));
   fscanf(fp_read,"%s%lf",str,&(c_l_eq[1]));
   fscanf(fp_read,"%s%lf",str,&(K[1]));
   fscanf(fp_read,"%s%lf",str,&(G[1]));
   fscanf(fp_read,"%s%lf",str,&(c_s_eq[2]));
   fscanf(fp_read,"%s%lf",str,&(c_l_eq[2]));
   fscanf(fp_read,"%s%lf",str,&(K[2]));
   fscanf(fp_read,"%s%lf",str,&(G[2]));
   //closing the file
   fclose(fp_read);
   //creating the row_vector
   c_s_min_c_l_eq[0]=c_s_eq[0]-c_l_eq[0];

   c_s_min_c_l_eq[1]=c_s_eq[1]-c_l_eq[1];

   c_s_min_c_l_eq[2]=c_s_eq[2]-c_l_eq[2];
   //=====================================================================//
   //=====================================================================//
   //opening the file for echoing back the parameters
   sprintf(fname,"check_parameters_D_11_%.1lf_D_22_%.1lf_D_33_%.1lf.dat",D_L[0][0],D_L[1][1],D_L[2][2]);
   if((fp_write=fopen(fname,"w"))==NULL)
   {
    printf("parameters can't be written\n");
    exit(1);
   }
   //reading the parameters
   fprintf(fp_write,"V\t%lf\n",V);
   fprintf(fp_write,"c_s_1\t%lf\n",c_s_eq[0]);
   fprintf(fp_write,"c_l_1\t%lf\n",c_l_eq[0]);
   fprintf(fp_write,"K_1\t%lf\n",K[0]);
   fprintf(fp_write,"G_1\t%.12lf\n",G[0]);
   fprintf(fp_write,"c_s_2\t%lf\n",c_s_eq[1]);
   fprintf(fp_write,"c_l_2\t%lf\n",c_l_eq[1]);
   fprintf(fp_write,"K_1\t%lf\n",K[1]);
   fprintf(fp_write,"G_1\t%.12lf\n",G[1]);
   fprintf(fp_write,"c_s_2\t%lf\n",c_s_eq[2]);
   fprintf(fp_write,"c_l_2\t%lf\n",c_l_eq[2]);
   fprintf(fp_write,"K_2\t%lf\n",K[2]);
   fprintf(fp_write,"G_2\t%.12lf\n",G[2]);
   //closing the file
   fclose(fp_write);
  //=====================================================================//
  //=====================================================================//
  //starting the omega loop
  for( omega = omega_min; omega <= omega_max; omega += omega_step)
  {
   //====================================================================//
   //computing omega_tilde and k_omega
   //====================================================================//
    k_omega     = K_omega(omega);
    omega_tilde = Omega_tilda(k_omega);
    b = Gamma*omega*omega/(mu_eq*(1-k));
    //so the amplification factor
    ampl_fac = omega_tilde*V*(-(b/G) + (1.0/omega_tilde)*(k_omega - (V/D_L[0][0])));

    //writing the file
    fprintf(fp,"%lf\t%.10lf\t%lf\t%lf\t%lf\n",omega,ampl_fac,k_omega[0],omega_tilde[0],b[0]);
  }
  fclose(fp);
}
