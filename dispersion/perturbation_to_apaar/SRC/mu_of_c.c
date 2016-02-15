void mu_fr_c(double *diff_pot,double *c,long phase_index)
{
	#ifdef BINARY

	if(phase_index==0)//this is the liquid
	{
		diff_pot[0]=c[0];
	}
	
	else if(phase_index==1)//this is the solid
	{
		diff_pot[0]=c[0]/K[0];
	}
	#endif	
		

	#ifndef BINARY
	long m,p;

	double d_c_d_mu[NUMCOMPONENTS-1][NUMCOMPONENTS-1];

	double inv_d_c_d_mu[NUMCOMPONENTS-1][NUMCOMPONENTS-1];		

	double c_eq[NUMCOMPONENTS-1],c_min_c_eq[NUMCOMPONENTS-1]; 

	double mu_min_mu_eq[NUMCOMPONENTS-1];

	//getting the d_c_d_mu's
	for(m=0;m<NUMCOMPONENTS-1;m++)
	{
		for(p=0;p<NUMCOMPONENTS-1;p++)
		{
			d_c_d_mu[m][p]=set_d_c_d_mu_bulk(mu_eq,phase_index,m,p);

			
		}
	}

	//d_mu_d_c's
	mat_inv(inv_d_c_d_mu,d_c_d_mu);

	//the equilibrium composition vector
	for(m=0;m<NUMCOMPONENTS-1;m++)
	{
		if(phase_index==1)
			c_eq[m]=c_s_eq[m];

		else if(phase_index==0)
			c_eq[m]=c_l_eq[m];
						
	}

	//the c_min_c_eq vector
	for(m=0;m<NUMCOMPONENTS-1;m++)
	{
		c_min_c_eq[m]=c[m]-c_eq[m];
	}	

	//forming the "mu_min_mu_eq" (see "matrix_operations.h")
	mat_vec_mult(mu_min_mu_eq,inv_d_c_d_mu,c_min_c_eq);

	//forming "mu"
	for(m=0;m<NUMCOMPONENTS-1;m++)
	{
		diff_pot[m]=mu_min_mu_eq[m]+mu_eq[m];
	}
	#endif

}			
