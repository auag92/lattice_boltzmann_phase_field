void c_fr_mu(double *diff_pot,double *c_s,double *c_l)
{
	#ifndef BINARY
	long m,n;	
	#endif

	#ifdef BINARY
	c_l[0]=diff_pot[0];

	c_s[0]=K[0]*c_l[0];
	#endif

	#ifndef BINARY
	for(m=0;m<NUMCOMPONENTS-1;m++)
	{
		c_s[m]=c_s_eq[m];
		
		c_l[m]=c_l_eq[m];

		for(n=0;n<NUMCOMPONENTS-1;n++)
		{	
			c_l[m]+=set_d_c_d_mu_bulk(diff_pot,0,m,n)*(diff_pot[n]-mu_eq[n]);

			c_s[m]+=set_d_c_d_mu_bulk(diff_pot,1,m,n)*(diff_pot[n]-mu_eq[n]);
		}
	}
	#endif	
}		
		
