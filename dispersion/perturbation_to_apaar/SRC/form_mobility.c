double form_M(struct var dep_var,long m,long n)
{
	long p;
		
	double val_l=0.0,val_s=0.0;

		
	for(p=0;p<NUMCOMPONENTS-1;p++)
	{
		val_l+=D_L[m][p]*set_d_c_d_mu_bulk(dep_var.mu,0,p,n); //see "setting d_c_d_mu_bulk.c"

		val_s+=D_S[m][p]*set_d_c_d_mu_bulk(dep_var.mu,1,p,n); //see "setting d_c_d_mu_bulk.c"
	}
			
	return(val_l*(1.0-g(dep_var.fi))+val_s*g(dep_var.fi));
}
