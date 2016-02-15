void comp_tau(int taskid)
{
	long m,n;//the component indices
	
	double prod_1[NUMCOMPONENTS-1];

	double D_L_times_d_c_d_mu_L[NUMCOMPONENTS-1][NUMCOMPONENTS-1],inv_D_L_times_d_c_d_mu_L[NUMCOMPONENTS-1][NUMCOMPONENTS-1];//they will participate in the calculation of tau for anti-trapping 

	double d_c_d_mu_L[NUMCOMPONENTS-1][NUMCOMPONENTS-1];

	

	//==========================================================================================//
	//computing d_c_d_mu for the liquid (see "setting d_c_d_mu_bulk.c")

	for(m=0;m<NUMCOMPONENTS-1;m++)
	{
		for(n=0;n<NUMCOMPONENTS-1;n++)
		{	
			d_c_d_mu_L[m][n]=set_d_c_d_mu_bulk(mu_eq,0,m,n);
		}
	}
	//==========================================================================================//
	 
	//==============================================================================//
	//I am going to print out the diffusivity matrices here (see "matrix_printer.c") 

	if(taskid==MASTER)
	{	
		printf("The d_c_d_mu_L matrix is\n");
		mat_print(d_c_d_mu_L);
	}
	//=============================================================================//


	//=================================================================================//
	//computing the tensors D_L_times_d_c_d_mu_L (see "matrix_matrix_multiplier.c")
	mat_mat_mult(D_L_times_d_c_d_mu_L,D_L,d_c_d_mu_L);	
	
	//===================================================================================//

	
	//===================================================================================//
	//inverting the matrix obtained (see "matrix_inverter.c")
	mat_inv(inv_D_L_times_d_c_d_mu_L,D_L_times_d_c_d_mu_L);	
	
	//=================================================================================//

	//=================================================================================//
	//computing tau	
			
	//first forming prod_1 (see 'matrix_vector_multiplier.c")
	mat_vec_mult(prod_1,inv_D_L_times_d_c_d_mu_L,c_s_min_c_l_eq);	

	//then forming tau (inner product of 2 vectors) (see "vector_inner_product.c")
	tau=vec_inn_prod(prod_1,c_s_min_c_l_eq)*(M_const+F_const)*epsilon;
}
