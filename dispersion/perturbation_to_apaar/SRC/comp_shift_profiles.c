void shift_profile(struct var *dep_var,struct domain_info worker_dom_decomp_info,struct worker_info worker)
{
	long i,j;

	long m;
	
	long array_index;

	long array_index_up;
	
	double mu_init[NUMCOMPONENTS-1];

	//===================================================================================//
	
	//I am going to concern myself with the real boundaries of the system in they direction 
	for(j=2;j<(worker_dom_decomp_info.rows-2);j++)
	{
		for(i=0;i<worker_dom_decomp_info.cols;i++)
		
		{
			//the current location
			array_index=i+worker_dom_decomp_info.cols*j;
	
			//the one on the right
			array_index_up=i+worker_dom_decomp_info.cols*(j+1);

			dep_var[array_index].fi=dep_var[array_index_up].fi;
	
			for(m=0;m<NUMCOMPONENTS-1;m++)
				dep_var[array_index].mu[m]=dep_var[array_index_up].mu[m];
		}
	}

	if(worker.num_y ==numworkers.y-1)//checking whether the worker occupies the topmost layer 
	{
		//recreating the profiles at the last layer (in the real domain)
		j=worker_dom_decomp_info.rows-3;
		
		for(i=0;i<worker_dom_decomp_info.cols;i++)
		{
			//the current location
			array_index=i+worker_dom_decomp_info.cols*j;

			dep_var[array_index].fi=fi_l;

			mu_fr_c(mu_init,c_l_init,0);	

			for(m=0;m<NUMCOMPONENTS-1;m++)
				dep_var[array_index].mu[m]=mu_init[m];
		}

	}

		
	//=========================================================================================//	
}
