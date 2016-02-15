void cr_temp_dep_var_for_up_below_recvs(struct var *dep_var,struct var *dep_var_fr_up,struct var *dep_var_fr_below,struct domain_info dom_decomp_info,struct domain_info worker_dom_decomp_info,struct worker_info worker)
{
	int i,j;
		
	long array_index_up,array_index_below,array_index;

	long array_index2cp;

	//================================================//
	//setting comp_fr_up
	for(j=(dom_decomp_info.rows+BUFFER_WIDTH);j<(dom_decomp_info.rows+2*BUFFER_WIDTH);j++)
	{
		for(i=0;i<worker_dom_decomp_info.cols;i++)
		{
			array_index = i+worker_dom_decomp_info.cols*j;
			
			array_index_up = i + worker_dom_decomp_info.cols * (j-(dom_decomp_info.rows+BUFFER_WIDTH));
			dep_var[array_index]=dep_var_fr_up[array_index_up];

		}
	}

	//setting comp_fr_below
	for(j=0;j<BUFFER_WIDTH;j++)
	{
		for(i=0;i<worker_dom_decomp_info.cols;i++)
		{
			array_index=i+worker_dom_decomp_info.cols*j;
	
			array_index_below= i + worker_dom_decomp_info.cols*j;

			dep_var[array_index]=dep_var_fr_below[array_index_below];
		}
	}

	//===============================================//
	

	//============================================================//
	#ifdef NOFLUX
	//vertical

	if(worker.num_y ==numworkers.y-1)
	{	
		//at the up end, outermost 
		j= worker_dom_decomp_info.rows-1;

		for(i=0;i<worker_dom_decomp_info.cols;i++)
		{
			array_index=i+worker_dom_decomp_info.cols*j;
			array_index2cp=i+worker_dom_decomp_info.cols*(j-3);
	
			dep_var[array_index]=dep_var[array_index2cp];
	
		}

		//at the up end, one inside 
		j= worker_dom_decomp_info.rows-2;

		for(i=0;i<worker_dom_decomp_info.cols;i++)
		{	
			array_index=i+worker_dom_decomp_info.cols*j;
			array_index2cp=i+worker_dom_decomp_info.cols*(j-1);

			dep_var[array_index]=dep_var[array_index2cp];
	
									
		}

	}

	if(worker.num_y==0)
	{
		//at the bottom end, outermost
		j=0;
		for(i=0;i<worker_dom_decomp_info.cols;i++)
		{
			array_index=i+worker_dom_decomp_info.cols*j;
			array_index2cp=i+worker_dom_decomp_info.cols*(j+3);

			dep_var[array_index]=dep_var[array_index2cp];

		
		}

		//at the bottom end, one inside
		j=1;
		for(i=0;i<worker_dom_decomp_info.cols;i++)
		{
			array_index=i+worker_dom_decomp_info.cols*j;
			array_index2cp=i+worker_dom_decomp_info.cols*(j+1);

			dep_var[array_index]=dep_var[array_index2cp];
		
		
		}
	}

	#endif	
	//==============================================================//	

}
