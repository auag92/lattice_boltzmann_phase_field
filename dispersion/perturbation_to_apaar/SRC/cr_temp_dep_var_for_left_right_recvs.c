void cr_temp_dep_var_for_left_right_recvs(struct var *dep_var,struct var *dep_var_fr_left,struct var *dep_var_fr_right, struct domain_info dom_decomp_info,struct domain_info worker_dom_decomp_info,struct worker_info worker)
{

	long i,j;

	long array_index,array_index2cp;

	long array_index_left,array_index_right;
	


	//========================================================================//

	//copying from dep_var_fr_left
	for(j=0;j<worker_dom_decomp_info.rows;j++)
	{
		for(i=0;i<BUFFER_WIDTH;i++)
		{
			array_index = i + worker_dom_decomp_info.cols*j;
			
			array_index_left = i  + BUFFER_WIDTH * j;

			dep_var[array_index]=dep_var_fr_left[array_index_left];

		}
	}

	//copying from comp_fr_right
	for(j=0;j<worker_dom_decomp_info.rows;j++)
	{
		for(i=(dom_decomp_info.cols+BUFFER_WIDTH);i<(dom_decomp_info.cols+2*BUFFER_WIDTH);i++)
		{
			array_index = i + worker_dom_decomp_info.cols*j;
			
			array_index_right = ( i - (dom_decomp_info.cols+BUFFER_WIDTH) ) + BUFFER_WIDTH * j;

			dep_var[array_index]=dep_var_fr_right[array_index_right];

		}
	}
	//=======================================================================//
	
	

	//imposing the no-flux boundary condition 

	//===================================================//
	#ifdef NOFLUX

	//==============================================//
	//horizontal

	
	//at the left end , outermost

	if(worker.num_x==0) //at the left-most boundary 
	{	
		i=0;
		for(j=0;j<worker_dom_decomp_info.rows;j++)
		{
			array_index=i + worker_dom_decomp_info.cols*j;

			array_index2cp=(i+3) + worker_dom_decomp_info.cols*j;
			
			dep_var[array_index]=dep_var[array_index2cp];
		
		}

		//at the left end , one inside
		i=1;
		for(j=0;j<worker_dom_decomp_info.rows;j++)
		{
			array_index=i + worker_dom_decomp_info.cols*j;
		
			array_index2cp=(i+1) + worker_dom_decomp_info.cols*j;
		
			dep_var[array_index]=dep_var[array_index2cp];

		}

	}


	if(worker.num_x ==numworkers.x-1) //at the right most boundary
	{
		//at the right end, the one at the outermost 
		i=worker_dom_decomp_info.cols-1;

		for(j=0;j<worker_dom_decomp_info.rows;j++)
		{	
			array_index=i + worker_dom_decomp_info.cols*j;
			array_index2cp=(i-3) + worker_dom_decomp_info.cols*j;

			dep_var[array_index]=dep_var[array_index2cp];

		}

		//at the right end, the one just inside the outermost 
		i=worker_dom_decomp_info.cols-2;

		for(j=0;j<worker_dom_decomp_info.rows;j++)
		{	
			array_index=i + worker_dom_decomp_info.cols*j;
			array_index2cp=(i-1) + worker_dom_decomp_info.cols*j;

			dep_var[array_index]=dep_var[array_index2cp];

		}
	}
	//================================================//
	#endif
	
}
