void cr_temp_dep_var_for_left_right_sends(struct var *dep_var,struct var *dep_var2left,struct var *dep_var2right,struct domain_info dom_decomp_info,struct domain_info worker_dom_decomp_info)
{
	long i,j;

		
	long array_index, array_index_left,array_index_right;

	
	//setting dep_var2left
	for(j=0;j<worker_dom_decomp_info.rows;j++)
	{
		for(i=BUFFER_WIDTH;i<2*BUFFER_WIDTH;i++)
		{
			array_index = i + worker_dom_decomp_info.cols*j;
			
			array_index_left = ( i - BUFFER_WIDTH) +BUFFER_WIDTH * j;

			dep_var2left[array_index_left]=dep_var[array_index];

		}
	}	 

	//setting comp2right
	for(j=0;j<worker_dom_decomp_info.rows;j++)
	{
		for(i=dom_decomp_info.cols;i<(dom_decomp_info.cols+BUFFER_WIDTH);i++)
		{
			array_index = i + worker_dom_decomp_info.cols*j;
			
			array_index_right = ( i - dom_decomp_info.cols ) + BUFFER_WIDTH * j;

			dep_var2right[array_index_right]=dep_var[array_index];

		}
	}
	
	

}
