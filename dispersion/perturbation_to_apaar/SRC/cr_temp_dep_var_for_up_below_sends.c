void cr_temp_dep_var_for_up_below_sends(struct var *dep_var,struct var *dep_var2up,struct var *dep_var2below,struct domain_info dom_decomp_info,struct domain_info worker_dom_decomp_info)
{
	int i,j;

	long array_index,array_index_up,array_index_below;

	//setting comp2up
	for(j=dom_decomp_info.rows;j<(dom_decomp_info.rows+BUFFER_WIDTH);j++)
	{
		for(i=0;i<worker_dom_decomp_info.cols;i++)
		{
			array_index=i+worker_dom_decomp_info.cols*j;
			
			array_index_up= i  + worker_dom_decomp_info.cols * (j-dom_decomp_info.rows);
			dep_var2up[array_index_up]=dep_var[array_index];

		}
	}

	//setting comp2below
	for(j=BUFFER_WIDTH;j<(2*BUFFER_WIDTH);j++)
	{
		for(i=0;i<worker_dom_decomp_info.cols;i++)
		{
			array_index=i+worker_dom_decomp_info.cols*j;
	
			array_index_below= i + worker_dom_decomp_info.cols*(j-BUFFER_WIDTH);

			dep_var2below[array_index_below]=dep_var[array_index];
		}
	}
}
