int ch_shift(struct var *dep_var,struct offset_info offset,struct domain_info dom_decomp_info,struct domain_info worker_dom_decomp_info)
{
	long global_y_bottom, global_y_top; 

	long i,j;

	int shift_on=0;

	long array_index;	

	//I am shifting only along the y direction
	global_y_bottom=offset.y;
		
	global_y_top=offset.y + dom_decomp_info.rows;		

	//checking whether the shift position falls under the purview of the current worker
	if(shift_length>global_y_top || shift_length<global_y_bottom)	
	{
		return (0);

	}
	
	else
	{
		j=shift_length-global_y_bottom+BUFFER_WIDTH;

		//considering the real boundaries of the system	
		for(i=2;i<(worker_dom_decomp_info.cols-2);i++)
		{	
		
			array_index=i+worker_dom_decomp_info.cols*j;

			if(dep_var[array_index].fi>0.5)
			{
				shift_on=1;
				
				break;
			}			
			
		}

		return (shift_on);
	}		
}
