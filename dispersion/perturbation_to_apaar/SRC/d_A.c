//the change in amplitude
void diff_amp(struct var *dep_var,double *boundary,long *peak_x,double *peak_y,long *trough_x,double *trough_y)
{
	long array_index,array_index_up;
	int i,j;
	
	
	double bound_inside,bound_outside, fi_inside,fi_outside,epsilon_inter;

	double fi_avg=0.5*(fi_s+fi_l);

	
	//locating the boundary
	//I have tacitly assumed that the cells are growing only in the vertical direction  
	for(i=0;i<array_length_x;i++) //the order of loops is interchanged
	{
		for(j=0;j<array_length_y-1;j++)
		{
			array_index=i+array_length_x*j;		

			array_index_up=i+array_length_x*(j+1);

			if(dep_var[array_index].fi>=fi_avg && dep_var[array_index_up].fi<fi_avg)

			{

				bound_inside=j;
				bound_outside=(j+1);

				fi_inside=dep_var[array_index].fi;
				fi_outside=dep_var[array_index_up].fi;

				epsilon_inter=2.0*(fi_avg-0.5*(fi_inside+fi_outside))/(fi_outside-fi_inside);

				boundary[i]=0.5*(bound_inside*(1.0-epsilon_inter)+bound_outside*(1.0+epsilon_inter));
				break;		
						
								
			}
			
		}
	}

	//calculating the peak and trough
	
	(*peak_y)=0;
	(*trough_y)=array_length_y;

	for(i=0;i<array_length_x;i++) //the order of loops is interchanged
	{
		if(boundary[i]>(*peak_y))
		{
			(*peak_y)=boundary[i];
			(*peak_x)=i;
		}

		if(boundary[i]<(*trough_y))
		{
			(*trough_y)=boundary[i];
			(*trough_x)=i;
		}		
	}
		
	
}
