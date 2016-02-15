void comp_diver_fi(struct grad_info *info,struct var *dep_var,struct domain_info worker_dom_decomp_info,long tsteps)
{
	long i,j;

	long m;

	long array_index;//the current location

	long array_index_below,array_index_left;//the adjacent locations dictated by the stencil

	double div_fi;		

	struct grad_info *grad_current,*grad_below,*grad_left;	

	double driving_force; 

	double d_a_d_fi_current[2],d_a_d_fi_left[2],d_a_d_fi_below[2]; 

	//I am going to stick to the single buffer layer all around boundaries of the domain 
	for(j=1;j<(worker_dom_decomp_info.rows-1);j++)
	{
		for(i=1;i<(worker_dom_decomp_info.cols-1);i++)
		{
			//the current location
			array_index=i+worker_dom_decomp_info.cols*j;

			array_index_below=i+worker_dom_decomp_info.cols*(j-1);
									
			array_index_left=(i-1)+worker_dom_decomp_info.cols*j;

			//=====================================================//

			grad_current=&info[array_index];

			grad_below=&info[array_index_below];

			grad_left=&info[array_index_left];

			//====================================================//
			
			//===================================================//
			
			//computing d_a_d_fi's (see "compute_d_a_d_fi.c")
			
			comp_d_a_d_fi(grad_current,d_a_d_fi_current);
			comp_d_a_d_fi(grad_left,d_a_d_fi_left);
			comp_d_a_d_fi(grad_below,d_a_d_fi_below);

			//computing div_fi

			div_fi=((d_a_d_fi_current[X]-d_a_d_fi_left[X])/dx)+((d_a_d_fi_current[Y]-d_a_d_fi_below[Y])/dy);
			//===================================================//
			
			//===================================================//	
			//computing the driving force
			driving_force=0.0;
			
			for(m=0;m<NUMCOMPONENTS-1;m++)
			{
				driving_force+=c_s_min_c_l_eq[m]*(dep_var[array_index].mu[m]-mu_eq[m]);
			}
			//===================================================//
			
			//===================================================//
			//computing d_fi_d_t
			//if(tsteps<=p_s_tsteps)//killing off the driving force 
			//{						
			//	grad_current->d_fi_d_t=(div_fi-((gamma_surf*9.0*d_f_d_fi(dep_var[array_index].fi))/epsilon)+0.0*driving_force*d_h_d_fi(dep_var[array_index].fi))/(tau*epsilon);
			//}
			//else
			//{						
				grad_current->d_fi_d_t=(div_fi-((gamma_surf*9.0*d_f_d_fi(dep_var[array_index].fi))/epsilon)+driving_force*d_h_d_fi(dep_var[array_index].fi))/(tau*epsilon);
			//}
			//==================================================//

		}
	}
}		
			
														
