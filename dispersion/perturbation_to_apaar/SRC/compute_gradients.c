void comp_grad(struct grad_info *info,struct var *dep_var,struct domain_info worker_dom_decomp_info)
{
	long i,j;//for iterating over the space

	long m,n;//for iterating over the no. of components 

	long array_index;//the current location

	long array_index_left,array_index_right;//for the central diffrence to find fi_x at the current location

	long array_index_up_left,array_index_up_right;//for the central diffrence to find fi_x at the adjacent location vertically up

	long array_index_up,array_index_below;//for the central diffrence to find fi_y at the current location

	long array_index_right_up,array_index_right_below;//for the central diffrence to find fi_y at the adjacent location
	
	struct grad_info *grad_current; 
	

	for(j=0;j<worker_dom_decomp_info.rows ;j++)
	{
		for(i=0;i<worker_dom_decomp_info.cols ;i++)
		{

			//=======================================================//
			//computing the various array indices

			array_index=i+worker_dom_decomp_info.cols*j;

			array_index_left=(i-1)+worker_dom_decomp_info.cols*j;
			
			array_index_right=(i+1)+worker_dom_decomp_info.cols*j;
	
			array_index_up_right=(i+1)+worker_dom_decomp_info.cols*(j+1);

			array_index_up_left=(i-1)+worker_dom_decomp_info.cols*(j+1);	

			array_index_up=i+worker_dom_decomp_info.cols*(j+1);

			array_index_below=i+worker_dom_decomp_info.cols*(j-1);

			array_index_right_up=(i+1)+worker_dom_decomp_info.cols*(j+1);

			array_index_right_below=(i+1)+worker_dom_decomp_info.cols*(j-1);
			//=====================================================//

			//=====================================================//
			grad_current=&info[array_index];			
			//====================================================//

			//================================================================================//
			//computing fi_x which has to be used in divergence along x also x component of mu gradient and M_mid along x
			if((i+1)<=worker_dom_decomp_info.cols-1)
			{
				//for fi
				grad_current->grad_fi[X]=(dep_var[array_index_right].fi-dep_var[array_index].fi)/dx;
				//for mu
				for(m=0;m<NUMCOMPONENTS-1;m++)
				{
					grad_current->grad_mu[X][m]=(dep_var[array_index_right].mu[m]-dep_var[array_index].mu[m])/dx;	
					//forming the mobility tensor
					for(n=0;n<NUMCOMPONENTS-1;n++)
					{
						grad_current->M_mid[X][m][n]=0.5*(form_M(dep_var[array_index],m,n)+form_M(dep_var[array_index_right],m,n));//see "form_mobility.c"
					}
				}	 
			}

			//============================================================================//


			//============================================================================//
			//computing fi_x which has to be used in divergence along y
			if((i+1)<=worker_dom_decomp_info.cols-1 && (i-1)>=0 && (j+1)<=worker_dom_decomp_info.rows-1)
			{
				//for fi
				grad_current->grad_fi_central[X]=0.5*(((dep_var[array_index_right].fi-dep_var[array_index_left].fi)/(2.0*dx))+((dep_var[array_index_up_right].fi-dep_var[array_index_up_left].fi)/(2.0*dx)));
			}
			//============================================================================//

	
			//============================================================================//
			//computing fi_y which has to be used in divergence along y
			if((j+1)<=worker_dom_decomp_info.rows-1)
			{
				//for fi
				grad_current->grad_fi[Y]=(dep_var[array_index_up].fi-dep_var[array_index].fi)/dy;
				//for mu
				for(m=0;m<NUMCOMPONENTS-1;m++)
				{
					grad_current->grad_mu[Y][m]=(dep_var[array_index_up].mu[m]-dep_var[array_index].mu[m])/dy;	
					//forming the mobility tensor
					for(n=0;n<NUMCOMPONENTS-1;n++)
					{
						grad_current->M_mid[Y][m][n]=0.5*(form_M(dep_var[array_index],m,n)+form_M(dep_var[array_index_up],m,n));//see "form_mobility.c"
					}
				}	 
			}
			//==========================================================================//	
		
			//==========================================================================//
			//computing fi_y which has to be used in divergence along x
			if((i+1)<=worker_dom_decomp_info.cols-1 && (j-1)>=0 && (j+1)<=worker_dom_decomp_info.rows-1)
			{
				//for fi
				grad_current->grad_fi_central[Y]=0.5*(((dep_var[array_index_up].fi-dep_var[array_index_below].fi)/(2.0*dy))+((dep_var[array_index_right_up].fi-dep_var[array_index_right_below].fi)/(2.0*dy)));
			}
			//=========================================================================//
		}
	}	
}			

		


