void comp_anti_trap(struct grad_info *info,struct var *dep_var,struct domain_info worker_dom_decomp_info)
{
	
	long i,j;

	long m;

	long array_index;//the current location

	long array_index_right,array_index_up;

	struct grad_info *grad_current,*grad_up,*grad_right;	

	double g_avg[2],fi_times_one_min_fi_avg[2],d_fi_d_t_avg[2],one_min_h_avg[2];

	double mod_grad_fi[2];

	double  c_s[NUMCOMPONENTS-1],c_l[NUMCOMPONENTS-1],c_l_min_c_s[NUMCOMPONENTS-1]; //the values of c_s and c_l corresponding to the value of mu at that particular location 

	double c_s_right[NUMCOMPONENTS-1],c_l_right[NUMCOMPONENTS-1],c_l_min_c_s_right[NUMCOMPONENTS-1]; //the values of c_s and c_l corresponding to the value of mu at a location which is just right to the point in question

	double c_s_up[NUMCOMPONENTS-1],c_l_up[NUMCOMPONENTS-1],c_l_min_c_s_up[NUMCOMPONENTS-1]; //the values of c_s and c_l corresponding to the value of mu at a location which is just above the point in question

	double c_l_min_c_s_avg[2][NUMCOMPONENTS-1];




	//I will stick to a single buffer layer all around but won't be using the values at the boundary 
	for(j=1;j<(worker_dom_decomp_info.rows-1);j++)
	{
		for(i=1;i<(worker_dom_decomp_info.cols-1);i++)
		{
			//====================================================//
			//the array indices
			array_index=i+worker_dom_decomp_info.cols*j;	

			array_index_right=(i+1)+worker_dom_decomp_info.cols*j;

			array_index_up=i+worker_dom_decomp_info.cols*(j+1);
	
			//=====================================================//

			//======================================================//
			grad_current=&info[array_index];			

			grad_right=&info[array_index_right];	

			grad_up=&info[array_index_up];			
	
			//====================================================//
			
			
			//=================================================================//
			//computing mod_grad_fi
			//this is second order accurate in between the nodes 

			mod_grad_fi[X]=sqrt(grad_current->grad_fi[X]*grad_current->grad_fi[X]+grad_current->grad_fi_central[Y]*grad_current->grad_fi_central[Y]);

			mod_grad_fi[Y]=sqrt(grad_current->grad_fi_central[X]*grad_current->grad_fi_central[X]+grad_current->grad_fi[Y]*grad_current->grad_fi[Y]);

			//=================================================================//

			//================================================================
					
			//getting c_s and c_l(see "c_of_mu.c")

			c_fr_mu(dep_var[array_index].mu,c_s,c_l);

			c_fr_mu(dep_var[array_index_right].mu,c_s_right,c_l_right);

			c_fr_mu(dep_var[array_index_up].mu,c_s_up,c_l_up);	
			

			//==================================================//
			
			//==================================================//
			//computing c_l_min_c_s
			for(m=0;m<NUMCOMPONENTS-1;m++)
			{ 		
				c_l_min_c_s[m]=c_l[m]-c_s[m];

				c_l_min_c_s_right[m]=c_l_right[m]-c_s_right[m];

				c_l_min_c_s_up[m]=c_l_up[m]-c_s_up[m];
			
			}
			//==================================================//

			//==================================================//
			//computing the average
			for(m=0;m<NUMCOMPONENTS-1;m++)
			{ 
				c_l_min_c_s_avg[X][m]=0.5*(c_l_min_c_s[m]+c_l_min_c_s_right[m]);	

				c_l_min_c_s_avg[Y][m]=0.5*(c_l_min_c_s[m]+c_l_min_c_s_up[m]);	
					
			}
			//==================================================//

			//=================================================================//
			//computing the average of g 

			g_avg[X]=0.5*(g(dep_var[array_index].fi)+g(dep_var[array_index_right].fi));

			g_avg[Y]=0.5*(g(dep_var[array_index].fi)+g(dep_var[array_index_up].fi));

			//=================================================================//

			//=================================================================//
			//computing the average of (1.0-h)
		
			one_min_h_avg[X]=1.0-0.5*(h(dep_var[array_index].fi)+h(dep_var[array_index_right].fi));

			one_min_h_avg[Y]=1.0-0.5*(h(dep_var[array_index].fi)+h(dep_var[array_index_up].fi));

			//=================================================================//

			//=================================================================//
			//computing the average of fi_times_one_min_fi

			fi_times_one_min_fi_avg[X]=0.5*(fi_times_one_min_fi(dep_var[array_index].fi)+fi_times_one_min_fi(dep_var[array_index_right].fi));

			fi_times_one_min_fi_avg[Y]=0.5*(fi_times_one_min_fi(dep_var[array_index].fi)+fi_times_one_min_fi(dep_var[array_index_up].fi));
			//=================================================================// 

			//=================================================================//
			//computing d_fi_d_t_avg
			d_fi_d_t_avg[X]=0.5*(grad_current->d_fi_d_t+grad_right->d_fi_d_t);

			d_fi_d_t_avg[Y]=0.5*(grad_current->d_fi_d_t+grad_up->d_fi_d_t);
			//=================================================================//
				
			//=================================================================//	
			//computing the anti-trapping current
			#ifdef ANTITRAPPING
			if(fabs(fi_times_one_min_fi_avg[X])>1e-12 )
			{
				for(m=0;m<NUMCOMPONENTS-1;m++)
				{ 
					grad_current->j_at[X][m]=-(epsilon*g_avg[X]*one_min_h_avg[X]*c_l_min_c_s_avg[X][m]*d_fi_d_t_avg[X]*(grad_current->grad_fi[X]/mod_grad_fi[X]))/(3.0*fi_times_one_min_fi_avg[X]);

					
				}
			}

			else
			{
				for(m=0;m<NUMCOMPONENTS-1;m++)
				{ 
					grad_current->j_at[X][m]=0.0;

				}
			}


			if( fabs(fi_times_one_min_fi_avg[Y])>1e-12)
			{
				for(m=0;m<NUMCOMPONENTS-1;m++)
				{ 
					

					grad_current->j_at[Y][m]=-(epsilon*g_avg[Y]*one_min_h_avg[Y]*c_l_min_c_s_avg[Y][m]*d_fi_d_t_avg[Y]*(grad_current->grad_fi[Y]/mod_grad_fi[Y]))/(3.0*fi_times_one_min_fi_avg[Y]);
				}
			}

			else
			{
				for(m=0;m<NUMCOMPONENTS-1;m++)
				{ 
					

					grad_current->j_at[Y][m]=0.0;
				}
			}
			#endif
	
			#ifdef NO_ANTITRAPPING
			for(m=0;m<NUMCOMPONENTS-1;m++)
			{ 
				grad_current->j_at[X][m]=0.0;

				grad_current->j_at[Y][m]=0.0;
			}
			#endif
			//=================================================================//
		}
	}

}	
