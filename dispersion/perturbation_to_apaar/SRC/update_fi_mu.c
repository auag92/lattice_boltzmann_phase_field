void upd_mu_fi(struct grad_info *info,struct var *dep_var,struct domain_info worker_dom_decomp_info,long tsteps)
{
	long i,j;

	long m,p;

	long array_index,array_index_below,array_index_left;

	struct grad_info *grad_current,*grad_below,*grad_left;

	double div_comp[NUMCOMPONENTS-1];

	double flux_in_c[2][NUMCOMPONENTS-1],flux_in_c_left[2][NUMCOMPONENTS-1],flux_in_c_below[2][NUMCOMPONENTS-1];
	
	double tot_flux[2][NUMCOMPONENTS-1],tot_flux_left[2][NUMCOMPONENTS-1],tot_flux_below[2][NUMCOMPONENTS-1];

	double  c_s[NUMCOMPONENTS-1],c_l[NUMCOMPONENTS-1],c_s_min_c_l[NUMCOMPONENTS-1]; //the values of c_s and c_l corresponding to the value of mu at that particular location 

	double RHS[NUMCOMPONENTS-1];
	
	double LHS[NUMCOMPONENTS-1][NUMCOMPONENTS-1],inv_LHS[NUMCOMPONENTS-1][NUMCOMPONENTS-1];

	double d_mu_d_t[NUMCOMPONENTS-1];

	double interp,interp_below,interp_left; 

	//I am going to concern myself with the real boundaries of the system 
	for(j=2;j<(worker_dom_decomp_info.rows-2);j++)
	{
		for(i=2;i<(worker_dom_decomp_info.cols-2);i++)
		{
			//======================================//
			//the locations
			array_index=i+worker_dom_decomp_info.cols*j;

			array_index_below=i+worker_dom_decomp_info.cols*(j-1);
									
			array_index_left=(i-1)+worker_dom_decomp_info.cols*j;
			//======================================//

			//=====================================================//

			grad_current=&info[array_index];

			grad_below=&info[array_index_below];

			grad_left=&info[array_index_left];

			//====================================================//

			//====================================================//
			//the fluxes in c
			
			for(m=0;m<NUMCOMPONENTS-1;m++)
			{ 

				flux_in_c[X][m]=0.0;

				flux_in_c[Y][m]=0.0;



				flux_in_c_left[X][m]=0.0;

				flux_in_c_left[Y][m]=0.0;



				flux_in_c_below[X][m]=0.0;

				flux_in_c_below[Y][m]=0.0;



				for(p=0;p<NUMCOMPONENTS-1;p++)
				{
					flux_in_c[X][m]+=grad_current->M_mid[X][m][p]*grad_current->grad_mu[X][p];
	
					flux_in_c[Y][m]+=grad_current->M_mid[Y][m][p]*grad_current->grad_mu[Y][p];



					flux_in_c_left[X][m]+=grad_left->M_mid[X][m][p]*grad_left->grad_mu[X][p];
	
					flux_in_c_left[Y][m]+=grad_left->M_mid[Y][m][p]*grad_left->grad_mu[Y][p];

					
					flux_in_c_below[X][m]+=grad_below->M_mid[X][m][p]*grad_below->grad_mu[X][p];
	
					flux_in_c_below[Y][m]+=grad_below->M_mid[Y][m][p]*grad_below->grad_mu[Y][p];


				}

				
				if(fabs(dep_var[array_index].fi)>0.0001 && fabs(dep_var[array_index].fi)<0.9999)
				{
					tot_flux[X][m]=-flux_in_c[X][m]+grad_current->j_at[X][m];

					tot_flux[Y][m]=-flux_in_c[Y][m]+grad_current->j_at[Y][m];
	
				}

				else
				{
					tot_flux[X][m]=-flux_in_c[X][m];

					tot_flux[Y][m]=-flux_in_c[Y][m];
				}


				if(fabs(dep_var[array_index_left].fi)>0.0001 && fabs(dep_var[array_index_left].fi)<0.9999)
				{
	
					tot_flux_left[X][m]=-flux_in_c_left[X][m]+grad_left->j_at[X][m];
	
					tot_flux_left[Y][m]=-flux_in_c_left[Y][m]+grad_left->j_at[Y][m];

				}

				else
				{
					tot_flux_left[X][m]=-flux_in_c_left[X][m];
	
					tot_flux_left[Y][m]=-flux_in_c_left[Y][m];
				}							

				if(fabs(dep_var[array_index_below].fi)>0.0001 && fabs(dep_var[array_index_below].fi)<0.9999)
				{

					tot_flux_below[X][m]=-flux_in_c_below[X][m]+grad_below->j_at[X][m];

					tot_flux_below[Y][m]=-flux_in_c_below[Y][m]+grad_below->j_at[Y][m];
				}

				else
				{
					tot_flux_below[X][m]=-flux_in_c_below[X][m];

					tot_flux_below[Y][m]=-flux_in_c_below[Y][m];
				}	

				//adding the noise in the composition fluxes

				//at the current location
				interp=dep_var[array_index].fi*dep_var[array_index].fi*(1.0-dep_var[array_index].fi)*(1.0-dep_var[array_index].fi);

				tot_flux[X][m]+=interp*noise_amplitude_2*(2.0*ran2(&SEED_2)-1.0);

				tot_flux[Y][m]+=interp*noise_amplitude_2*(2.0*ran2(&SEED_2)-1.0);

				//at the left location
				interp_left=dep_var[array_index_left].fi*dep_var[array_index_left].fi*(1.0-dep_var[array_index_left].fi)*(1.0-dep_var[array_index_left].fi);

				tot_flux_left[X][m]+=interp_left*noise_amplitude_2*(2.0*ran2(&SEED_2)-1.0);

				tot_flux_left[Y][m]+=interp_left*noise_amplitude_2*(2.0*ran2(&SEED_2)-1.0);
						
				//at the below location
				interp_below=dep_var[array_index_below].fi*dep_var[array_index_below].fi*(1.0-dep_var[array_index_below].fi)*(1.0-dep_var[array_index_below].fi);

				tot_flux_below[X][m]+=interp_below*noise_amplitude_2*(2.0*ran2(&SEED_2)-1.0);

				tot_flux_below[Y][m]+=interp_below*noise_amplitude_2*(2.0*ran2(&SEED_2)-1.0);							

				//computing the divergences

				div_comp[m]=((tot_flux[X][m]-tot_flux_left[X][m])/dx)+((tot_flux[Y][m]-tot_flux_below[Y][m])/dy);
				

			}
			//====================================================================//

			//getting c_s and c_l
			c_fr_mu(dep_var[array_index].mu,c_s,c_l);

			//forming the RHS
			for(m=0;m<NUMCOMPONENTS-1;m++)
			{ 		
				c_s_min_c_l[m]=c_s[m]-c_l[m];

				RHS[m]=-div_comp[m]-c_s_min_c_l[m]*d_h_d_fi(dep_var[array_index].fi)*grad_current->d_fi_d_t;
						
					
						
			}	

			//forming the LHS	
			for(m=0;m<NUMCOMPONENTS-1;m++)
			{ 
				for(p=0;p<NUMCOMPONENTS-1;p++)
				{
					LHS[m][p]=set_d_c_d_mu_bulk(dep_var[array_index].mu,0,m,p)*(1.0-h(dep_var[array_index].fi)) + set_d_c_d_mu_bulk(dep_var[array_index].mu,1,m,p)*h(dep_var[array_index].fi);	
				}
			}
		
			//inverting the LHS 	
			mat_inv(inv_LHS, LHS);

			//computing "d_mu_d_t" and updating 'mu'
			for(m=0;m<NUMCOMPONENTS-1;m++)
			{
				d_mu_d_t[m]=0.0;
				
				for(p=0;p<NUMCOMPONENTS-1;p++)
				{
					d_mu_d_t[m]+=inv_LHS[m][p]*RHS[p];
				}
				
				//if(tsteps>p_s_tsteps)
				dep_var[array_index].mu[m]+=(dt*d_mu_d_t[m]);
			}

			//updating 'fi'
			dep_var[array_index].fi+=(dt*grad_current->d_fi_d_t);

		}

	}
}	
						
