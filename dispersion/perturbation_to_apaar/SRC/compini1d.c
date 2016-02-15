void compini(struct var *dep_var)
{

	long i,j;

	long m;
	
	long array_index;
	
	FILE *fp;
	char fname[300];

	
	double distance;

	double mu_init[NUMCOMPONENTS-1];

	//====================================================================//
	sprintf(fname,"../DATA/initial_profile.dat");
	if((fp=fopen(fname,"w"))==NULL)
	{
		printf("Initial profile can't be written\n");
		exit(1);
	}
	//====================================================================//

	//====================================================================//
	

	
	//==================================================================//	

	//creating a one D profile		
	for(j=0;j< array_length_y;j++)
	{
		for(i=0;i<array_length_x;i++)
		{
			array_index=i+array_length_x*j;	

			//computing the distance from the corner

			if(j<=solid_phase_length-BUFFER_WIDTH)//this region is the solid
			{
				dep_var[array_index].fi=fi_s;
			
				mu_fr_c(mu_init,c_s_init,1);					

				for(m=0;m<NUMCOMPONENTS-1;m++)
				{
					dep_var[array_index].mu[m]=mu_init[m];
				}
					
			}
		
			else//this region is the liquid
			{
				dep_var[array_index].fi=fi_l;

				mu_fr_c(mu_init,c_l_init,0);					

				for(m=0;m<NUMCOMPONENTS-1;m++)
				{
					dep_var[array_index].mu[m]=mu_init[m];
				}
			}

			#ifdef BINARY 
			
			fprintf(fp,"%ld\t%ld\t%lf\t%lf\n",i,j,dep_var[array_index].fi,dep_var[array_index].mu[0]);
			#endif 	
	
			#ifdef TERNARY 
			
			fprintf(fp,"%ld\t%ld\t%lf\t%lf\t%lf\n",i,j,dep_var[array_index].fi,dep_var[array_index].mu[0],dep_var[array_index].mu[1]);
			#endif 	

			#ifdef QUARTERNARY 
			
			fprintf(fp,"%ld\t%ld\t%lf\t%lf\t%lf\t%lf\n",i,j,dep_var[array_index].fi,dep_var[array_index].mu[0],dep_var[array_index].mu[1],dep_var[array_index].mu[2]);
			#endif 	
			
		}
		fprintf(fp,"\n");
	}

	fclose(fp);

	
}
