void write_out(struct var *dep_var,double time)
{
	long i,j;//indices to iterate over the space

	long m;//index to go over the no. of components

	long array_index;

	double c_s[NUMCOMPONENTS-1],c_l[NUMCOMPONENTS-1],c[NUMCOMPONENTS-1];	

	char fname[300];

	printf("File writing will commence at time=%lf\n",time);


	//for the 2D files
	FILE *fp;

	sprintf(fname,"../DATA/sup_sat_%.5lf_%.5lf_%.5lf_w_l_%ld/phase_field_profiles_%.5lf.dat",c_l_eq[0]-c_l_init[0],c_l_eq[1]-c_l_init[1],c_l_eq[2]-c_l_init[2],sine_wave_length,time);
	if((fp=fopen(fname,"w"))==NULL)
	{
		printf("phase field profile can't be written\n");
		exit(1);
	}

	//for 1d files
	#ifdef ONE_D
	FILE *fp_1d,*fp_shift;

	sprintf(fname,"../DATA/sup_sat_%.5lf_%.5lf_%.5lf_w_l_%ld/one_D_phase_field_profiles_%.5lf.dat",c_l_eq[0]-c_l_init[0],c_l_eq[1]-c_l_init[1],c_l_eq[2]-c_l_init[2],sine_wave_length,time);
	if((fp_1d=fopen(fname,"w"))==NULL)
	{
		printf("one D phase field profile can't be written\n");
		exit(1);
	}

	sprintf(fname,"../DATA/sup_sat_%.5lf_%.5lf_%.5lf_w_l_%ld/shift.dat",c_l_eq[0]-c_l_init[0],c_l_eq[1]-c_l_init[1],c_l_eq[2]-c_l_init[2],sine_wave_length);
	if((fp_shift=fopen(fname,"a"))==NULL)
	{
		printf("shift profile can't be written\n");
		exit(1);
	}

	fprintf(fp_shift,"%lf\t%ld\n",time,shift_count);
	fclose(fp_shift);
	#endif
	
	//I am going to stick to the real boundaries of the domain 
	for(j=0;j<array_length_y;j++)
	{
		for(i=0;i<array_length_x;i++)
		{
			//the current location
			array_index=i+array_length_x*j;

			//================================================================//
			//================================================================//
			//computing c_s and c_l

			//at the current location
			//==================================================//
			//getting c_s and c_l(see "c_of_mu.c")
			//for the current location
			c_fr_mu(dep_var[array_index].mu,c_s,c_l);
			//==================================================//

			//==================================================//
			//computing 'c'
			for(m=0;m<NUMCOMPONENTS-1;m++)
			{
				c[m]=h(dep_var[array_index].fi)*c_s[m]+(1.0-h(dep_var[array_index].fi))*c_l[m];
			}
			//====================================================//

			//==================================================================//
			//==================================================================//		

			#ifdef BINARY
			fprintf(fp,"%ld\t%ld\t%lf\t%lf\t%lf\n",i,j,dep_var[array_index].fi,dep_var[array_index].mu[0],c[0]);
			#endif
			
			#ifdef TERNARY
			fprintf(fp,"%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\n",i,j,dep_var[array_index].fi,dep_var[array_index].mu[0],dep_var[array_index].mu[1],c[0],c[1]);
			#endif
	
			#ifdef QUARTERNARY
			fprintf(fp,"%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",i,j,dep_var[array_index].fi,dep_var[array_index].mu[0],dep_var[array_index].mu[1],dep_var[array_index].mu[2],c[0],c[1],c[2]);
			#endif
			
			#ifdef ONE_D
			if(i==array_length_x/2)
			{
			  #ifdef BINARY
			  fprintf(fp_1d,"%ld\t%lf\t%lf\t%lf\n",j,dep_var[array_index].fi,dep_var[array_index].mu[0],c[0]);
			  #endif
			
			  #ifdef TERNARY
			  fprintf(fp_1d,"%ld\t%lf\t%lf\t%lf\t%lf\t%lf\n",j,dep_var[array_index].fi,dep_var[array_index].mu[0],dep_var[array_index].mu[1],c[0],c[1]);
			  #endif
	
			  #ifdef QUARTERNARY
			  fprintf(fp_1d,"%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",j,dep_var[array_index].fi,dep_var[array_index].mu[0],dep_var[array_index].mu[1],dep_var[array_index].mu[2],c[0],c[1],c[2]);
			  #endif
			}
			#endif
			


		}
		fprintf(fp,"\n");
	}

	fclose(fp);
	
	#ifdef ONE_D
	fclose(fp_1d);
	#endif
	

}
