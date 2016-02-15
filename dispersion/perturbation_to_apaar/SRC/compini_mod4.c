void compini(struct var *dep_var)
{
  
	#ifndef QUARTERNARY
	printf("equilibrium profile can't be read\n")
	exit(1);
	#endif

	long i,j;
	
	long j_read;

	long m;
	
	long array_index;
	
	FILE *fp,*fp_bound,*fp_read,*fp_check,*fp_1d;
	char fname[300];

	
	//double mu_init[NUMCOMPONENTS-1];

		
	struct read_var {
	  double fi;
	  double mu[NUMCOMPONENTS-1];
	  double c[NUMCOMPONENTS-1];
	} *read_dep_var;
	
	
	long read_int_left,read_int_right;//this for storing the points just to the left and right of fi=0.5
	double fi_int_left,fi_int_right;//this is for storing the values of 'fi' at those points
	
	
	struct dist {
	  long cr_int_left;
	  long cr_int_right;
	  double cr_int;
	} *dist_local;
		
	long count;
	
	
	
	//====================================================================//
	//opening the file for reading in the equilibrium profiles
	sprintf(fname,"../INPAR/one_D_phase_field_profiles_%.5lf.dat",eq_time);
	if((fp_read=fopen(fname,"r"))==NULL)
	{
		printf("Equilibrium profile can't be read\n");
		exit(1);
	}
	
	sprintf(fname,"../DATA/sup_sat_%.5lf_%.5lf_%.5lf_w_l_%ld/check.dat",c_l_eq[0]-c_l_init[0],c_l_eq[1]-c_l_init[1],c_l_eq[2]-c_l_init[2],sine_wave_length);
	if((fp_check=fopen(fname,"w"))==NULL)
	{
		printf("check profile can't be written\n");
		exit(1);
	}
	
	sprintf(fname,"../DATA/sup_sat_%.5lf_%.5lf_%.5lf_w_l_%ld/initial_profile.dat",c_l_eq[0]-c_l_init[0],c_l_eq[1]-c_l_init[1],c_l_eq[2]-c_l_init[2],sine_wave_length);
	if((fp=fopen(fname,"w"))==NULL)
	{
		printf("Initial profile can't be written\n");
		exit(1);
	}

	sprintf(fname,"../DATA/sup_sat_%.5lf_%.5lf_%.5lf_w_l_%ld/initial_boundary_profile.dat",c_l_eq[0]-c_l_init[0],c_l_eq[1]-c_l_init[1],c_l_eq[2]-c_l_init[2],sine_wave_length);
	if((fp_bound=fopen(fname,"w"))==NULL)
	{
		printf("Initial boundary profile can't be written\n");
		exit(1);
	}
	
	sprintf(fname,"../DATA/sup_sat_%.5lf_%.5lf_%.5lf_w_l_%ld/one_D_initial_profile.dat",c_l_eq[0]-c_l_init[0],c_l_eq[1]-c_l_init[1],c_l_eq[2]-c_l_init[2],sine_wave_length);
	if((fp_1d=fopen(fname,"w"))==NULL)
	{
		printf("1d initial profile can't be written\n");
		exit(1);
	}
	//====================================================================//
	
	//=====================================================================//
	//allocating space for the read_dep var array
	
	if((read_dep_var=(struct read_var *)malloc(eq_node_y*sizeof(struct read_var)))==NULL)
	{
		printf("Space creation failed for the dep_var array\n");
		exit(1);
	}
	
	
	//============================================================================================================================//
	//reading the equilibrium profiles and computing the maximum/minimum in composition and storing their location
	for(j=0;j<eq_node_y;j++)
	{
	  fscanf(fp_read,"%ld%lf%lf%lf%lf%lf%lf%lf",&j_read,&(read_dep_var[j].fi),&(read_dep_var[j].mu[0]),&(read_dep_var[j].mu[1]),&(read_dep_var[j].mu[2]),&(read_dep_var[j].c[0]),&(read_dep_var[j].c[1]),&(read_dep_var[j].c[2]));
	  
	  fprintf(fp_check,"%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",j_read,read_dep_var[j].fi,read_dep_var[j].mu[0],read_dep_var[j].mu[1], 
		  read_dep_var[j].mu[2],read_dep_var[j].c[0],read_dep_var[j].c[1],read_dep_var[j].c[2]);
	  
	  //extracting the point corresponfing to fi=0.5
	  if(j>0)
	  {
	    if(read_dep_var[j-1].fi>=0.5 && read_dep_var[j].fi<0.5)
	    {
	      read_int_left=j-1;
	      read_int_right=j;
	      
	      fi_int_left=read_dep_var[j-1].fi;
	      fi_int_right=read_dep_var[j].fi;
	    }
	  }
	  
	  	  
	}
	
	fclose(fp_read);
	fclose(fp_check);
	
	
	printf("read:loc_left=%ld val_left=%lf\n",read_int_left,fi_int_left);
	printf("read:loc_right=%ld val_right=%lf\n",read_int_right,fi_int_right);
	//=============================================================================================================================//
	
	//=============================================================================================================================//
	//allocating space for the dist_local array
	if((dist_local=(struct dist *)malloc(array_length_x*sizeof(struct dist)))==NULL)
	{
		printf("Space creation failed for the dist_local array\n");
		exit(1);
	}
	//==============================================================================================================================//
	
	
	//==============================================================================================================================//
	//creating the sinusoidal perturbation at the interface
	for(i=0;i<array_length_x;i++)
	{
	  dist_local[i].cr_int_left=(solid_phase_length+BUFFER_WIDTH)+(long)(floor(sine_wave_amplitude*sin((M_PI/2.0)+(2.0*M_PI*no_of_waves*(i-BUFFER_WIDTH)/nodes_x))));
		
	  dist_local[i].cr_int_right=(solid_phase_length+BUFFER_WIDTH)+(long)(ceil(sine_wave_amplitude*sin((M_PI/2.0)+(2.0*M_PI*no_of_waves*(i-BUFFER_WIDTH)/nodes_x))));
	  
	  if(dist_local[i].cr_int_left==dist_local[i].cr_int_right)
	    dist_local[i].cr_int_right++;
		
	  dist_local[i].cr_int=(solid_phase_length+BUFFER_WIDTH)+(sine_wave_amplitude*sin((M_PI/2.0)+(2.0*M_PI*no_of_waves*(i-BUFFER_WIDTH)/nodes_x)));
		
	  fprintf(fp_bound,"%ld\t%ld\t%ld\n",i,dist_local[i].cr_int_left,dist_local[i].cr_int_right);
	}		

	fclose(fp_bound);
	//==============================================================================================================================//

	//==============================================================================================================================//
	//I am going to consider the real domain 
	//the order of the loops have been interchanged	
	for(i=0;i<array_length_x;i++)
	{	
	  
		//creating the liquid
		count=read_int_right;
		  
		for(j=dist_local[i].cr_int_right;j<array_length_y;j++)//initialising the liquid part  
		{
		    array_index=i+array_length_x*j;

		    dep_var[array_index].fi=read_dep_var[count].fi;
			
		    for(m=0;m<NUMCOMPONENTS-1;m++)
		    {
		      dep_var[array_index].mu[m]=read_dep_var[count].mu[m];
		    }
			
		    if(count<eq_node_y-1)
		      count++;

		}
		
		//creating the solid		  
		count=read_int_left;
		
		for(j=dist_local[i].cr_int_left;j>=0;j--)  
		{
			array_index=i+array_length_x*j;

			dep_var[array_index].fi=read_dep_var[count].fi;

			//mu_fr_c(mu_init,c_s_init,1);					

			for(m=0;m<NUMCOMPONENTS-1;m++)
			{
				dep_var[array_index].mu[m]=read_dep_var[count].mu[m];
			}
			
			if(count>0)
			count--;
		}
	    }
	    
	   
	//==================================================================//	

	for(j=0;j<array_length_y;j++)
	{
		for(i=0;i<array_length_x;i++)
		{
			array_index=i+array_length_x*j;		

			#ifdef BINARY 
			
			fprintf(fp,"%ld\t%ld\t%lf\t%lf\n",i,j,dep_var[array_index].fi,dep_var[array_index].mu[0]);
			#endif 	
	
			#ifdef TERNARY 
			
			fprintf(fp,"%ld\t%ld\t%lf\t%lf\t%lf\n",i,j,dep_var[array_index].fi,dep_var[array_index].mu[0],dep_var[array_index].mu[1]);
			#endif 	

			#ifdef QUARTERNARY 
			
			fprintf(fp,"%ld\t%ld\t%lf\t%lf\t%lf\t%lf\n",i,j,dep_var[array_index].fi,dep_var[array_index].mu[0],dep_var[array_index].mu[1],dep_var[array_index].mu[2]);
			#endif 	
			
			if(i==array_length_x/2)
			  fprintf(fp_1d,"%ld\t%lf\t%lf\t%lf\t%lf\n",j,dep_var[array_index].fi,dep_var[array_index].mu[0],dep_var[array_index].mu[1],dep_var[array_index].mu[2]);
			
		}
		fprintf(fp,"\n");
	}

	fclose(fp);
	fclose(fp_1d);

	free(dist_local);
	free(read_dep_var);
}
