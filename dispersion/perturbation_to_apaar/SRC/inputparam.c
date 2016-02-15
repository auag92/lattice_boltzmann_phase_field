//this is a routine to read the system parameters

void input(long taskid)
{
	char str[100],fname[1000],dirname[1000];
	FILE *fpread,*fpwrite;
	//opening the files for the reading and writing of system parameters
	if((fpread=fopen("../INPAR/read.dat","r"))==NULL)
	{
		printf("System parameters can't be read\n");
		exit(1);
	}
	//reading from one and writing in the other
	//fscanf(fpread,"%s%ld",str,&nodes_x);
	fscanf(fpread,"%s%ld",str,&nodes_y);
	fscanf(fpread,"%s%lf",str,&d_sp);
	fscanf(fpread,"%s%lf",str,&fi_s);
	fscanf(fpread,"%s%lf",str,&fi_l);
	fscanf(fpread,"%s%lf",str,&gamma_surf);
	fscanf(fpread,"%s%lf",str,&ratio);
	fscanf(fpread,"%s%lf",str,&delta);	
	fscanf(fpread,"%s%lf",str,&dt);
	fscanf(fpread,"%s%lf",str,&total_time);
	//fscanf(fpread,"%s%ld",str,&p_s_tsteps);	
	fscanf(fpread,"%s%ld",str,&output_tsteps_interval);
	fscanf(fpread,"%s%lf",str,&rad_init);
	fscanf(fpread,"%s%ld",str,&shift_length);
	fscanf(fpread,"%s%ld",str,&solid_phase_length);	
	//fscanf(fpread,"%s%ld",str,&SEED_1);
	fscanf(fpread,"%s%ld",str,&SEED_2);	
	//fscanf(fpread,"%s%ld",str,&noise_amplitude_1);		
	fscanf(fpread,"%s%lf",str,&noise_amplitude_2);			
	fscanf(fpread,"%s%ld",str,&numworkers.x);
	fscanf(fpread,"%s%ld",str,&numworkers.y);
	//fscanf(fpread,"%s%ld",str,&start_time);
	//fscanf(fpread,"%s%ld",str,&start_shift_count);
	fscanf(fpread,"%s%ld",str,&sine_wave_amplitude);
	fscanf(fpread,"%s%lf",str,&eq_time);
	fscanf(fpread,"%s%ld",str,&eq_node_y);	
	
	fclose(fpread);

	if(taskid==numworkers.x*numworkers.y)//restricting the echo_back to the master
	{

	  sprintf(dirname,"../DATA/sup_sat_%.5lf_%.5lf_%.5lf_w_l_%ld",c_l_eq[0]-c_l_init[0],c_l_eq[1]-c_l_init[1],c_l_eq[2]-c_l_init[2],sine_wave_length);
	  if (mkdir(dirname, S_IRWXU) != 0)
	  {
		printf("error in creating directory\n");
		exit(1);
	  }


	  sprintf(fname,"../DATA/sup_sat_%.5lf_%.5lf_%.5lf_w_l_%ld/param.dat",c_l_eq[0]-c_l_init[0],c_l_eq[1]-c_l_init[1],c_l_eq[2]-c_l_init[2],sine_wave_length);
	  if((fpwrite=fopen(fname,"w"))==NULL)
	  {
		printf("System parameters can't be written\n");
		exit(1);
	  }

	
	  //fprintf(fpwrite,"nodes_x\t%ld\n",nodes_x);
	  fprintf(fpwrite,"nodes_y\t%ld\n",nodes_y);
	  fprintf(fpwrite,"d_sp\t%lf\n",d_sp);
	  fprintf(fpwrite,"fi_s\t%lf\n",fi_s);
	  fprintf(fpwrite,"fi_l\t%.12lf\n",fi_l);
	  fprintf(fpwrite,"gamma\t%lf\n",gamma_surf);
	  fprintf(fpwrite,"ratio\t%lf\n",ratio);
	  fprintf(fpwrite,"delta\t%lf\n",delta);
	  fprintf(fpwrite,"dt\t%lf\n",dt);
	  fprintf(fpwrite,"total_time\t%lf\n",total_time);
	  //fprintf(fpwrite,"p_s_tsteps\t%ld\n",p_s_tsteps);
	  fprintf(fpwrite,"output_time_interval\t%ld\n",output_tsteps_interval);
	  fprintf(fpwrite,"rad_init\t%lf\n",rad_init);
	  fprintf(fpwrite,"shift_length\t%ld\n",shift_length);
	  fprintf(fpwrite,"solid_phase_length\t%ld\n",solid_phase_length);	
	  //fprintf(fpwrite,"SEED_1\t%ld\n",SEED_1);
	  fprintf(fpwrite,"SEED_2\t%ld\n",SEED_2);
	  //fprintf(fpwrite,"noise_amplitude_1\t%ld\n",noise_amplitude_1);
	  fprintf(fpwrite,"noise_amplitude_2\t%lf\n",noise_amplitude_2);		
	  fprintf(fpwrite,"numworkers_x\t%ld\n",numworkers.x);
	  fprintf(fpwrite,"numworkers_y\t%ld\n",numworkers.y);
	  //fprintf(fpwrite,"start_time\t%ld\n",start_time);
	  //fprintf(fpwrite,"start_shift_count\t%ld\n",start_shift_count);
	  fprintf(fpwrite,"sine_wave_amplitude\t%ld\n",sine_wave_amplitude);
	  fprintf(fpwrite,"eq_time\t%lf\n",eq_time);
	  fprintf(fpwrite,"eq_node_y\t%ld\n",eq_node_y);	    
	
		fclose(fpwrite);
	}
}				
