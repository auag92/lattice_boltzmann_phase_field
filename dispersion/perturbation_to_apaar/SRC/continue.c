void cont(struct var *dep_var,long time)
{
	#ifndef QUARTERNARY
	printf("can't continue without quarternary\n");
	exit(1);
	#endif	
		
	char fname[300];
        long i,j;
        long i_read,j_read;
	
        double t=(double)time;

        long array_index;

	double c[NUMCOMPONENTS-1];	

        FILE *fp,*fp_read;

        sprintf(fname,"../INPAR/RESTART/phase_field_profiles_%.5lf.dat",t);
        if((fp_read=fopen(fname,"r"))==NULL)
        {
                printf("profile can't be opened for restarting simulation\n");
                exit(1);
        }

        //opening files for writing the array
        sprintf(fname,"../DATA/restart_phase_field_profiles_%.5lf.dat",t);
        if((fp=fopen(fname,"w"))==NULL)
        {
                printf("Initial profile can't be written\n");
                exit(1);
        }

	for(j=0;j<array_length_y;j++)
	{	
		for(i=0;i<array_length_x;i++)
		{
			array_index=i+array_length_x*j;


			//==========================================================//		
			//the quarternary thing is hard coded
                       	fscanf(fp_read,"%ld%ld%lf%lf%lf%lf%lf%lf%lf",&i_read,&j_read,&(dep_var[array_index].fi),&(dep_var[array_index].mu[0]),&(dep_var[array_index].mu[1]),&(dep_var[array_index].mu[2]),&(c[0]),&(c[1]),&(c[2]));

			fprintf(fp,"%ld\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",i,j,dep_var[array_index].fi,dep_var[array_index].mu[0],dep_var[array_index].mu[1],dep_var[array_index].mu[2],c[0],c[1],c[2]);
			//==========================================================//

                        		

		}
		fprintf(fp,"\n");
			
        }
                
        fclose(fp_read);
        fclose(fp);
}