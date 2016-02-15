void compini(struct var *dep_var)
{

	long i,j;

	long m;
	
	long array_index;
	
	FILE *fp;
	char fname[300];

	double dist_1,dist_2,dist_3,dist_4,dist_5,dist_6,dist_7,dist_8,dist_9,dist_10;

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
	
	//printf("in compini_multi\n");
	//printf("x1=%ld\ty1=%ld\tz1=%ld\n",x_1,y_1,z_1);
	//printf("x2=%ld\ty2=%ld\tz2=%ld\n",x_2,y_2,z_2);
	//printf("x3=%ld\ty3=%ld\tz3=%ld\n",x_3,y_3,z_3);
	//printf("x4=%ld\ty4=%ld\tz4=%ld\n",x_4,y_4,z_4);
	//printf("x5=%ld\ty5=%ld\tz5=%ld\n",x_5,y_5,z_5);	
	
	//====================================================================//
	//====================================================================//
	//I am going to consider the real domain 
	for(j=0;j<array_length_y;j++)
	{	
		for(i=0;i<array_length_x;i++)
		{
			array_index=i+array_length_x*j;

			//computing the distances
			dist_1=sqrt((i-x_1)*(i-x_1)*dx*dx+(j-y_1)*(j-y_1)*dy*dy);
			
			dist_2=sqrt((i-x_2)*(i-x_2)*dx*dx+(j-y_2)*(j-y_2)*dy*dy);

			dist_3=sqrt((i-x_3)*(i-x_3)*dx*dx+(j-y_3)*(j-y_3)*dy*dy);
	
			dist_4=sqrt((i-x_4)*(i-x_4)*dx*dx+(j-y_4)*(j-y_4)*dy*dy);
	
			dist_5=sqrt((i-x_5)*(i-x_5)*dx*dx+(j-y_5)*(j-y_5)*dy*dy);
	
			dist_6=sqrt((i-x_6)*(i-x_6)*dx*dx+(j-y_6)*(j-y_6)*dy*dy);
			
			dist_7=sqrt((i-x_7)*(i-x_7)*dx*dx+(j-y_7)*(j-y_7)*dy*dy);

			dist_8=sqrt((i-x_8)*(i-x_8)*dx*dx+(j-y_8)*(j-y_8)*dy*dy);
	
			dist_9=sqrt((i-x_9)*(i-x_9)*dx*dx+(j-y_9)*(j-y_9)*dy*dy);
	
			dist_10=sqrt((i-x_10)*(i-x_10)*dx*dx+(j-y_10)*(j-y_10)*dy*dy);
				
			if(dist_1<=rad_init || dist_2<=rad_init || dist_3<=rad_init || dist_4<=rad_init || dist_5<=rad_init || dist_6<=rad_init || dist_7<=rad_init || dist_8<=rad_init || dist_9<=rad_init || dist_10<=rad_init)
			{
				dep_var[array_index].fi=fi_s;

				mu_fr_c(mu_init,c_s_init,1);		

				for(m=0;m<NUMCOMPONENTS-1;m++)
				{
					dep_var[array_index].mu[m]=mu_init[m];
				}			

			}

			else
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
			
			fprintf(fp,"%ld\t%ld\t%lf\t%lf\t%lf\t%lf\n",i,j,dep_var[array_index].fi,dep_var[array_index].mu[0],dep_var[array_index].mu[1]);
			#endif 	

			#ifdef QUARTERNARY 
			
			fprintf(fp,"%ld\t%ld\t%lf\t%lf\t%lf\t%lf\n",i,j,dep_var[array_index].fi,dep_var[array_index].mu[0],dep_var[array_index].mu[1],dep_var[array_index].mu[2]);
			#endif 	

		}
		fprintf(fp,"\n");
	}
	
	//==================================================================//	

	fclose(fp);
}
