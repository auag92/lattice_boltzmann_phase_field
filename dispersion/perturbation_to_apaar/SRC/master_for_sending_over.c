void cr_arr_for_send_over(struct var *temp_dep_var,struct domain_info dom_decomp_info,struct offset_info offset, struct var *dep_var,long worker_number)
{
	int i,j;
	long temp_array_index,array_index;

	//char fname[300];
	//FILE *fp_ch;

//	sprintf(fname,"../DATA/master_ph_inter_profile_taskid_%d.dat",worker_number);
//	if((fp_ch=fopen(fname,"w"))==NULL)
//	{
//		printf("check profile can't be written\n");
//		exit(1);
//	}

	for(j=0;j<dom_decomp_info.rows;j++)
	{	
		for(i=0;i<dom_decomp_info.cols;i++)
		{
			temp_array_index = i + dom_decomp_info.cols * j ;
			
			array_index = (offset.x + i) + nodes_x* (offset.y + j) ;

			temp_dep_var[temp_array_index] = dep_var[array_index];

			//fprintf(fp_ch,"%d\t%d\t%d\t%lf\n",i,j,k,temp_ph_inter[temp_array_index].comp);

				
		}
		//fprintf(fp_ch, "\n");
	}
	
	//fclose(fp_ch);	
} 							
