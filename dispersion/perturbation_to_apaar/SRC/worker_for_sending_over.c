void cr_arr_for_send_over_fr_wrkr(struct var *temp_dep_var,struct domain_info dom_decomp_info,struct domain_info worker_dom_decomp_info,struct var *dep_var,int taskid)
{
	int i,j;
	long temp_array_index,array_index;

	//char fname[300];
	//FILE *fp_ch;

	//sprintf(fname,"../DATA/from_worker_send_ph_inter_profile_taskid_%d.dat",taskid);
	//if((fp_ch=fopen(fname,"w"))==NULL)
	//{
	//	printf("check profile can't be written\n");
	//	exit(1);
	//}

	for(j=BUFFER_WIDTH;j<(BUFFER_WIDTH+dom_decomp_info.rows);j++)
	{	
		for(i=BUFFER_WIDTH;i<(BUFFER_WIDTH+dom_decomp_info.cols);i++)
		{
			temp_array_index = (i-BUFFER_WIDTH) + dom_decomp_info.cols * (j-BUFFER_WIDTH) ;
			
			array_index = i + worker_dom_decomp_info.cols * j ;

			temp_dep_var[temp_array_index] = dep_var[array_index];

		}
		//fprintf(fp_ch, "\n");
	}
	
	//fclose(fp_ch);
} 			
