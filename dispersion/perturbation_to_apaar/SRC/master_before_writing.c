void cr_arr_after_recv_from_worker(struct var *dep_var,struct var *temp_dep_var,struct offset_info offset,struct domain_info dom_decomp_info,int source)
{
	int i,j;
	
	long array_index, temp_array_index;

	//char fname[300];
	//FILE *fp_ch;

	//sprintf(fname,"../DATA/before_writing_master_ph_inter_profile_taskid_%d.dat",source);
	//if((fp_ch=fopen(fname,"w"))==NULL)
	//{
	//	printf("check profile can't be written\n");
	//	exit(1);
	//}


	for(j=0;j<dom_decomp_info.rows;j++)
	{
		for(i=0;i<dom_decomp_info.cols;i++)
		{
			temp_array_index=i+dom_decomp_info.cols*j;

			array_index = (offset.x	+ i) + nodes_x * (offset.y + j) ;
	
			dep_var[array_index] = temp_dep_var[temp_array_index];
			
			//fprintf(fp_ch,"%d\t%d\t%d\t%lf\t%lf\n",i,j,k,ph_inter[array_index].comp,temp_ph_inter[temp_array_index].comp);

				
		}
		//fprintf(fp_ch, "\n");
	}
	
	//fclose(fp_ch);	
}	
