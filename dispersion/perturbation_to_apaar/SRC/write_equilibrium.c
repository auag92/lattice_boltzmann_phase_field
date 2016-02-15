void write_eq()
{
	char fname[1000];

	FILE *fp_eq;

	sprintf(fname,"../DATA/sup_sat_%.5lf_%.5lf_%.5lf_w_l_%ld/equilibrium_settings.dat",c_l_eq[0]-c_l_init[0],c_l_eq[1]-c_l_init[1],c_l_eq[2]-c_l_init[2],sine_wave_length);
	if((fp_eq=fopen(fname,"w"))==NULL)
	{
		printf("equilibrium settings can't be written\n");
		exit(1);
	}	

	//=============================================//
	#ifdef BINARY
	
	//writing the files
	fprintf(fp_eq,"BINARY\n");
	fprintf(fp_eq,"mu_eq=%lf \n",mu_eq[0]);		
	fprintf(fp_eq,"c_s_eq=%lf c_l_eq=%lf c_s_min_c_l_eq=%lf\n",c_s_eq[0],c_l_eq[0],c_s_min_c_l_eq[0]);
	fprintf(fp_eq,"c_s_init=%lf c_l_init=%lf\n",c_s_init[0],c_l_init[0]);	
	#endif
	//=============================================//


	//=============================================//
	#ifdef TERNARY
	
	//writing the files
	fprintf(fp_eq,"TERNARY\n");
	fprintf(fp_eq,"mu_eq_1=%lf\n",mu_eq[0]);	
	fprintf(fp_eq,"mu_eq_2=%lf\n",mu_eq[1]);	
	fprintf(fp_eq,"c_s_eq_1=%lf c_l_eq_1=%lf c_s_min_c_l_eq_1=%lf\n",c_s_eq[0],c_l_eq[0],c_s_min_c_l_eq[0]);	
	fprintf(fp_eq,"c_s_eq_2=%lf c_l_eq_2=%lf c_s_min_c_l_eq_2=%lf\n",c_s_eq[1],c_l_eq[1],c_s_min_c_l_eq[1]);	
	fprintf(fp_eq,"c_s_init_1=%lf c_l_init_1=%lf\n",c_s_init[0],c_l_init[0]);	
	fprintf(fp_eq,"c_s_init_2=%lf c_l_init_2=%lf\n",c_s_init[1],c_l_init[1]);	
	
	#endif
	//===========================================//

	//=============================================//
	#ifdef QUARTERNARY
	
	//writing the files
	fprintf(fp_eq,"QUARTERNARY\n");
	fprintf(fp_eq,"mu_eq_1=%lf\n",mu_eq[0]);	
	fprintf(fp_eq,"mu_eq_2=%lf\n",mu_eq[1]);	
	fprintf(fp_eq,"mu_eq_3=%lf\n",mu_eq[2]);	
	fprintf(fp_eq,"c_s_eq_1=%lf c_l_eq_1=%lf c_s_min_c_l_eq_1=%lf\n",c_s_eq[0],c_l_eq[0],c_s_min_c_l_eq[0]);	
	fprintf(fp_eq,"c_s_eq_2=%lf c_l_eq_2=%lf c_s_min_c_l_eq_2=%lf\n",c_s_eq[1],c_l_eq[1],c_s_min_c_l_eq[1]);	
	fprintf(fp_eq,"c_s_eq_3=%lf c_l_eq_3=%lf c_s_min_c_l_eq_3=%lf\n",c_s_eq[2],c_l_eq[2],c_s_min_c_l_eq[2]);	
	fprintf(fp_eq,"c_s_init_1=%lf c_l_init_1=%lf\n",c_s_init[0],c_l_init[0]);	
	fprintf(fp_eq,"c_s_init_2=%lf c_l_init_2=%lf\n",c_s_init[1],c_l_init[1]);	
	fprintf(fp_eq,"c_s_init_3=%lf c_l_init_3=%lf\n",c_s_init[2],c_l_init[2]);
	fprintf(fp_eq,"D_11=%lf D_12=%lf D_13=%lf D_21=%lf D_22=%lf D_23=%lf D_31=%lf D_32=%lf D_33=%lf\n",D_L[0][0],D_L[0][1],D_L[0][2],D_L[1][0],D_L[1][1],D_L[1][2],D_L[2][0],D_L[2][1],D_L[2][2]);
	
	#endif
	//===========================================//
	
	fclose(fp_eq);
}
