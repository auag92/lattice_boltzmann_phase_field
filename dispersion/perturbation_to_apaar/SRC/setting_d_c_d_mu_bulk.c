double set_d_c_d_mu_bulk(double *diff_pot,long phase_index,long comp_index_i,long comp_index_j)
{

	//=========================================================//
	//here we are going to code in the derivatives d_c_d_mu

	//they cab be functions of mu
	//=====================================//	
	#ifdef BINARY

	//for the solid
	d_c_d_mu_eq[1][0][0]=K[0];
	

	//for the liquid
	d_c_d_mu_eq[0][0][0]=1.0;
	
	#endif
	//====================================//


	//===================================//
	#ifdef TERNARY
	
	//for the solid	
	//d_c_d_mu_eq[1][0][0]=0.0489/0.2527;

	//d_c_d_mu_eq[1][1][1]=0.2527/0.2527;

	//d_c_d_mu_eq[1][0][1]=-0.0244/0.2527;

	//d_c_d_mu_eq[1][1][0]=-0.0244/0.2527;


	d_c_d_mu_eq[1][0][0]=0.0196;

	d_c_d_mu_eq[1][1][1]=0.0099;

	d_c_d_mu_eq[1][0][1]=-0.0098;

	d_c_d_mu_eq[1][1][0]=-0.0098;

	//for the liquid
	//d_c_d_mu_eq[0][0][0]=0.0489/0.2527;

	//d_c_d_mu_eq[0][1][1]=0.0216/0.2527;

	//d_c_d_mu_eq[0][0][1]=-0.0244/0.2527;

	//d_c_d_mu_eq[0][1][0]=-0.0244/0.2527;


		
	d_c_d_mu_eq[0][0][0]=0.09;

	d_c_d_mu_eq[0][1][1]=0.04749;

	d_c_d_mu_eq[0][0][1]=-0.045;

	d_c_d_mu_eq[0][1][0]=-0.045;
	
	#endif	

	//==============================//

	//===================================//
	#ifdef QUARTERNARY
	
	//for the solid	
	d_c_d_mu_eq[1][0][0]=0.000555;

	d_c_d_mu_eq[1][0][1]=-0.000013;
	
	d_c_d_mu_eq[1][0][2]=-0.000099;

	d_c_d_mu_eq[1][1][0]=-0.000013;	

	d_c_d_mu_eq[1][1][1]=0.023264;

	d_c_d_mu_eq[1][1][2]=-0.004229;

	d_c_d_mu_eq[1][2][0]=-0.000099;	

	d_c_d_mu_eq[1][2][1]=-0.004229;

	d_c_d_mu_eq[1][2][2]=0.145959;


	//for the liquid
	d_c_d_mu_eq[0][0][0]=0.003190;

	d_c_d_mu_eq[0][0][1]=-0.000029;
	
	d_c_d_mu_eq[0][0][2]=-0.000630;

	d_c_d_mu_eq[0][1][0]=-0.000029;	

	d_c_d_mu_eq[0][1][1]=0.009115;

	d_c_d_mu_eq[0][1][2]=-0.001811;

	d_c_d_mu_eq[0][2][0]=-0.000630;	

	d_c_d_mu_eq[0][2][1]=-0.001811;

	d_c_d_mu_eq[0][2][2]=0.158130;

	
	#endif	

	//==============================//
	
	//===================================================//

	return(d_c_d_mu_eq[phase_index][comp_index_i][comp_index_j]);

}
