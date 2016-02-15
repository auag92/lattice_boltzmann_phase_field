void comp_d_a_d_fi(struct grad_info *grad,double *d_a_d_fi)
{
	

	double mod_grad_fi[2],mod_grad_fi_raised_2_6[2],mod_grad_fi_sq[2],fi_x_4_pl_fi_y_4[2],mod_grad_fi_raised_2_4[2];
	
	double a_c[2];
	
	double fi_x_cubed,fi_y_cubed;

	

	//==========================================================================//
	//forming d_a_d_fi_x

	//==========================================================================//	
	mod_grad_fi[X]=sqrt(grad->grad_fi[X]*grad->grad_fi[X]+grad->grad_fi_central[Y]*grad->grad_fi_central[Y]);
	mod_grad_fi[Y]=sqrt(grad->grad_fi_central[X]*grad->grad_fi_central[X]+grad->grad_fi[Y]*grad->grad_fi[Y]);
	//=========================================================================//


	//=========================================================================//
	mod_grad_fi_sq[X]=mod_grad_fi[X]*mod_grad_fi[X];

	mod_grad_fi_sq[Y]=mod_grad_fi[Y]*mod_grad_fi[Y];	
	//=========================================================================//
	

	//=========================================================================//
	mod_grad_fi_raised_2_4[X]=mod_grad_fi_sq[X]*mod_grad_fi_sq[X];

	mod_grad_fi_raised_2_4[Y]=mod_grad_fi_sq[Y]*mod_grad_fi_sq[Y];
	//=========================================================================//
	


	//=========================================================================//
	mod_grad_fi_raised_2_6[X]=mod_grad_fi_sq[X]*mod_grad_fi_raised_2_4[X];

	mod_grad_fi_raised_2_6[Y]=mod_grad_fi_sq[Y]*mod_grad_fi_raised_2_4[Y];
	//=========================================================================//	

	//=========================================================================//
	fi_x_cubed=grad->grad_fi[X]*grad->grad_fi[X]*grad->grad_fi[X];

	fi_y_cubed=grad->grad_fi[Y]*grad->grad_fi[Y]*grad->grad_fi[Y];
	//=========================================================================//

	//=========================================================================//
	fi_x_4_pl_fi_y_4[X]=grad->grad_fi[X]*grad->grad_fi[X]*grad->grad_fi[X]*grad->grad_fi[X]+grad->grad_fi_central[Y]*grad->grad_fi_central[Y]*grad->grad_fi_central[Y]*grad->grad_fi_central[Y];	

	fi_x_4_pl_fi_y_4[Y]=grad->grad_fi_central[X]*grad->grad_fi_central[X]*grad->grad_fi_central[X]*grad->grad_fi_central[X]+grad->grad_fi[Y]*grad->grad_fi[Y]*grad->grad_fi[Y]*grad->grad_fi[Y];	
	//=========================================================================//

	//=========================================================================//
	a_c[X]=1.0-delta*(3.0-4.0*(fi_x_4_pl_fi_y_4[X]/mod_grad_fi_raised_2_4[X]));
			
	a_c[Y]=1.0-delta*(3.0-4.0*(fi_x_4_pl_fi_y_4[Y]/mod_grad_fi_raised_2_4[Y]));
	//=========================================================================//


	//==============================================================================//

	//forming d_a_d_fi_x

	if(fabs(grad->grad_fi[X])<1e-12 && fabs(grad->grad_fi_central[Y])<1e-12)
		d_a_d_fi[X]=0.0;			

	else
	{
		d_a_d_fi[X]=2.0*gamma_surf*epsilon*a_c[X]*a_c[X]*grad->grad_fi[X];

		d_a_d_fi[X]+=32.0*gamma_surf*epsilon*delta*mod_grad_fi_sq[X]*a_c[X]*((fi_x_cubed/mod_grad_fi_raised_2_4[X])-((fi_x_4_pl_fi_y_4[X]*grad->grad_fi[X])/mod_grad_fi_raised_2_6[X]));	
	}

	//forming d_a_d_fi_y
	
	if(fabs(grad->grad_fi_central[X])<1e-12 && fabs(grad->grad_fi[Y])<1e-12)
		d_a_d_fi[Y]=0.0;	

	else
	{
		d_a_d_fi[Y]=2.0*gamma_surf*epsilon*a_c[Y]*a_c[Y]*grad->grad_fi[Y];

		d_a_d_fi[Y]+=32.0*gamma_surf*epsilon*delta*mod_grad_fi_sq[Y]*a_c[Y]*((fi_y_cubed/mod_grad_fi_raised_2_4[Y])-((fi_x_4_pl_fi_y_4[Y]*grad->grad_fi[Y])/mod_grad_fi_raised_2_6[Y]));	
	}

	//===============================================================================//
}		
