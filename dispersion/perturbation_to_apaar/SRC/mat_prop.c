void init_mat_prop()
{
//========================================================//
//========================================================//
//setting the diffusivity tensor 

//============================//
#ifdef BINARY

//for the solid
#ifdef ANTITRAPPING
D_S[0][0]=0.0;	
#endif

#ifdef NO_ANTITRAPPING

	#ifdef SOLID_DIFFUSIVITY
	D_S[0][0]=1.0;
	#endif	

	#ifdef NO_SOLID_DIFFUSIVITY
	D_S[0][0]=0.0;
	#endif	
		
#endif

//for the liquid
D_L[0][0]=1.0;
#endif
//=============================//	


//=============================//
#ifdef TERNARY
//for the solid
#ifdef ANTITRAPPING
D_S[0][0]=0.0;
D_S[0][1]=0.0;
D_S[1][0]=0.0;
D_S[1][1]=0.0;	
#endif

#ifdef NO_ANTITRAPPING
	#ifdef SOLID_DIFFUSIVITY
	D_S[0][0]=1.0;
	D_S[0][1]=0.0;
	D_S[1][0]=0.0;
	D_S[1][1]=1.0;
	#endif

	#ifdef NO_SOLID_DIFFUSIVITY
	D_S[0][0]=0.0;
	D_S[0][1]=0.0;
	D_S[1][0]=0.0;
	D_S[1][1]=0.0;
	#endif

			
#endif

//for the liquid
D_L[0][0]=1.0;
D_L[0][1]=0.0;
D_L[1][0]=0.0;
D_L[1][1]=1.0;	
#endif		
//============================//

//=============================//
#ifdef QUARTERNARY
//for the solid
#ifdef ANTITRAPPING
D_S[0][0]=0.0;
D_S[0][1]=0.0;
D_S[0][2]=0.0;
D_S[1][0]=0.0;
D_S[1][1]=0.0;	
D_S[1][2]=0.0;
D_S[2][0]=0.0;
D_S[2][1]=0.0;	
D_S[2][2]=0.0;
#endif

#ifdef NO_ANTITRAPPING
	#ifdef SOLID_DIFFUSIVITY
	D_S[0][0]=1.0;
	D_S[0][1]=0.0;
	D_S[0][2]=0.0;
	D_S[1][0]=0.0;
	D_S[1][1]=1.0;
	D_S[1][2]=0.0;
	D_S[2][0]=0.0;
	D_S[2][1]=0.0;
	D_S[2][2]=1.0;
	#endif

	#ifdef NO_SOLID_DIFFUSIVITY
	D_S[0][0]=0.0;
	D_S[0][1]=0.0;
	D_S[0][2]=0.0;
	D_S[1][0]=0.0;
	D_S[1][1]=0.0;
	D_S[1][2]=0.0;
	D_S[2][0]=0.0;
	D_S[2][1]=0.0;
	D_S[2][2]=0.0;
	#endif

			
#endif

//for the liquid
D_L[0][0]=0.5;
D_L[0][1]=0.0;
D_L[0][2]=0.0;	
D_L[1][0]=0.0;
D_L[1][1]=0.8;
D_L[1][2]=0.0;
D_L[2][0]=0.0;
D_L[2][1]=0.0;
D_L[2][2]=1.0;				
#endif		
//============================//


//=========================================================//
//=========================================================//

//=========================================================//
//=========================================================//
//the equilibrium values

#ifdef BINARY
//the information is in terms of B

//the mu's
//the equilibrium values 
mu_eq[0]=1.0;
		
//initialising the partition coefficient
K[0]=0.5;

//the compositions

//for the solid
c_s_eq[0]=K[0]*mu_eq[0]; 

//for the liquid	
c_l_eq[0]=mu_eq[0];

//creating the vector c_s_min_c_l_eq
c_s_min_c_l_eq[0]=c_s_eq[0]-c_l_eq[0];

//initialising the compositions

c_l_init[0]=0.7;

c_s_init[0]=0.5;   

#endif
//===============================================//

//===============================================//

#ifdef TERNARY 

//the information appearing below is in terms of A and B
//the mu's
//the equilibrium values 
mu_eq[0]=1.0;

mu_eq[1]=1.0;
			
//the compositions

//for the solid
c_s_eq[0]= 0.98;

c_s_eq[1]= 0.01;
	
//for the liquid	
c_l_eq[0]=0.9;
	
c_l_eq[1]=0.05;	
	
//creating the row_vector
c_s_min_c_l_eq[0]=c_s_eq[0]-c_l_eq[0];
		
c_s_min_c_l_eq[1]=c_s_eq[1]-c_l_eq[1];


//initialising the partition coefficient
K[0]=c_s_eq[0]/c_l_eq[0];

K[1]=c_s_eq[1]/c_l_eq[1];


//initialising the compositions in the system
c_l_init[0]=0.95;

c_l_init[1]=0.03;	

c_s_init[0]=0.98;

c_s_init[1]= 0.01;	



#endif 
//=============================================//

//===============================================//

#ifdef QUARTERNARY 

//the information appearing below is in terms of A, B and C
//the mu's
//the equilibrium values 
mu_eq[0]=1.0;

mu_eq[1]=1.0;

mu_eq[2]=1.0;
			
//the compositions

//for the solid
c_s_eq[0]=0.000556;

c_s_eq[1]=0.023831;

c_s_eq[2]=0.177447;
	
//for the liquid	
c_l_eq[0]=0.0032;
	
c_l_eq[1]=0.0092; 	

c_l_eq[2]=0.1969; 	
	
//creating the row_vector
c_s_min_c_l_eq[0]=c_s_eq[0]-c_l_eq[0];
		
c_s_min_c_l_eq[1]=c_s_eq[1]-c_l_eq[1];

c_s_min_c_l_eq[2]=c_s_eq[2]-c_l_eq[2];


//initialising the partition coefficient
K[0]=c_s_eq[0]/c_l_eq[0];

K[1]=c_s_eq[1]/c_l_eq[1];

K[2]=c_s_eq[2]/c_l_eq[2];


//initialising the compositions in the system
c_l_init[0]=0.001878;

c_l_init[1]=0.0165155;	

c_l_init[2]=0.1871735;	

c_s_init[0]=0.000556;

c_s_init[1]=0.023831;	

c_s_init[2]=0.177447;	


#endif 
//=============================================//

//========================================================//
}



