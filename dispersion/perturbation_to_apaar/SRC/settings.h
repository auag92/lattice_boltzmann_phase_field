//some code settings which will determine the type of simulation to be run

//======================================================================================//
//this will have the information about the no. of componenets in my system
#define QUARTERNARY//this tag is important for using the right  matrix inversion formula 

#ifdef BINARY
#define NUMCOMPONENTS 2//this is the explicit no. of components
#endif

#ifdef TERNARY
#define NUMCOMPONENTS 3//this is the explicit no. of components
#endif

#ifdef QUARTERNARY
#define NUMCOMPONENTS 4//this is the explicit no. of components
#endif

//======================================================================================//

//======================================================================================//
#define ANTITRAPPING//with "ANTITRAPPING" I am going to turn on 2 things: a variable mobility, a non-zero anti-trapping current ; the options are "NO_ANTITRAPPING" and "ANTITRAPPING"
#define NO_SOLID_DIFFUSIVITY//this will specify the diffusivity in the solid ; the options are "SOLID_DIFFUSIVITY" and "NO_SOLID_DIFFFUSIVITY". This will be active only when there is NO_ANTITRAPPING  
#define NOFLUX//this is related to the boundary condition (the options can be "NOFLUX" or "PERIODIC")
#define ONE_D//this controls the initial conditions 
//=======================================================================================//

//=======================================================================================//
//setting whether checking mode is on
//#define CHECK_OFF//the other option is CHECK_OFF
//=======================================================================================//
 
//========================================================================================//
//here I am going to put in some code constants 

//========================================================================================//
#define BUFFER_WIDTH 2//the width of the buffer layer
//some constants related to the anti-trapping current 
#define M_const 0.018252
#define F_const 0.147715 
//========================================================================================//


//========================================================================================//
//initialising the directions to the numerical values
#define	X 0
#define Y 1
//========================================================================================//

//=======================================================================================//
//setting the no. of waves
#define no_of_waves 0.5
//=======================================================================================//

//so we are going to declare some structure types
//========================================================================================//
//declaring some other structure types
//a structure for the dependent varables 
struct var {
	double fi;//I am going consider a 2 - phase equilibrium

	double mu[NUMCOMPONENTS-1];//the no. of diffusion potentials I require	
};	


//a structure to hold the gradients in different fields (I am going to restrict to a 2D problem)
struct grad_info {
	double grad_fi[2];
	double grad_fi_central[2];
	double grad_mu[2][NUMCOMPONENTS-1];
	double M_mid[2][NUMCOMPONENTS-1][NUMCOMPONENTS-1];//the mobility	
	double j_at[2][NUMCOMPONENTS-1];
	double d_fi_d_t;  	
};




//=====================================================================================================//


//=====================================================================================================//
//allocating the matrices which are dependent on the no. of componenents
double D_S[NUMCOMPONENTS-1][NUMCOMPONENTS-1],D_L[NUMCOMPONENTS-1][NUMCOMPONENTS-1];//the diffusivity matrices

double c_s_min_c_l_eq[NUMCOMPONENTS-1];

double c_s_eq[NUMCOMPONENTS-1],c_l_eq[NUMCOMPONENTS-1];//they are going to contain the phase diagram information

double mu_eq[NUMCOMPONENTS-1];

double c_l_init[NUMCOMPONENTS-1],c_s_init[NUMCOMPONENTS-1];	

double d_c_d_mu_eq[2][NUMCOMPONENTS-1][NUMCOMPONENTS-1]; //the fist index will denote the phase information
//I am going to restrict myself to a 2 phase equilibrium

double K[NUMCOMPONENTS-1];

//======================================================================================================//
  	

//======================================================================================================//
//some system constants to be use by all routines
 
double gamma_surf,epsilon,delta,tau;
double fi_s,fi_l;

double dt,total_time;//I can also get rid of the 'total time' thing and run my simulations as long it is required to get an equilibrium profile   

long output_tsteps_interval,total_tsteps;

//long p_s_tsteps;
	
long nodes_x,nodes_y,array_length_x,array_length_y,tot_array_length;
	
double dx,dy,d_sp;//the last quantity will be used for reading the d-spacings    

double rad_init;//the initial radius of the solid nucleus at the centre/corner of the domain (this is what is needed for studyng dendrites)

double ratio; //the ratio of epsilon to dx (important for verification of anti-trapping) 

//int shift_factor,solid_phase_factor;//for directional simulations

long shift_length,solid_phase_length;//for directional simulations 

long eq_node_y;

long sine_wave_length;

long sine_wave_amplitude;

double eq_time;

long SEED_1,SEED_2;//the wave length of the initial perturbed interface

long noise_amplitude_1;//the amplitude of the initial perturbed interface

double noise_amplitude_2;

long start_time,start_shift_count,shift_count;

long x_1,x_2,x_3,x_4,x_5,x_6,x_7,x_8,x_9,x_10;

long y_1,y_2,y_3,y_4,y_5,y_6,y_7,y_8,y_9,y_10;

//======================================================================================================//
