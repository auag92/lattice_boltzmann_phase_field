//a finite difference code for the 2D grand potential model
//I am going to use a 9-point stencil in my bid to implement anisotropy
//I need a buffer width which is 1 layer thick to implement the 9-point stencil
//I need to incorporate the anti-trapping current and put the interface under diffusive control

///=================================================================================================//
//I am going to make this code general enough to tackle a multicomponent situation starting from binary to higher up
//==================================================================================================// 
 

//including the header files
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include <sys/stat.h>


//inlcuding files 
//============================================================//
//the files which specify the constants

#include"settings.h"//it will determine the type of simulation to be run and initialise the variables

#include"const_mpi.h"//it will contain the constants related to "mpi" 

#include"mat_prop.c"//initialising the material properties
//==================================================================================//

#include"matrix_operations.h"//it will contain different matrix operation routines

#include"functions.h"//contains the routines to compute different interpolants

//the other routines
#include"inputparam.c"//to read some other input parameters 

#include"write_equilibrium.c"//to write back the equilibrium settings for confirmation

#include"setting_d_c_d_mu_bulk.c"

#include"mu_of_c.c"

#include"ran2.c"

//==================================================================//
//for the creation of initial composition 
#ifdef MULTI
  #include"input_multi.c"

  #include"compini_multi.c"//creating the multi particles
#endif

#ifdef ONE_D
  #include"compini_mod4.c"
#endif
//==================================================================//



#include"continue.c"


//===========================================//
//some routines related to the MPI 

#include"build_datatypes1.c"

#include"build_datatypes2.c"

#include"master_for_sending_over.c"

#include"master_before_writing.c"

#include"worker_transfer_from_temp.c"

#include"worker_for_sending_over.c"

#include"cr_temp_dep_var_for_left_right_sends.c"
#include"cr_temp_dep_var_for_left_right_recvs.c"
#include"neighbour_exchange_dep_var_left_right.c"

#include"cr_temp_dep_var_for_up_below_sends.c"
#include"cr_temp_dep_var_for_up_below_recvs.c"
#include"neighbour_exchange_dep_var_up_below.c"

#include"vert_exchange_dep_var.c"

#include"exchange_dep_var.c"
//=========================================//

#include"compute_tau.c"

#include"form_mobility.c"

#include"compute_gradients.c"

#include"compute_d_a_d_fi.c"

#include"divergence_fi.c" 

#include"c_of_mu.c"

#include"compute_anti_trapping_currents.c"

#include"update_fi_mu.c"

#include"check_for_shifting.c"

#include"comp_shift_profiles.c"

#include"d_A.c"

#include"write_output.c"

//==========================================//


main(int argc, char *argv[])
{
	//========COMMON REGION FOR ALL THE PROCESSES BEGINS=============================================//
	//===============================================================================================//

	//================================================================================================//
	//=====the arrays and variables which are required on all of the nodes============================//
	struct var *dep_var;//I am going to solve for 'mu' and 'fi'; I will get back 'c' from 'mu'
	
	struct grad_info *info;  

	double time;	

	long tsteps;

	long i,j;
	
	int shift_flag,tot_shift_flag;

	//============================================================================================//

	


	//===========================================================================================//	
	//the quantities relevant to MPI are
	int taskid; /* this task's unique id */
	int numtasks; /* number of tasks */
  	int destination, source; /* to - from for message send-receive */
 	int rc;// variable for storing the error code 
	
	//==========================================================================================//	

	//==========================================================================================//
	//some variables related to 3D domain decomposition 	
	struct domain_info dom_decomp_info,min_number,extra,worker_dom_decomp_info;
	
	struct neighbour_info neighbour;

	struct offset_info offset;
	
	struct worker_info worker;
	//=========================================================================================//
	
	//========================================================================================//
	//declaring an array for distributing work to the workers and also for receiving
	struct var *temp_dep_var;	  
	//========================================================================================//
		
	//===========================================================================================//	
	//===========================================================================================//	
	/* First, find out my taskid and how many tasks are running */
 	rc = MPI_Init(&argc,&argv);
  	rc|= MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  	rc|= MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

	
	if (rc != 0)
    		printf ("error initializing MPI and obtaining task ID information\n");

	//===========================================================================================//
	//============================================================================================//


	//====================================================================================//
	//initialising the material property settings (see "mat_prop.c")
	init_mat_prop();
	//===================================================================================//	

	//===================================================================================//
	//obtaining the sine wave length

	sine_wave_length=atol(argv[1]);
	//===================================================================================//	
					
	
	//===================================================================================//
	//reading the system parameters (see "inputparam.c")
	input(taskid);
	//===================================================================================//

	//===================================================================================//
	//===================================================================================//
	//computing some of the MPI rekated quantities 

	//computing the total no. of workers
	numworkers.tot=numworkers.x*numworkers.y;

	/* one of the processors is completely reserved for coordination */
	//we are going to call it the master
  	MASTER=numworkers.tot; //this is because I am going to allocate workers from 0
	
	if(taskid==MASTER)
		printf("The task id of the master is=%d\n",MASTER); 
	

	if(taskid==MASTER)	
		printf("the no. of workers=%ld\n",numworkers.tot);		
	//==========================================================================================//


	//==========================================================================================//
	//building the derived datatypes for discretization information 
	//see "build_datatypes1.c" 
	//this function is called by all processes at this point	
	build_derived_type1(&dom_decomp_info, &offset,&neighbour,&worker,&MPI_DOMAIN_INFO, &MPI_NEIGHBOUR_INFO, &MPI_OFFSET_INFO,&MPI_WORKER_INFO);
	//===========================================================================================//
	//===========================================================================================//
	

	//===================================================================================//
	//computing some derived qunatities unrelated to MPI
	//===================================================================================//
	//initialising the dx and dy (I am going to stick to a regular grid)
	dx=d_sp;
	dy=d_sp;
	//===================================================================================//

	//===================================================================================//
	//computing the value of epsilon
	epsilon=ratio*d_sp;
	if(taskid==MASTER)	
		printf("The value of epsilon for this current simulation=%lf\n",epsilon);
	//===================================================================================//

	//===================================================================================//
	//echoing back the property matrices 

	//==============================================================================//
	//I am going to print out the diffusivity matrices here (see "matrix_printer.c") 
		
	if(taskid==MASTER)	
	{
		printf("The D_S matrix is\n");
		mat_print(D_S);

		printf("The D_L matrix is\n");
		mat_print(D_L);
	}
	//=============================================================================//	 

	//=============================================================================//
	//writing out the equilibrium property values (write_equilibrium.c")
	if(taskid==MASTER)	
		write_eq();
	//=============================================================================//
		
		
	//===================================================================================//
	//calculating the no. of time-steps the simulation is going to run in the
	
	total_tsteps=(long)(total_time/dt);
	if(taskid==MASTER)	
		printf("The total no. of time steps the simulation is going to run is =%ld\n",total_tsteps); 
	//===================================================================================//
	
	//===================================================================================//
	//computing the no. of nodes along x
	nodes_x=(long)(ceil(sine_wave_length/(2*dx)));
	
	printf("The no. of nodes in the x direction=%ld\n",nodes_x);

	//===================================================================================//

	//==================================================================================//
	//computing the total length of the array (including the buffer width)
	//as all computations will be done with the workers (so I will create the buffer layers in the workers) 
	array_length_x=nodes_x;
	array_length_y=nodes_y;

	tot_array_length=array_length_x*array_length_y;
	//==================================================================================//	

	//===================================================================================//
	//computing tau

	//computing "tau" under anti-trapping conditions (see "compute_tau.c")
	comp_tau(taskid);	
	
	if(taskid==MASTER)	
		printf("The value of tau for this current simulation=%lf\n",tau);	

	//====================================================================================//

	//==================================================================================//
	//setting the solid phase length 
	//solid_phase_length=nodes_y/solid_phase_factor;  

	//if(taskid==MASTER)	
		//printf("solid_phase_length=%ld\n",solid_phase_length);	

	//setting the shift length
	//shift_length=nodes_y/shift_factor;
	
	//if(taskid==MASTER)	
		//printf("shift_length=%ld\n",shift_length);
		
 	//==================================================================================//		
		
	//========================================================================================//

	//==============================================================================================//
	//=======COMMON CODE FOR ALL PROCESSES ENDS=====================================================//

	//==============================================================================================//
	//=====coding for the master node===============================================================//
	/* the MASTER node subdivides the problem by grid decomposition,
in this case by rows. The MASTER also initializes the starting values, sets the
boundary conditions on the problem , and receives results from the
nodes after TIME_STEPS */

	

	if(taskid==MASTER)
	{
		
	 
		//===================================================================================//
		//allocating memory to the array that would be relevant to the master 

		if((dep_var=(struct var *)malloc(tot_array_length*sizeof(struct var)))==NULL)
		{
			printf("Space creation failed for the dep_var array\n");
			exit(1);
		}


		//======================================================================================//


		//=======================================================================================//
		//Building new datatypes (see "build_datatypes2.c")====================================================================================//
		//this is being called by the master alone here
		build_derived_type2(dep_var,&MPI_VAR); 
		//=======================================================================================//
		
		#ifdef MULTI
		//=============================================================================//
		//reading the positions (see "input_multi.c")
		mult_input();
		//=============================================================================//
		#endif
		
		//===================================================================================//
		//creating the initial profile (see "compini_multi.c")
		//I am going to set the inital 'mu' and 'fi'profile and other functions
		//if(start_time==0)
        	//{
        	        printf("The code starts from zero\n");
	

        	        //creating the initial composition profile (see "compini_mod4.c")

			compini(dep_var);


	                tsteps=0;

			//=============================================================================//
			//setting the shift count to zero
			shift_count=0;
			//=============================================================================//

	        //}
	
		//else
	        //{
	          //      printf("The code starts from %ld \n",start_time);
	               
			//reading the latest profile 
		//	cont(dep_var,start_time);
	
	        //        tsteps=(long)(start_time/dt);
	
	        //        printf("tsteps=%ld\n",tsteps);

			//=============================================================================//
			//setting the shift count to the supplied value
		//	shift_count=start_shift_count;
			//=============================================================================//
			
	        //}
		//======================================================================================//
		//======================================================================================//	    
		
	
		//=====================================================================================//
		
		//==================================================================================//
		//==================================================================================//
		/* Distribute work to workers. Must first figure out how many rows and columns to
         send and what to do with extra rows and columns. */

		//tackling the rows
      		min_number.rows = nodes_y/numworkers.y; //the no. of rows every worker must handle
    	 	extra.rows = nodes_y%numworkers.y; //the extra rows which can't uniformly distributed to all the worker nodes
		
		//tackling the columns
      		min_number.cols = nodes_x/numworkers.x; //the no. of columns every worker must handle
    	 	extra.cols = nodes_x%numworkers.x; //the extra columns which can't uniformly distributed to all the worker nodes

		printf("min_number.rows=%ld\textra.rows=%ld\n",min_number.rows,extra.rows);
		printf("min_number.cols=%ld\textra.cols=%ld\n",min_number.cols,extra.cols);

		//initialising the offset
		//this carries the information about the starting memory location of the array 
		//important for putting back the array portions correctly  

		offset.x = 0;
		offset.y = 0;

		
		//looping over all the workers 
		//contrary to what I did previously taskid=0 corresponds to a worker instaed of a master
			
		for(j=0;j<numworkers.y;j++)
		{
	     		for (i=0; i<numworkers.x; i++)
	        	{
          
				//===============================================================//
				//===============================================================//

          			/* The following is a particularly compact and efficient  way of distributing the 
        	     		grid. It assures that the number of rows received by each node is only up to one more
        	     		than any other node. It can be read as:
        	     		if (i<=extra_rows) number_rows = min_number_rows + 1
        	    		else number_rows = min_number_rows                   
        	  		*/

				//doing the rows first
        	  		dom_decomp_info.rows = ((j+1) <= extra.rows) ? min_number.rows+1 : min_number.rows;

				//doing the columns
				dom_decomp_info.cols = ((i+1) <= extra.cols) ? min_number.cols+1 : min_number.cols;
				
				//===============================================================//
				//===============================================================//


				//==============================================================//
				//==============================================================//
				//figuring out the rank/worker_number for the process corresponding to the combination of 'i' and 'j'
				worker.num_x=i;
		
				worker.num_y=j;

				worker.number = i+numworkers.x*j;	

				//printf("worker_number=%d\n",worker.number);

				//===============================================================//
				//===============================================================//
		
				//===============================================================//
				/* Tell each worker who its neighbors are, since they must exchange
            			data with each other. */
				
				//for points not on the boundary
				neighbour.left = (i-1)+numworkers.x*j;
				neighbour.right = (i+1)+numworkers.x*j;
				neighbour.up = i+numworkers.x*(j+1);
				neighbour.below = i+numworkers.x*(j-1);
					

				//employing pbc
				//along x
				if(i==0) //blocks which are at the left-most boundary of the system 
					neighbour.left = (numworkers.x-1)+numworkers.x*j;

				if(i==numworkers.x-1)//blocks which are at the right most boundary of the system
					neighbour.right = numworkers.x*j;

				//along y
				if(j==0) //blocks at the bottom most boundary of the system
					neighbour.below = i+numworkers.x*(numworkers.y-1);
				
				if(j==numworkers.y-1) //blocks at the top most boundary of the system
					neighbour.up = i;	

				//printf("after neighbour calculation in the master\n");
				//===============================================================//
				//===============================================================//
		
				//===============================================================//
				//===============================================================//
				//sending over the array to the workers
				//===============================================================//
				//allocating memory to the temporary dependent variable array
				if((temp_dep_var=(struct var *)malloc(dom_decomp_info.rows*dom_decomp_info.cols*sizeof(struct var)))==NULL)
				{
					printf("Space creation failed for the temp_dep_var array for the master\n");
					exit(1);
				}		

				//=============================================================//
				//=============================================================//

				//=============================================================//
				//creating an array which can be directly sent over (see "master_for_sending_over.c")
				cr_arr_for_send_over(temp_dep_var,dom_decomp_info,offset,dep_var,worker.number);
				//=============================================================//

				//============================================================//
				/* Now send startup information to each worker Note that this 
          			information is "tagged" as the BEGIN message*/
         			destination = worker.number;

				//printf("destination=%d\n",destination);
          					
				/* Send the required information to each node */
				MPI_Send(&tsteps, 1, MPI_LONG, destination, BEGIN_TSTEPS, MPI_COMM_WORLD);
				
          			MPI_Send(&worker, 1, MPI_WORKER_INFO, destination, BEGIN_WORKER_INFO, MPI_COMM_WORLD);

	          		MPI_Send(&offset, 1, MPI_OFFSET_INFO, destination, BEGIN_OFFSET_INFO, MPI_COMM_WORLD);

				MPI_Send(&dom_decomp_info, 1, MPI_DOMAIN_INFO, destination, BEGIN_DOMAIN_INFO, MPI_COMM_WORLD);
	
				MPI_Send(&neighbour, 1, MPI_NEIGHBOUR_INFO, destination, BEGIN_NEIGHBOUR_INFO, MPI_COMM_WORLD);

				//printf("destination=%d\n",destination);
					
          			MPI_Send(temp_dep_var,dom_decomp_info.rows*dom_decomp_info.cols, MPI_VAR, destination,BEGIN_DEP_VAR,MPI_COMM_WORLD);

				//printing out the sent information
				printf("Sent to= %d offset_x= %ld offset_y= %ld number_rows= %ld number_columns=%ld neighbour_left= %ld neighbour_right= %ld neighbour_up =%ld neighbour_below =%ld\n",        destination,offset.x,offset.y,dom_decomp_info.rows,dom_decomp_info.cols,	neighbour.left,neighbour.right,neighbour.up,neighbour.below);

				/* increment the offset by the number_rows so the next node will
            			know where its grid begins */
         			offset.x += dom_decomp_info.cols; //updating the x-offset within the i loop
				
				//freeing up the temporary array
				free(temp_dep_var);
					

			} //closing the loop for i

			offset.x =0; //re setting the offset.x
			offset.y += dom_decomp_info.rows;

		}//closing the loop for j
			
		//=======================================================================================//
		//=======================================================================================//

		//===============================================================================//
		//from now on the master coordinates but the workers work 
		for(;tsteps<=total_tsteps;tsteps++)
		{

			//==================================================================//		
			time=tsteps*dt;
	
			//printf("time=%lf\n",time);

			//getchar();
			//==================================================================//

			//======================================================================//
			// checking whether there is a need for shifting 

			tot_shift_flag=0;
	
			//getting the shift information from the workers
			for(j=0;j<numworkers.y;j++)
			{
	     			for(i=0;i<numworkers.x;i++)
	        		{

					source = i+(numworkers.x)*j;

					MPI_Recv(&shift_flag, 1, MPI_INT, source, SHIFT_INFO, MPI_COMM_WORLD,&status);
				
					tot_shift_flag+=shift_flag;
			
				}
			}

			//sending the shift command back to the workers 
			for(j=0;j<numworkers.y;j++)
			{
	     			for(i=0;i<numworkers.x;i++)
	        		{

					destination = i+numworkers.x*j;
				
					MPI_Send(&tot_shift_flag, 1, MPI_INT, destination, TOT_SHIFT_INFO, MPI_COMM_WORLD);
				}
			}

			if(tot_shift_flag>0)
			  shift_count++;	

			

			//=============================================================================//

			//opening the files for writing the composition and eta fields in one and the stress fields in  another
			if(tsteps%output_tsteps_interval==0)
			{
				
				//receiving the computed values from the workers
				
				for(j=0;j<numworkers.y;j++)
				{
	     				for(i=0;i<numworkers.x;i++)
	        			{
       						 source = i+numworkers.x*j;

						
						//============================================//
						//some basic information
         	
       						MPI_Recv(&offset, 1, MPI_OFFSET_INFO, source, END_OFFSET_INFO, MPI_COMM_WORLD,&status);
			
						MPI_Recv(&dom_decomp_info, 1, MPI_DOMAIN_INFO, source, END_DOMAIN_INFO, MPI_COMM_WORLD,&status);
						//============================================//

						//============================================//
						//allocating memory to the temporary arrays
						//the array for the dependent variables 
						if((temp_dep_var=(struct var *)malloc(dom_decomp_info.rows*dom_decomp_info.cols*sizeof(struct var)))==NULL)
						{
							printf("Space creation failed for the temp_dep_var array for the master\n");
							exit(1);
						}			

						//===========================================//	

			
						//===========================================//	
						//receiving the phase field 
						MPI_Recv(temp_dep_var,dom_decomp_info.rows*dom_decomp_info.cols,MPI_VAR,source,END_DEP_VAR,MPI_COMM_WORLD,&status);			
						//===========================================//
						//===========================================//


						//=============================================//
						//reconstructing the domain after copying from the master (see "master_before_writing.c")
						cr_arr_after_recv_from_worker(dep_var,temp_dep_var, offset,dom_decomp_info,source);
						//=============================================//
			
						//=============================================//
						//freeing up the temporary arrays
						free(temp_dep_var);
						//=============================================//
	
					}//ending the neghbour loop along x

				}//ending the neighbour loop along y

				
				//=============================================================//
				//computing the growth rate of the amplitude (see "d_A.c")
				//diff_amp(dep_var,boundary,&peak_x,&peak_y,&trough_x,&trough_y);
	
				//growth_amplitude=(peak_y-trough_y)*d_sp;	

				//=============================================================//

				//writing out the files (see "write_output.c")
				write_out(dep_var,time); 
			}//the 'if' ends here

			

		}//the time loop in the master ends here	

		//============================================================//
		//freeing up the arrays which have been declared in the master (the other temporary arrays are 

		free(dep_var);
	
		
		//============================================================//

	}//The MASTER code ends
	//===================================================================================//
	
	
	//==========================================================================================//
	//===starting worker code===================================================================// 	

	if (taskid < MASTER)
    	{
			
      		/************************* worker code**********************************/

		//=========================================================================//
		//=========================================================================//
		//initial receive from the master==========================================//
		
		//we are going to get the worker no., and the no. of rows information
		source = MASTER;
	
		MPI_Recv(&tsteps, 1, MPI_LONG, source, BEGIN_TSTEPS, MPI_COMM_WORLD,
               &status);
	
	     	MPI_Recv(&worker, 1, MPI_WORKER_INFO, source, BEGIN_WORKER_INFO, MPI_COMM_WORLD,
               &status);

		MPI_Recv(&offset, 1, MPI_OFFSET_INFO, source, BEGIN_OFFSET_INFO, MPI_COMM_WORLD,
               &status);

		MPI_Recv(&dom_decomp_info, 1, MPI_DOMAIN_INFO, source, BEGIN_DOMAIN_INFO, MPI_COMM_WORLD,
               &status);
      
		MPI_Recv(&neighbour, 1, MPI_NEIGHBOUR_INFO, source, BEGIN_NEIGHBOUR_INFO, MPI_COMM_WORLD,
               &status);

		
		
		//===================================================================================//

		//===================================================================================//
		//allocating the memory to the temporary dep_var array
		if((temp_dep_var=(struct var *)malloc(dom_decomp_info.rows*dom_decomp_info.cols*sizeof(struct var)))==NULL)
		{
			printf("Space creation failed for the temp_dep_var array for the worker\n");
			exit(1);
		}			
		//======================================================================================//
		//======================================================================================//
	
		//===================================================================================//
		//===================================================================================//
		//creating domains along with ghost nodes
	
		//the ghost nodes dimensions
		worker_dom_decomp_info.rows = dom_decomp_info.rows+2*BUFFER_WIDTH;
		worker_dom_decomp_info.cols = dom_decomp_info.cols+2*BUFFER_WIDTH;
		
			
		//allocating memory in the workers
		//the dependent variable array
		if((dep_var=(struct var *)malloc(worker_dom_decomp_info.rows*worker_dom_decomp_info.cols*sizeof(struct var)))==NULL)
		{
			printf("Space creation failed for the ph_inter array for the worker\n");
			exit(1);
		}

		if((info=(struct grad_info *)malloc(worker_dom_decomp_info.rows*worker_dom_decomp_info.cols*sizeof(struct grad_info)))==NULL)
		{
			printf("Space creation failed for the info array\n");
			exit(1);
		}	

		
		//====================================================================================//
		//====================================================================================//
	
		//=======================================================================================//
		//Building new datatypes (see "build_datatypes2.c")						
		//====================================================================//
		//this is being called by the master alone here
		build_derived_type2(dep_var,&MPI_VAR); 
		//=======================================================================================//


		//======================================================================================//
		//======================================================================================//
		//receiving the different arrays

		//the dependent variable
		MPI_Recv(temp_dep_var, dom_decomp_info.rows*dom_decomp_info.cols, MPI_VAR, source, BEGIN_DEP_VAR,MPI_COMM_WORLD, &status);

		printf("at=%d\tsource=%d\n",taskid,source);

		//=======================================================================================//
		//=======================================================================================//

		//====================================================================================//
		//====================================================================================//
		//transferring the contents from the temporary array to the the ones with buffers layers (see "worker_transfer_from_temp.c")
		cr_worker_transfer_from_temp(temp_dep_var,dom_decomp_info,worker_dom_decomp_info,dep_var,taskid);
		//====================================================================================//
		//====================================================================================//

		//====================================================================================//
		//====================================================================================//
		//creating the temporary arrays for doing the sends to the neighbouring blocks for phase field and interpolants
		struct var *dep_var2left,*dep_var2right,*dep_var2up,*dep_var2below;

		//allocating memory to the temporary arrays for sends to the neighbouring blocks
		if((dep_var2left=(struct var *)malloc(worker_dom_decomp_info.rows*BUFFER_WIDTH*sizeof(struct var)))==NULL)
		{
			printf("Space creation failed for the dep_var2left array for the worker\n");
			exit(1);
		}

		if((dep_var2right=(struct var *)malloc(worker_dom_decomp_info.rows*BUFFER_WIDTH*sizeof(struct var)))==NULL)
		{
			printf("Space creation failed for the dep_var2right array for the worker\n");
			exit(1);
		}
		if((dep_var2up=(struct var *)malloc(worker_dom_decomp_info.cols*BUFFER_WIDTH*sizeof(struct var)))==NULL)
		{
			printf("Space creation failed for the dep_var2up array for the worker\n");
			exit(1);
		}
		if((dep_var2below=(struct var *)malloc(worker_dom_decomp_info.cols*BUFFER_WIDTH*sizeof(struct var)))==NULL)
		{
			printf("Space creation failed for the dep_var2below array for the worker\n");
			exit(1);
		}

		//====================================================================================//
		//====================================================================================//
		//creating the temporary arrays for doing the receives from the neighbouring blocks
		struct var *dep_var_fr_left,*dep_var_fr_right,*dep_var_fr_up,*dep_var_fr_below;
		
		//allocating memory to the temporary arrays for receives from the neighbouring blocks
		if((dep_var_fr_left=(struct var *)malloc(worker_dom_decomp_info.rows*BUFFER_WIDTH*sizeof(struct var)))==NULL)
		{
			printf("Space creation failed for the dep_var_fr_left array for the worker\n");
			exit(1);
		}

		if((dep_var_fr_right=(struct var *)malloc(worker_dom_decomp_info.rows*BUFFER_WIDTH*sizeof(struct var)))==NULL)
		{
			printf("Space creation failed for the dep_var_fr_right array for the worker\n");
			exit(1);
		}
		if((dep_var_fr_up=(struct var *)malloc(worker_dom_decomp_info.cols*BUFFER_WIDTH*sizeof(struct var)))==NULL)
		{
			printf("Space creation failed for the dep_var_fr_up array for the worker\n");
			exit(1);
		}
		if((dep_var_fr_below=(struct var *)malloc(worker_dom_decomp_info.cols*BUFFER_WIDTH*sizeof(struct var)))==NULL)
		{
			printf("Space creation failed for the dep_var_fr_below array for the worker\n");
			exit(1);
		}
	

		//=====================================================================================//
		//=====================================================================================//

		//=====================================================================================//
		//computing the exchange lengths for left-right
		exchange_length_left_right= worker_dom_decomp_info.rows*BUFFER_WIDTH;
		
		//computing the exchange lengths for up-below
		exchange_length_up_below= worker_dom_decomp_info.cols*BUFFER_WIDTH;

		//====================================================================================//


		//starting the time iteration loop in the workers 
		for(;tsteps<=total_tsteps;tsteps++)
		{

			//==================================================================//
			time=tsteps*dt;
	
			//printf("time=%lf\n",time);

			//getchar();
			//=================================================================//

			//=================================================================//		
			//we have to exchange the phase field at the buffer layers
			//see "exchange_dep_var.c"
			ex_dep_var(dep_var,dep_var2left,dep_var2right,dep_var2up,dep_var2below, dep_var_fr_left,dep_var_fr_right,dep_var_fr_up,dep_var_fr_below,dom_decomp_info,worker_dom_decomp_info, worker, neighbour, numworkers);
			//=================================================================//


		
			//==================================================================//
				
			//computing the gradients of 'fi' and 'mu'(see "compute_gradients.c")
			comp_grad(info,dep_var,worker_dom_decomp_info);
		
			//computing the divergence (see "divergence_fi.c")
			comp_diver_fi(info,dep_var,worker_dom_decomp_info,tsteps);

			//computing the anti-trapping fluxes  (see "compute_anti_trapping_currents.c")
			comp_anti_trap(info,dep_var,worker_dom_decomp_info); 

			//updating 'fi' and 'mu' (see "update_fi_mu.c")

			upd_mu_fi(info,dep_var,worker_dom_decomp_info,tsteps);

			//=======================================================================//
			//the shifting check (see "check_for_shifting.c")
			shift_flag=ch_shift(dep_var,offset,dom_decomp_info,worker_dom_decomp_info);	
		

			//sending over the shift flags to the master
			MPI_Send(&shift_flag, 1, MPI_INT, MASTER, SHIFT_INFO, MPI_COMM_WORLD);

			//getting the total shift flags from the master
 			MPI_Recv(&tot_shift_flag, 1, MPI_INT, MASTER, TOT_SHIFT_INFO, MPI_COMM_WORLD,&status);

			//doing the shift
			if(tot_shift_flag>0)
			{
				//doing neighbour exchange in the vertical direction
				//see "vert_exchange_dep_var.c"
				vert_ex_dep_var(dep_var,dep_var2up,dep_var2below,dep_var_fr_up,dep_var_fr_below,dom_decomp_info,worker_dom_decomp_info, worker, neighbour, numworkers);

				//shifting the profiles (see "comp_shift_profiles.c")
				shift_profile(dep_var,worker_dom_decomp_info,worker);
			}
			//=============================================================================//
			
			//=============================================================================//
	
			//sending over the computed values to the master 
			if(tsteps%output_tsteps_interval==0)
			{
			
					
				//moving the arrays into a no-ghost node setting (see "worker_for_sending_over.c") 
				cr_arr_for_send_over_fr_wrkr(temp_dep_var,dom_decomp_info,worker_dom_decomp_info,dep_var,taskid);	

				//printf("after array rearrange from the worker = %d\n",taskid);		

				MPI_Send(&offset, 1, MPI_OFFSET_INFO, MASTER, END_OFFSET_INFO, MPI_COMM_WORLD);

	      			MPI_Send(&dom_decomp_info, 1, MPI_DOMAIN_INFO, MASTER, END_DOMAIN_INFO, MPI_COMM_WORLD);

	     			MPI_Send(temp_dep_var, dom_decomp_info.rows*dom_decomp_info.cols, MPI_VAR, MASTER, END_DEP_VAR, MPI_COMM_WORLD);

				
				
			}//the 'if' for writing the files	
			//=======================================================================//

		
			
									

	

		}//the 'time' loop in the workers

		//freeing up the arrays in the workers
		free(dep_var);
		free(info);
		
		free(temp_dep_var);

		free(dep_var2left);
		free(dep_var2right);
		free(dep_var2up);
		free(dep_var2below);
		

		free(dep_var_fr_left);
		free(dep_var_fr_right);
		free(dep_var_fr_up);
		free(dep_var_fr_below);
		
		
	}//the worker code ends	
	/*gracefully exit MPI */
	printf("Exiting mpi by taskid=%d\n",taskid);
		
  	MPI_Finalize();	
	
}
