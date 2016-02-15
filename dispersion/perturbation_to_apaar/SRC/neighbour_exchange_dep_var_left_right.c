void exchange_neigh_dep_var_left_right(struct worker_info worker,struct neighbour_info neighbour,struct var *dep_var2left,struct var *dep_var2right,struct var *dep_var_fr_left,struct var *dep_var_fr_right)
{

		if (worker.num_x %2) //a process with an odd no. as a taskid along x 
		{
			if (worker.num_x != 0) //not at the left_most boundary of the system  and exchanging information with the lower ranked process along the x-direction which has an even taskid
			{

				//the send
				
			
				MPI_Send(dep_var2left, exchange_length_left_right, MPI_VAR, neighbour.left, NGHBOR_LEFT, MPI_COMM_WORLD);

				
				
				//the receive 

				

				MPI_Recv(dep_var_fr_left, exchange_length_left_right, MPI_VAR, neighbour.left, NGHBOR_RIGHT, MPI_COMM_WORLD, &status);
				

				
						 
			}
	
			if (worker.num_x !=numworkers.x-1) //not at the right most boundary and exchanging information with the higher ranked process along the x-direction which has an even rank   
			{

				//the send
				
			
				MPI_Send(dep_var2right, exchange_length_left_right, MPI_VAR, neighbour.right, NGHBOR_RIGHT, MPI_COMM_WORLD);

				//the receive
				
			
				MPI_Recv(dep_var_fr_right, exchange_length_left_right, MPI_VAR, neighbour.right, NGHBOR_LEFT, MPI_COMM_WORLD, &status);

				
			  	
			}
		} 

		else //when the taskid is divisible by 2 
		{
			if (worker.num_x !=numworkers.x-1) //not at the right most boundary  and so exchanging information with the higher rank process along the x which has an odd rank  
			{
				//the order of receieves and sends here is dictated by the odd ranked processes
				
				//the receive
				
			
				MPI_Recv(dep_var_fr_right, exchange_length_left_right, MPI_VAR, neighbour.right, NGHBOR_LEFT, MPI_COMM_WORLD, &status);

				

				//the send
				
			
				MPI_Send(dep_var2right, exchange_length_left_right, MPI_VAR, neighbour.right, NGHBOR_RIGHT, MPI_COMM_WORLD);

			
			  
			}

			if(worker.num_x != 0)   //not at the left most boundary and so exchanging information with the lower rank process along the x which has an odd rank  
			{
				//the order of receieves and sends here is dictated by the odd ranked processes
				//the receive
				
				
		
				MPI_Recv(dep_var_fr_left, exchange_length_left_right, MPI_VAR, neighbour.left, NGHBOR_RIGHT, MPI_COMM_WORLD, &status);
			
						 


			 	//the send
				
			
				MPI_Send(dep_var2left,exchange_length_left_right, MPI_VAR, neighbour.left, NGHBOR_LEFT, MPI_COMM_WORLD);

				

			}
		}
		//===================================================================================//
		//===================================================================================//	

		//===================================================================================//
		//===================================================================================//	
		//exchanging information for the boundary elements 
		//this is done separately because we are not sure whether 'worker.num_x' is an even or odd no. 						
		//suppose "worker.num_x" is an odd no. then the above coordination breaks down
			
		if(worker.num_x==0) //at the left-most boundary 
		{
				
				//the first send
				
			
				MPI_Send(dep_var2left,exchange_length_left_right , MPI_VAR, neighbour.left, NGHBOR_LEFT, MPI_COMM_WORLD);

			
				//the receive
				
				
		
				MPI_Recv(dep_var_fr_left, exchange_length_left_right, MPI_VAR, neighbour.left, NGHBOR_RIGHT, MPI_COMM_WORLD, &status);
			
				
		}

		if(worker.num_x ==numworkers.x-1) //at the right most boundary
		{
				//the receive
				
			
				MPI_Recv(dep_var_fr_right,exchange_length_left_right, MPI_VAR, neighbour.right, NGHBOR_LEFT, MPI_COMM_WORLD, &status);

			

				//the send
				
			
				MPI_Send(dep_var2right, exchange_length_left_right, MPI_VAR, neighbour.right, NGHBOR_RIGHT, MPI_COMM_WORLD);

				
		}


}
