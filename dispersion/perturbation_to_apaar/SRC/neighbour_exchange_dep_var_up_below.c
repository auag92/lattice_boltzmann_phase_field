void exchange_neigh_dep_var_up_below(struct worker_info worker,struct neighbour_info neighbour,struct var *dep_var2up,struct var *dep_var2below,struct var *dep_var_fr_up,struct var *dep_var_fr_below)
{
		
		
		if (worker.num_y %2) //a process with an odd no. as a taskid along y 
		{
			if (worker.num_y != 0) //not at the lower most boundary of the system  and exchanging information with the lower ranked process along the y-direction which has an even taskid
			{

				//the send
				
			
				MPI_Send(dep_var2below, exchange_length_up_below, MPI_VAR, neighbour.below, NGHBOR_BELOW, MPI_COMM_WORLD);

				
				
				//the receive 

				

				MPI_Recv(dep_var_fr_below, exchange_length_up_below, MPI_VAR, neighbour.below, NGHBOR_UP, MPI_COMM_WORLD, &status);
				

				
						 
			}
	
			if (worker.num_y !=numworkers.y-1) //not at the upper most boundary and exchanging information with the higher ranked process along the y-direction which has an even rank   
			{

				//the send
				
			
				MPI_Send(dep_var2up, exchange_length_up_below, MPI_VAR, neighbour.up, NGHBOR_UP, MPI_COMM_WORLD);

				//the receive
				
			
				MPI_Recv(dep_var_fr_up, exchange_length_up_below, MPI_VAR, neighbour.up, NGHBOR_BELOW, MPI_COMM_WORLD, &status);

				
			  	
			}
		} 

		else //when the taskid is divisible by 2 
		{
			if (worker.num_y !=numworkers.y-1) //not at the uppermost boundary  and so exchanging information with the higher rank process along the y which has an odd rank  
			{
				//the order of receieves and sends here is dictated by the odd ranked processes
				
				//the receive
				
			
				MPI_Recv(dep_var_fr_up, exchange_length_up_below, MPI_VAR, neighbour.up, NGHBOR_BELOW, MPI_COMM_WORLD, &status);

				

				//the send
				
			
				MPI_Send(dep_var2up, exchange_length_up_below, MPI_VAR, neighbour.up, NGHBOR_UP, MPI_COMM_WORLD);

			
			  
			}

			if(worker.num_y != 0)   //not at the lowermost boundary and so exchanging information with the lower rank process along the y which has an odd rank  
			{
				//the order of receieves and sends here is dictated by the odd ranked processes
				//the receive
				
				
		
				MPI_Recv(dep_var_fr_below, exchange_length_up_below, MPI_VAR, neighbour.below, NGHBOR_UP, MPI_COMM_WORLD, &status);
			
						 


			 	//the send
				
			
				MPI_Send(dep_var2below,exchange_length_up_below, MPI_VAR, neighbour.below, NGHBOR_BELOW, MPI_COMM_WORLD);

				

			}
		}
		//===================================================================================//
		//===================================================================================//	

		//===================================================================================//
		//===================================================================================//	
		//exchanging information for the boundary elements 
		//this is done separately because we are not sure whether 'worker.num_y' is an even or odd no. 						
		//suppose "worker.num_y" is an odd no. then the above coordination breaks down
			
		if(worker.num_y==0) //at the bottom-most boundary 
		{
				
				//the first send
				
			
				MPI_Send(dep_var2below,exchange_length_up_below , MPI_VAR, neighbour.below, NGHBOR_BELOW, MPI_COMM_WORLD);

			
				//the receive
				
				
		
				MPI_Recv(dep_var_fr_below, exchange_length_up_below, MPI_VAR, neighbour.below, NGHBOR_UP, MPI_COMM_WORLD, &status);
			
				
		}

		if(worker.num_y ==numworkers.y-1) //at the upper most boundary
		{
				//the receive
				
			
				MPI_Recv(dep_var_fr_up,exchange_length_up_below, MPI_VAR, neighbour.up, NGHBOR_BELOW, MPI_COMM_WORLD, &status);

			

				//the send
				
			
				MPI_Send(dep_var2up, exchange_length_up_below, MPI_VAR, neighbour.up, NGHBOR_UP, MPI_COMM_WORLD);

				
		}


}
