//building the derived datatype for u_n 
void build_derived_type1(struct domain_info *dom_decomp_info, struct offset_info *offset,struct neighbour_info *neighbour,struct worker_info *worker,MPI_Datatype *MPI_DOMAIN_INFO, MPI_Datatype *MPI_NEIGHBOUR_INFO, MPI_Datatype *MPI_OFFSET_INFO,MPI_Datatype *MPI_WORKER_INFO)
{

	
	//===================================================================================//
	//===================================================================================//
	//creating the datatype MPI_DOMAIN_INFO

	int block_lengths_DOMAIN_INFO[2];
	MPI_Aint displacements_DOMAIN_INFO[2];
	MPI_Datatype typelists_DOMAIN_INFO[2];
	MPI_Aint start_address_DOMAIN_INFO;
	MPI_Aint fin_address_DOMAIN_INFO;
	//fill block lengths
	block_lengths_DOMAIN_INFO[0] = 1; //we are looking at a single long
	block_lengths_DOMAIN_INFO[1] = 1; //we are looking at a single long
	
	//Fill Typelists
	typelists_DOMAIN_INFO[0] = MPI_LONG;
	typelists_DOMAIN_INFO[1] = MPI_LONG;
	
	//Calculate Displacements
	displacements_DOMAIN_INFO[0] = 0;

	MPI_Get_address(&((*dom_decomp_info).cols),&start_address_DOMAIN_INFO);
	MPI_Get_address(&((*dom_decomp_info).rows),&fin_address_DOMAIN_INFO);
	

	displacements_DOMAIN_INFO[1] = fin_address_DOMAIN_INFO - start_address_DOMAIN_INFO;

		
	//Create and Commit new mpi type
	MPI_Type_create_struct(2,block_lengths_DOMAIN_INFO,displacements_DOMAIN_INFO,typelists_DOMAIN_INFO,MPI_DOMAIN_INFO);
	MPI_Type_commit(MPI_DOMAIN_INFO);

	//===================================================================================//
	//===================================================================================//



	//========================================================================================//
	//========================================================================================//

	//creating the datatype MPI_NEIGHBOUR_INFO

	int block_lengths_NEIGHBOUR_INFO[4];
	MPI_Aint displacements_NEIGHBOUR_INFO[4];
	MPI_Datatype typelists_NEIGHBOUR_INFO[4];
	
	MPI_Aint address1_NEIGHBOUR_INFO;
	MPI_Aint address2_NEIGHBOUR_INFO;
	MPI_Aint address3_NEIGHBOUR_INFO;
	MPI_Aint address4_NEIGHBOUR_INFO;
	
	//fill block lengths
	block_lengths_NEIGHBOUR_INFO[0] = 1; //we are looking at a single long
	block_lengths_NEIGHBOUR_INFO[1] = 1; //we are looking at a single long
	block_lengths_NEIGHBOUR_INFO[2] = 1; //we are looking at a single long
	block_lengths_NEIGHBOUR_INFO[3] = 1; //we are looking at a single long
	
	//Fill Typelists
	typelists_NEIGHBOUR_INFO[0] = MPI_LONG;
	typelists_NEIGHBOUR_INFO[1] = MPI_LONG;
	typelists_NEIGHBOUR_INFO[2] = MPI_LONG;
	typelists_NEIGHBOUR_INFO[3] = MPI_LONG;
	
	//Calculate Displacements
	displacements_NEIGHBOUR_INFO[0] = 0;

	MPI_Get_address(&((*neighbour).left),&address1_NEIGHBOUR_INFO);
	MPI_Get_address(&((*neighbour).right),&address2_NEIGHBOUR_INFO);
	MPI_Get_address(&((*neighbour).up),&address3_NEIGHBOUR_INFO);
	MPI_Get_address(&((*neighbour).below),&address4_NEIGHBOUR_INFO);


	displacements_NEIGHBOUR_INFO[1] = address2_NEIGHBOUR_INFO - address1_NEIGHBOUR_INFO;
	displacements_NEIGHBOUR_INFO[2] = address3_NEIGHBOUR_INFO - address1_NEIGHBOUR_INFO;
	displacements_NEIGHBOUR_INFO[3] = address4_NEIGHBOUR_INFO - address1_NEIGHBOUR_INFO;
	
	//Create and Commit new mpi type
	MPI_Type_create_struct(4,block_lengths_NEIGHBOUR_INFO,displacements_NEIGHBOUR_INFO, typelists_NEIGHBOUR_INFO,MPI_NEIGHBOUR_INFO);
	MPI_Type_commit(MPI_NEIGHBOUR_INFO);

	//===================================================================================//
	//===================================================================================//



	//===================================================================================//
	//===================================================================================//

	//creating the datatype MPI_OFFSET_INFO

	int block_lengths_OFFSET_INFO[2];
	MPI_Aint displacements_OFFSET_INFO[2];
	MPI_Datatype typelists_OFFSET_INFO[2];
	MPI_Aint start_address_OFFSET_INFO;
	MPI_Aint fin_address_OFFSET_INFO;
	//fill block lengths
	block_lengths_OFFSET_INFO[0] = 1; //we are looking at a single long
	block_lengths_OFFSET_INFO[1] = 1; //we are looking at a single long
	
	//Fill Typelists
	typelists_OFFSET_INFO[0] = MPI_LONG;
	typelists_OFFSET_INFO[1] = MPI_LONG;
	
	//Calculate Displacements
	displacements_OFFSET_INFO[0] = 0;

	MPI_Get_address(&((*offset).x),&start_address_OFFSET_INFO);
	MPI_Get_address(&((*offset).y),&fin_address_OFFSET_INFO);

	displacements_OFFSET_INFO[1] = fin_address_OFFSET_INFO - start_address_OFFSET_INFO;

	
	//Create and Commit new mpi type
	MPI_Type_create_struct(2,block_lengths_OFFSET_INFO,displacements_OFFSET_INFO,typelists_OFFSET_INFO,MPI_OFFSET_INFO);
	MPI_Type_commit(MPI_OFFSET_INFO);


	//===================================================================================//
	//===================================================================================//


	//===================================================================================//
	//===================================================================================//

	//creating the datatype MPI_WORKER_INFO

	int block_lengths_WORKER_INFO[3];
	MPI_Aint displacements_WORKER_INFO[3];
	MPI_Datatype typelists_WORKER_INFO[3];
	MPI_Aint address1_WORKER_INFO;
	MPI_Aint address2_WORKER_INFO;
	MPI_Aint address3_WORKER_INFO;
	
	
	//fill block lengths
	block_lengths_WORKER_INFO[0] = 1; //we are looking at a single long
	block_lengths_WORKER_INFO[1] = 1; //we are looking at a single long
	block_lengths_WORKER_INFO[2] = 1; //we are looking at a single long
	
	//Fill Typelists
	typelists_WORKER_INFO[0] = MPI_LONG;
	typelists_WORKER_INFO[1] = MPI_LONG;
	typelists_WORKER_INFO[2] = MPI_LONG;
	
	//Calculate Displacements
	displacements_WORKER_INFO[0] = 0;

	MPI_Get_address(&((*worker).num_x),&address1_WORKER_INFO);
	MPI_Get_address(&((*worker).num_y),&address2_WORKER_INFO);
	MPI_Get_address(&((*worker).number),&address3_WORKER_INFO);

	displacements_WORKER_INFO[1] = address2_WORKER_INFO - address1_WORKER_INFO;
	displacements_WORKER_INFO[2] = address3_WORKER_INFO - address1_WORKER_INFO;
	
	//Create and Commit new mpi type
	MPI_Type_create_struct(3,block_lengths_WORKER_INFO,displacements_WORKER_INFO,typelists_WORKER_INFO,MPI_WORKER_INFO);
	MPI_Type_commit(MPI_WORKER_INFO);

	//===================================================================================//
	//===================================================================================//

	
}
