void build_derived_type2(struct var *dep_var, MPI_Datatype *MPI_VAR)
{
	//===========================================================================================//
	//===========================================================================================//
	//creating the datatype MPI_VAR

	int block_lengths_var[2];
	MPI_Aint displacements_var[2];
	MPI_Datatype typelists_var[2];
	MPI_Aint start_address_var;
	MPI_Aint fin_address_var;

	//fill block lengths
	block_lengths_var[0] = 1; //we are looking at a single double
	block_lengths_var[1] = NUMCOMPONENTS-1; //we are looking at an array of the type "double"
	
	//Fill Typelists
	typelists_var[0] = MPI_DOUBLE;
	typelists_var[1] = MPI_DOUBLE;
	
	//Calculate Displacements
	displacements_var[0] = 0;

	MPI_Get_address(&((*dep_var).fi),&start_address_var);
	MPI_Get_address(&((*dep_var).mu[0]),&fin_address_var);

	displacements_var[1] = fin_address_var - start_address_var;

	//Create and Commit new mpi type
	MPI_Type_create_struct(2,block_lengths_var,displacements_var,typelists_var,MPI_VAR);
	MPI_Type_commit(MPI_VAR);
	//==============================================================================================//
	//==============================================================================================//


}
