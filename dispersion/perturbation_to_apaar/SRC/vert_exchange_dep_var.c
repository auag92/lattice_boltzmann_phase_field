void vert_ex_dep_var(struct var *dep_var,struct var *dep_var2up,struct var *dep_var2below,struct var *dep_var_fr_up,struct var *dep_var_fr_below,struct domain_info dom_decomp_info,struct domain_info worker_dom_decomp_info,struct worker_info worker,struct neighbour_info neighbour,struct no_of_workers numworkers)
{
	//================================================================================================//
	//=====================UP-BELOW===================================================================//

	//=============================================================================================//
	//setting up the temporary arrays for up below send 	
	//see "cr_temp_dep_var_for_up_below_sends.c"
	cr_temp_dep_var_for_up_below_sends(dep_var,dep_var2up,dep_var2below,dom_decomp_info,worker_dom_decomp_info);	
 	//=============================================================================================//
	

	//=============================================================================================//
	//exchanging the neighbour information for phase_field and interpolants array (see "neighbour_exchange_dep_var_up_below.c") UP-BELOW exchange

	exchange_neigh_dep_var_up_below(worker,neighbour,dep_var2up,dep_var2below,dep_var_fr_up,dep_var_fr_below);	
	//==============================================================================================//
	
	//================================================================================//
	//populating the worker from the temporary arrays for receives (see "cr_temp_dep_var_for_up_below_recvs.c")
	cr_temp_dep_var_for_up_below_recvs(dep_var,dep_var_fr_up,dep_var_fr_below,dom_decomp_info,worker_dom_decomp_info,worker); 
	//================================================================================//

	//==============================================================================================//
	//===========================================================================================//	
}	
