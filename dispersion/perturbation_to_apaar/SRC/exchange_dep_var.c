void ex_dep_var(struct var *dep_var,struct var *dep_var2left,struct var *dep_var2right,struct var *dep_var2up,struct var *dep_var2below, struct var *dep_var_fr_left,struct var *dep_var_fr_right,struct var *dep_var_fr_up,struct var *dep_var_fr_below,struct domain_info dom_decomp_info,struct domain_info worker_dom_decomp_info,struct worker_info worker,struct neighbour_info neighbour,struct no_of_workers numworkers)
{
	
	
	//==============================================================================================//
	//============LEFT_RIGHT========================================================================//

	//=============================================================================================//
	//setting up the temporary arrays for left right send 	
	//see "cr_temp_dep_var_for_left_right_sends.c"
	cr_temp_dep_var_for_left_right_sends(dep_var,dep_var2left,dep_var2right,dom_decomp_info,worker_dom_decomp_info);	
 	//=============================================================================================//

	//=============================================================================================//
	//exchanging the neighbour information for phase_field and interpolants array (see "neighbour_exchange_dep_var_left_right.c") LEFT-RIGHT exchange

	exchange_neigh_dep_var_left_right(worker,neighbour,dep_var2left,dep_var2right,dep_var_fr_left,dep_var_fr_right);	
	//==============================================================================================//

	
	//================================================================================//
	//populating the worker from the temporary arrays for receives (see "cr_temp_dep_var_for_left_right_recvs.c")
	cr_temp_dep_var_for_left_right_recvs(dep_var,dep_var_fr_left,dep_var_fr_right,dom_decomp_info,worker_dom_decomp_info,worker); 
	//================================================================================//

	//==============================================================================================//
	//==============================================================================================//

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
	//==============================================================================================//



}
		
