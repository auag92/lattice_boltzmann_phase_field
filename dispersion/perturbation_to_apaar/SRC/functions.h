//computing g
double g(double fi)
{
	return(fi);

	//return(fi*fi*(3.0-2.0*fi));
}

//computing h
double h(double fi)
{
	return(fi*fi*(3.0-2.0*fi));
}

//computing f
double f(double fi)
{
	return(fi*fi*(1.0-fi)*(1.0-fi));
}

//computing d_h_d_fi
double d_h_d_fi(double fi)
{
	return(6.0*fi*(1.0-fi));
}

//computing d_f_d_fi
double d_f_d_fi(double fi)
{
	return(2.0*fi*(1.0-fi)*(1.0-2.0*fi));
}

//computing fi_times_one_min_fi
double fi_times_one_min_fi(double fi)
{
	return(fi*(1.0-fi));
}



