//this routine will contain the different matrix operations

//printing out a square matrix
void mat_print(double (*mat)[NUMCOMPONENTS-1])
{
	long m,n;
	
	for(m=0;m<NUMCOMPONENTS-1;m++) 	
	{
		for(n=0;n<NUMCOMPONENTS-1;n++) 	
		{ 	
			printf("%lf ",mat[m][n]);
		}
		printf("\n");
	}
}

//the product of 2 square matrices
void mat_mat_mult(double (*prod)[NUMCOMPONENTS-1],double (*mat_1)[NUMCOMPONENTS-1],double (*mat_2)[NUMCOMPONENTS-1])
{
	long m,n,p;	

	for(m=0;m<NUMCOMPONENTS-1;m++) 	
	{
		for(n=0;n<NUMCOMPONENTS-1;n++) 	
		{ 	
			prod[m][n]=0.0;

			for(p=0;p<NUMCOMPONENTS-1;p++)
			{
				prod[m][n]+=mat_1[m][p]*mat_2[p][n];	
			}
		}
	}
}

void mat_inv(double (*inv_mat)[NUMCOMPONENTS-1],double (*mat)[NUMCOMPONENTS-1])//as the formulae are hard coded we don't need to pass the dimensional equations 
{
	//I am going to code the explicit formulae
	
	#ifdef BINARY
	if(fabs(mat[0][0])<1e-16)
	{
		printf("The matrix can't be inverted\n");
		exit(1);
	}
	
	inv_mat[0][0]=1.0/mat[0][0];
	#endif

	#ifdef TERNARY
	double det=mat[0][0]*mat[1][1]-mat[1][0]*mat[0][1];

	if(fabs(det)<1e-16)
	{
		printf("The matrix can't be inverted\n");
		exit(1);
	}

	inv_mat[0][0]=mat[1][1]/det;

	inv_mat[1][1]=mat[0][0]/det;
		
	inv_mat[0][1]=-mat[0][1]/det;

	inv_mat[1][0]=-mat[1][0]/det;	
	#endif

	#ifdef QUARTERNARY
	double det=mat[0][0]*mat[1][1]*mat[2][2]+mat[1][0]*mat[2][1]*mat[0][2]+mat[2][0]*mat[0][1]*mat[1][2]-mat[0][0]*mat[2][1]*mat[1][2]-mat[2][0]*mat[1][1]*mat[0][2]-mat[1][0]*mat[0][1]*mat[2][2];
	
	if(fabs(det)<1e-16)
	{
		printf("The matrix can't be inverted\n");
		exit(1);
	}

	inv_mat[0][0]= (mat[2][2]*mat[1][1]-mat[2][1]*mat[1][2])/det;
	
	inv_mat[0][1]=-(mat[2][2]*mat[0][1]-mat[2][1]*mat[0][2])/det;

	inv_mat[0][2]= (mat[1][2]*mat[0][1]-mat[1][1]*mat[0][2])/det;

	inv_mat[1][0]=-(mat[2][2]*mat[1][0]-mat[2][0]*mat[1][2])/det;	

	inv_mat[1][1]= (mat[2][2]*mat[0][0]-mat[2][0]*mat[0][2])/det;
	
	inv_mat[1][2]=-(mat[1][2]*mat[0][0]-mat[1][0]*mat[0][2])/det;

	inv_mat[2][0]= (mat[2][1]*mat[1][0]-mat[2][0]*mat[1][1])/det;	

	inv_mat[2][1]=-(mat[2][1]*mat[0][0]-mat[2][0]*mat[0][1])/det;

	inv_mat[2][2]= (mat[1][1]*mat[0][0]-mat[1][0]*mat[0][1])/det;
		
	#endif
}
//a matrix-vector product
void mat_vec_mult(double prod[NUMCOMPONENTS-1],double (*mat)[NUMCOMPONENTS-1],double vec[NUMCOMPONENTS-1])
{
		
	long n,p;	
	
	for(n=0;n<NUMCOMPONENTS-1;n++) 	
	{ 	
		prod[n]=0.0;

		for(p=0;p<NUMCOMPONENTS-1;p++)
		{
			prod[n]+=mat[n][p]*vec[p];	
		}
	}
}

//the vector inner product	
double vec_inn_prod(double vec_1[NUMCOMPONENTS-1],double vec_2[NUMCOMPONENTS-1])
{
	double val=0.0;

	long p;
	
	for(p=0;p<NUMCOMPONENTS-1;p++)
	{
		val+=vec_1[p]*vec_2[p];
	}
	
	return(val);
}	
