//The constants related to the mpi are mentioned here 
//the parameters related to initial transfer from master to the workers
#define BEGIN_WORKER_INFO 1
#define BEGIN_OFFSET_INFO 2
#define BEGIN_DOMAIN_INFO 3 
#define BEGIN_NEIGHBOUR_INFO 4
#define BEGIN_DEP_VAR 5


//the parameters relevant to buffer layer exchange
#define NGHBOR_LEFT 6
#define NGHBOR_RIGHT 7
#define NGHBOR_UP 8
#define NGHBOR_BELOW 9

//the parameters relevant to transfer from workers to master 
#define END_OFFSET_INFO 10
#define END_DOMAIN_INFO 11
#define END_DEP_VAR 12
#define SHIFT_INFO 13
#define TOT_SHIFT_INFO 14

#define BEGIN_TSTEPS 15

//=======================================================================================================//
//creating some structure types
struct no_of_workers {
	long x;
	long y;	
	long tot;
	};

struct no_of_workers numworkers; //the variable containing the no. of workers 

struct offset_info {
	long x;
	long y;
	};

struct domain_info {
	long cols;
	long rows;
	};

struct neighbour_info {
	long left;
	long right;
	long up;	
	long below;
	};

struct worker_info {
	long num_x;
	long num_y; 
	long number;	
	};
//========================================================================================================//
//some datatypes for exchanging with MPI
MPI_Datatype MPI_DOMAIN_INFO, MPI_NEIGHBOUR_INFO, MPI_OFFSET_INFO,MPI_WORKER_INFO;

MPI_Datatype MPI_VAR; //the datatypes which will be sent
//========================================================================================================//

long exchange_length_left_right,exchange_length_up_below;

MPI_Status status; //status variable

int MASTER;//the taskid corresponding to the master
