#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "phi_definitions.c"
#include "lbm_definitions.c"
#include "file_io.c"
#include "mpi.h"
#include "time.h"
//****************************************************************************//
#define MASTER  0
#define NONE    0
#define BEGIN   111
#define LTAG    222
#define RTAG    333
#define WRITE   444
#define ERROR   555
#define BREAK   666
//****************************************************************************//
void allocate_rows();
void allocate_memory();
void free_mem();
void start_clock();
void end_clock();
//*****************************MPI variables**********************************//
int numtasks, numworkers, taskid, rank, dest;
int averow, extra, offset;
int left_node, right_node;
int offset_ax, offset_ay;
int source, msgtype;
int rows;
MPI_Status status;
//*****************************************************************************//
clock_t start_t, end_t;
double total_t;
//****************************************************************************//
int t;
//****************************************************************************//
int main(int argc, char *argv[]){


  MPI_Init(&argc,&argv);
  start_clock();
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  numworkers = numtasks - 1;
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

  allocate_rows();
  allocate_memory();

  if(taskid == MASTER){ // This code fragment runs in the master porcesses
    printf("This is the Master.");

    //***********main solver loop**********************************************//
      for (t = 1; t <= tsteps; t++) {
    //**********File output routine********************************************//
        if(t%savet == 0){
          //receive from worker and print
        }
        printf("iteration no. %d \n",t);
      }
  } else { // This code fragment runs in the worker processes
      printf("This is worker no. :%d with rows:%d\n", taskid, rows);

    //************main solver loop*********************************************//
      for (t = 1; t <= tsteps; t++) {

        phi_solver();
    //**********File output routine**************************************************//
        if(t%savet == 0){
          //send to master
        }
      }
  }
  free_mem();
  end_clock();
  MPI_Finalize();
  return(0);
}
void allocate_rows() {
  if (taskid == MASTER) {
    averow    =   Mx/numworkers;
    extra     =   Mx%numworkers;
    for ( rank=1; rank <= (numworkers); rank++) {
      rows      =   (rank <= extra) ? averow+1 : averow;
      dest      =   rank;
      MPI_Send(&rows, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
    }
  }
  else{
    source =  MASTER;
    MPI_Recv(&rows, 1, MPI_INT, source, BEGIN, MPI_COMM_WORLD, &status);
    start = 1;
    if((taskid ==1) || (taskid == numworkers)) {
      end = rows - 1;
    } else {
      end = rows;
    }
  }
}
void allocate_memory(){

  if(taskid == MASTER){
    phi_old   = (double *)malloc(Mx*My*sizeof(double));
    mu_old    = (double *)malloc(Mx*My*sizeof(double));
    conc      = (double *)malloc(Mx*My*sizeof(double));
  }
  else if((taskid ==1) || (taskid == numworkers)) {
    phi_old   = (double *)malloc((rows+1)*My*sizeof(double));
    phi_new   = (double *)malloc((rows+1)*My*sizeof(double));
    mu_old    = (double *)malloc((rows+1)*My*sizeof(double));
    mu_new    = (double *)malloc((rows+1)*My*sizeof(double));
    conc      = (double *)malloc((rows+1)*My*sizeof(double));
    lap_phi   = (double *)malloc((rows+1)*My*sizeof(double));
    lap_mu    = (double *)malloc((rows+1)*My*sizeof(double));
    dphi_now  = (double *)malloc(Mx*4*sizeof(double));
    dphi_next = (double *)malloc(Mx*4*sizeof(double));
  } else {
    phi_old   = (double *)malloc((rows+2)*My*sizeof(double));
    phi_new   = (double *)malloc((rows+2)*My*sizeof(double));
    mu_old    = (double *)malloc((rows+2)*My*sizeof(double));
    mu_new    = (double *)malloc((rows+2)*My*sizeof(double));
    conc      = (double *)malloc((rows+2)*My*sizeof(double));
    lap_phi   = (double *)malloc((rows+2)*My*sizeof(double));
    lap_mu    = (double *)malloc((rows+2)*My*sizeof(double));
    dphi_now  = (double *)malloc(Mx*4*sizeof(double));
    dphi_next = (double *)malloc(Mx*4*sizeof(double));
  }
}
void free_mem(){
  if(taskid == MASTER) {
    free(phi_old);
    free(mu_old);
    free(conc);
  }else{
    free(phi_old);
    free(phi_new);
    free(mu_new);
    free(phi_new);
    free(lap_phi);
    free(lap_mu);
    free(dphi_now);
    free(dphi_next);
  }
}
void start_clock(){
  if(taskid == MASTER){
    start_t = clock();
  }
  printf("Starting of the program, start_t = %ld\n", start_t);
}

void end_clock(){
  if(taskid == MASTER) {
    end_t = clock();
    printf("End of the big loop, end_t = %ld\n", end_t);
    total_t = (double)((end_t - start_t) / CLOCKS_PER_SEC);
    printf("Total time taken by CPU: %f\n", total_t  );
    printf("Exiting the program...\n");
  }
}
void sendtomaster(int taskid, int mesh, double *P) {
  dest = MASTER;

  MPI_Send(&offset,       1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&rows,         1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&left_node,    1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  MPI_Send(&right_node,   1,    MPI_INT,         dest, WRITE, MPI_COMM_WORLD);
  if (taskid == 1) {
    MPI_Send(&P[0],     rows*mesh, MPI_DOUBLE,    dest, WRITE, MPI_COMM_WORLD);
  } else {
    MPI_Send(&P[mesh], rows*mesh, MPI_DOUBLE,     dest, WRITE, MPI_COMM_WORLD);
  }
}
void receivefrmworker() {
  int rank;
  for (rank=1; rank <= numworkers; rank++) {
    source = rank;
    MPI_Recv(&offset,             1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&rows,               1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&left_node,          1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&right_node,         1,             MPI_INT,       source,   WRITE,  MPI_COMM_WORLD, &status);
    MPI_Recv(&P[offset*pmesh],    rows*pmesh,    MPI_DOUBLE,    source,   WRITE,  MPI_COMM_WORLD, &status);
  }
}
void mpiexchange(int taskid) {
  if ((taskid%2) == 0) {
    if (taskid != (numworkers)) {
      MPI_Send(&P[end*pmesh],       pmesh, MPI_DOUBLE,  right_node, LTAG,      MPI_COMM_WORLD);
      source  = right_node;
      msgtype = RTAG;
      MPI_Recv(&P[(end+1)*pmesh],  pmesh, MPI_DOUBLE,   source,     msgtype,   MPI_COMM_WORLD, &status);
    }
    MPI_Send(&P[start*pmesh],      pmesh, MPI_DOUBLE,   left_node,   RTAG,     MPI_COMM_WORLD);
    source  = left_node;
    msgtype = LTAG;
    MPI_Recv(&P[0],                pmesh, MPI_DOUBLE,   source,     msgtype,   MPI_COMM_WORLD, &status);
  } else {
    if (taskid != 1) {
       source  = left_node;
       msgtype = LTAG;
       MPI_Recv(&P[0],             pmesh, MPI_DOUBLE,   source,      msgtype,  MPI_COMM_WORLD, &status);
       MPI_Send(&P[start*pmesh],   pmesh, MPI_DOUBLE,   left_node,   RTAG,     MPI_COMM_WORLD);
    }
    if (taskid != numworkers) {
      source  = right_node;
      msgtype = RTAG;
      MPI_Recv(&P[(end+1)*pmesh],  pmesh, MPI_DOUBLE, source,      msgtype,  MPI_COMM_WORLD, &status);
      MPI_Send(&P[(end)*pmesh],    pmesh, MPI_DOUBLE, right_node,  LTAG,     MPI_COMM_WORLD);
    }
  }
}
