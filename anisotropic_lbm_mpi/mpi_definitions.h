#define MASTER  0
#define NONE    0
#define BEGIN   999
#define LTAG    777
#define RTAG    666
#define WRITE   555
#define ERROR   888
#define BREAK   111
//------------------------------------------------
int numtasks, numworkers, taskid, rank, dest;
int averow, extra, offset;
int left_node, right_node;
int start, end;
int offset_ax, offset_ay;
int source, msgtype;
int rows;
MPI_Status status;
