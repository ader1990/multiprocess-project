#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>


#define SIZE 1024	/* assumption: SIZE a multiple of number of nodes */
#define MAXNUM 15	
#define DIFFLIMIT 0.00001*SIZE
#define FROM_MASTER 1	/* setting a message type */
#define FROM_WORKER 2	/* setting a message type */
#define DEBUG 0		/* 1 = debug on, 0 = debug off */

MPI_Status status;
static double m[SIZE+2][SIZE+2];


static void
print_matrix(void)
{
    int i, j;

    for (i = 0; i < SIZE+2; i++) {
        for (j = 0; j < SIZE+2; j++)
	    printf(" %7.2f", m[i][j]);
	printf("\n");
    }
	printf("\n");
}


static void
init_matrix(void)
{
    int i, j;

    for (i = 1; i < SIZE+1; i++) {
        for (j = 1; j < SIZE+1; j++) {
			m[i][j] = (rand() % MAXNUM) + 1.0;
        }
	}
	
    m[0][0] = m[1][1];
    m[0][SIZE+1] = m[1][SIZE];
    m[SIZE+1][0] = m[SIZE][1];
    m[SIZE+1][SIZE+1] = m[SIZE][SIZE];
    /* fix the top and bottom rows */
    for (i = 1; i < SIZE+1; i++) {
		m[0][i] = m[1][i];
		m[SIZE+1][i] = m[SIZE][i];
    }
    /* fix the left and right columns */
    for (i = 1; i < SIZE+1; i++) {
		m[i][0] = m[i][1];
		m[i][SIZE+1] = m[i][SIZE];
    }

    printf("done \n\n");
	//print_matrix();
}

int
main(int argc, char **argv)
{
    int myrank, nproc;
    int rows; /* amount of work per node (rows per worker) */
    int mtype; /* message type: send/recv between master and workers */
    int dest, src, offsetrow;
    double start_time, end_time;
    int i, j, k;
	double largestDiff = -1.00;
	double largestDiffWorker = -1;
	int iteration = 0;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	
	printf("node = %d \n",myrank);
	
	offsetrow = rows;
	do {
		if (myrank == 0) {
		
			iteration++;
			
			/* Master task */
			if (iteration == 1) {
				printf("SIZE = %d, number of nodes = %d\n", SIZE, nproc);
				init_matrix();
				start_time = MPI_Wtime();
			}
			rows = SIZE / nproc;
			/* Send part of matrix a and the whole matrix b to workers */
			mtype = FROM_MASTER;
			offsetrow = rows;
			
			/* let master do its part of the work */

			for (dest = 1; dest < nproc; dest++) {
				if (DEBUG)
					printf("   sending %d rows to task %d\n",rows,dest);
				MPI_Send(&offsetrow, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
				MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
				MPI_Send(&m[offsetrow][0], (rows+2)*(SIZE+2), MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);

				offsetrow += rows;
			}
			
			for (i = 1; i < rows+1; i++) {

				for (j = 1; j < SIZE+1; j++) {
					double tmp = (m[i-1][j]+m[i+1][j]+m[i][j-1]+m[i][j+1])/4;

					if (i == 1 && j == 1) {
						largestDiffWorker = fabs(tmp - m[i][j]);
					}
					if (fabs(tmp - m[i][j])>largestDiffWorker) {
						largestDiffWorker = fabs(tmp - m[i][j]);
					}
					
					m[i][j] = tmp;
				}
			}
			
			
			/* collect the results from all the workers */
			mtype = FROM_WORKER;
			for (src = 1; src < nproc; src++) {
				MPI_Recv(&offsetrow, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
				MPI_Recv(&rows, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
				MPI_Recv(&m[offsetrow+1][0], rows*(SIZE+2), MPI_DOUBLE, src, mtype, MPI_COMM_WORLD, &status);
				if (DEBUG)
				printf("   recvd %d rows from task %d, offset = %d\n",
					   rows, src, offsetrow);
			}
			
			if (iteration%100 == 0) {
				printf("____iteration=%d________largestdiff=%7.5f\n",iteration,largestDiff);
			}
			//print_matrix();
		} else {
		/* Worker tasks */
		
			/* Receive data from master */
			mtype = FROM_MASTER;
			MPI_Recv(&offsetrow, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&m[offsetrow][0],(rows+2)*(SIZE+2), MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &status);
			if (DEBUG) 
				printf ("Rank=%d, offset=%d, row =%d, m[offset][0]=%7.2f\n ",
					myrank, offsetrow, rows, m[offsetrow][0]);
			
			
			/* do the workers part of the calculation */
			for (i=offsetrow+1; i<offsetrow+rows+1; i++) {

				for (j=1; j<SIZE+1; j++) {

					double tmp = (m[i-1][j]+m[i+1][j]+m[i][j-1]+m[i][j+1])/4;
					
					if (i == offsetrow+1 && j == 1) {
						largestDiffWorker = fabs(tmp - m[i][j]);
					}
					if (fabs(tmp - m[i][j])>largestDiffWorker) {
						largestDiffWorker = fabs(tmp - m[i][j]);
					}
					m[i][j] = tmp;

				}
			}
			if (DEBUG)
				printf ("Rank=%d, offset=%d, row =%d, c[offset][0]=%e, largestDiffWorker=%7.2f\n",
					myrank, offsetrow, rows, m[offsetrow][0], largestDiffWorker);
			
			/* send the results to the master */
			mtype = FROM_WORKER;
			MPI_Send(&offsetrow, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
			MPI_Send(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
			MPI_Send(&m[offsetrow+1][0], rows*(SIZE+2),  MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD);	
		}

		MPI_Bcast(&iteration,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Allreduce(&largestDiffWorker,&largestDiff,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);



	}while(largestDiff > DIFFLIMIT && iteration < 10000);
	if (myrank == 0) {
		end_time = MPI_Wtime();
		printf("Execution time on %2d nodes: %f\n", nproc, end_time-start_time);
		//print_matrix();
	}
	if (iteration) {
		printf("___-_iteration=%d____-____largestdiff=%7.2f\n",iteration,largestDiff);
	}
	MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
	//print_matrix();
    return 0;
}
