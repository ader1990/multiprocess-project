/***************************************************************************
 *
 * MPI-version of row-wise Matrix-Matrix multiplication
 * 
 *             File : matmul_mpi.c
 *        Author(s) : Carduner damien
 *          Created : 2014-04-24
 *    Last Modified : 2014-04-26
 * Last Modified by : Håkan Grahn
 * 
 ***************************************************************************/

/* 
 * Compile with:
 * mpicc -o mm matmul_mpi.c 
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 1024	/* assumption: SIZE a multiple of number of nodes */
			/* Hint: use small sizes when testing, e.g., SIZE 8 */
#define FROM_MASTER 1	/* setting a message type */
#define FROM_WORKER 2	/* setting a message type */
#define DEBUG 1		/* 1 = debug on, 0 = debug off */

MPI_Status status;
MPI_Datatype rowtype;
MPI_Datatype columntype;
MPI_Datatype resulttype;

static double a[SIZE][SIZE];
static double b[SIZE][SIZE];
static double c[SIZE][SIZE];

static void
init_matrix(void)
{
    int i, j;

    for (i = 0; i < SIZE; i++)
        for (j = 0; j < SIZE; j++) {
			/* Simple initialization, which enables us to easily check
			 * the correct answer. Each element in c will have the same 
			 * value as SIZE after the matmul operation.
			 */
			a[i][j] = 1.0;
			b[i][j] = 1.0;
        }
}

static void
print_matrix(void)
{
    int i, j;

    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++)
	    printf(" %7.2f", c[i][j]);
	printf("\n");
    }
	printf("\n");
}

int
main(int argc, char **argv)
{
    int myrank, nproc;
    int rows, columns; /* amount of work per node (rows per worker) */
    int mtype; /* message type: send/recv between master and workers */
    int dest, src, offsetrow, offsetcolumn;
    double start_time, end_time;
    int i, j, k;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	if (nproc == 4 || nproc == 2) {
		rows = SIZE / 2;
	} else if (nproc == 1) {
		rows = SIZE;
	}
	
	if (nproc == 4) {
		columns = SIZE / 2;
	} else if (nproc == 2 || nproc == 1) {
		columns = SIZE;
	} 
	MPI_Type_contiguous(SIZE*rows,MPI_DOUBLE,&rowtype);
	MPI_Type_commit(&rowtype);
	MPI_Type_vector(SIZE,columns,SIZE,MPI_DOUBLE,&columntype);
	MPI_Type_commit(&columntype);
	MPI_Type_vector(rows,columns,SIZE,MPI_DOUBLE,&resulttype);
	MPI_Type_commit(&resulttype);

    if (myrank == 0) {
	/* Master task */

		/* Initialization */
		// printf("SIZE = %d, number of nodes = %d\n", SIZE, nproc);
		init_matrix();
		start_time = MPI_Wtime();

				/* Send part of matrix a and the whole matrix b to workers */

		
		mtype = FROM_MASTER;
		offsetrow = 0;
		offsetcolumn = 0;
		for (dest = 1; dest < nproc; dest++) {
		
		
			if (DEBUG)
				printf("   sending %d rows and %d columns to task %d\n",rows,columns,dest);
			offsetrow = (offsetrow+rows)%SIZE;	
			if (dest == 2) {
				offsetcolumn = (offsetcolumn+columns)%SIZE;
			}	
			
			MPI_Send(&offsetrow, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
			MPI_Send(&offsetcolumn, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
			MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
			MPI_Send(&columns, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
			MPI_Send(&a[offsetrow][0], 1, rowtype, dest, mtype, MPI_COMM_WORLD);
			MPI_Send(&b[0][offsetcolumn], 1, columntype, dest, mtype, MPI_COMM_WORLD);
			

		}
		printf(" ---- SEND ----- Execution time on %2d nodes: %5.2f\n", myrank, MPI_Wtime()-start_time);
		/* let master do its part of the work */
		for (i = 0; i < rows; i++) {
			for (j = 0; j < columns; j++) {
				c[i][j] = 0;
				for (k = 0; k < SIZE; k++)
				{
					c[i][j] += a[i][k] * b[k][j];
				}
			}
		}
		printf("---- algo ----- Execution time on %2d nodes: %5.2f\n",  myrank, MPI_Wtime()-start_time);

		/* collect the results from all the workers */
		mtype = FROM_WORKER;
		for (src = 1; src < nproc; src++) {
			MPI_Recv(&offsetrow, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&offsetcolumn, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&rows, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&columns, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&c[offsetrow][offsetcolumn], 1, resulttype, src, mtype, MPI_COMM_WORLD, &status);
			if (DEBUG)
			printf("   recvd %d rows and %d columns from task %d, offsetrow = %d, offsetcolumn = %d\n",
				   rows, columns, src, offsetrow, offsetcolumn);
		}
		printf(" ---- RECV ----- Execution time on %2d nodes: %5.2f\n", myrank, MPI_Wtime()-start_time);
		end_time = MPI_Wtime();

		printf("Execution time on %2d nodes: %f\n", nproc, end_time-start_time);
		//if (DEBUG)
			/* Prints the resulting matrix c */
			//print_matrix();
    } else {
	/* Worker tasks */


		/* Receive data from master */
		mtype = FROM_MASTER;
		MPI_Recv(&offsetrow, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&offsetcolumn, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&columns, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &status);	
		
		
		MPI_Recv(&a[offsetrow][0],1, rowtype, 0, mtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&b[0][offsetcolumn], 1, columntype, 0, mtype, MPI_COMM_WORLD, &status);
		if (DEBUG)
			printf ("Rank=%d, offsetrow=%d,offsetcolumn=%d, row =%d, column=%d, a[offsetrow][0]=%e, b[0][offsetcolumn]=%e\n",
				myrank, offsetrow, offsetcolumn, rows, columns,  a[offsetrow][0], b[0][offsetcolumn]);

		/* do the workers part of the calculation */
		for (i=offsetrow; i<offsetrow+rows; i++) {
			for (j=offsetcolumn; j<offsetcolumn+columns; j++) {
				c[i][j] = 0.0;
				for (k=0; k<SIZE; k++){
					c[i][j] = c[i][j] + a[i][k] * b[k][j];
				}
			}
		}
		if (DEBUG)
			printf ("Rank=%d, offsetrow=%d,  offsetcolumn=%d,row =%d, column=%d, c[offsetrow][0]=%e\n",
				myrank, offsetrow,offsetcolumn, rows, columns, a[offsetrow][0]);

		/* send the results to the master */
		mtype = FROM_WORKER;
		MPI_Send(&offsetrow, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
		MPI_Send(&offsetcolumn, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
		MPI_Send(&rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
		MPI_Send(&columns, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
		MPI_Send(&c[offsetrow][offsetcolumn], 1, resulttype, 0, mtype, MPI_COMM_WORLD);

    }

    MPI_Finalize();
    return 0;
}
