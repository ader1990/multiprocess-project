
<!-- saved from url=(0551)https://files.itslearning.com/File/Download/GetFile.aspx?FileName=matmul_mpi.c&Path=P6iJSPMT%2famjsdJNqoWTD3v560KVMi6JPfkFyywxOjiqqLRcu7kqoyZHH5FKhEi8mGMZyy2Bioh1RBH%2ftZnLHgd6cOGvdzNab1hnL9Cn%2fH4sMj4OX%2bPBmPBuK8n8KnlX916Ycreaba%2fMGa588vcKNgYC7040lzo0%2bsCb7fdjuOs%3d&MimeType=c&Domain=www.itslearning.com&TimeStamp=635356895146407685&Unicode=True&Signature=axjYr9n3j2tbewWMbbg0bANqsKAg67%2fCTbxLHVzuXLfYwb%2boieTY9FzRvL4WeZ%2b45fPsZSmO8SIJrp3%2fbgFQ%2f5a%2bJwBa77Qsu3LN9d4Qu7wzaeB2PQVNi%2fNdQquxBnMwXpxMIrRbMqKbzX7t8Co3JzlAV90AbzHlml7Q%2bAPdXag%3d -->
<html><head><meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></head><body><pre style="word-wrap: break-word; white-space: pre-wrap;">/***************************************************************************
 *
 * MPI-version of row-wise Matrix-Matrix multiplication
 * 
 *             File : matmul_mpi.c
 *        Author(s) : Håkan Grahn
 *          Created : 2009-01-30
 *    Last Modified : 2009-01-30
 * Last Modified by : Håkan Grahn
 * 
 * © 2009 by Håkan Grahn, Blekinge Institute of Technology.
 * All Rights Reserved
 ***************************************************************************/

/* 
 * Compile with:
 * mpicc -o mm matmul_mpi.c 
 */

#include &lt;stdio.h&gt;
#include &lt;stdlib.h&gt;
#include &lt;mpi.h&gt;

#define SIZE 1024	/* assumption: SIZE a multiple of number of nodes */
			/* Hint: use small sizes when testing, e.g., SIZE 8 */
#define FROM_MASTER 1	/* setting a message type */
#define FROM_WORKER 2	/* setting a message type */
#define DEBUG 0		/* 1 = debug on, 0 = debug off */

MPI_Status status;

static double a[SIZE][SIZE];
static double b[SIZE][SIZE];
static double c[SIZE][SIZE];

static void
init_matrix(void)
{
    int i, j;

    for (i = 0; i &lt; SIZE; i++)
        for (j = 0; j &lt; SIZE; j++) {
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

    for (i = 0; i &lt; SIZE; i++) {
        for (j = 0; j &lt; SIZE; j++)
	    printf(" %7.2f", c[i][j]);
	printf("\n");
    }
}

int
main(int argc, char **argv)
{
    int myrank, nproc;
    int rows; /* amount of work per node (rows per worker) */
    int mtype; /* message type: send/recv between master and workers */
    int dest, src, offset;
    double start_time, end_time;
    int i, j, k;

    MPI_Init(&amp;argc, &amp;argv);
    MPI_Comm_size(MPI_COMM_WORLD, &amp;nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &amp;myrank);

    if (myrank == 0) {
	/* Master task */

	/* Initialization */
	printf("SIZE = %d, number of nodes = %d\n", SIZE, nproc);
	init_matrix();
	start_time = MPI_Wtime();

	/* Send part of matrix a and the whole matrix b to workers */
	rows = SIZE / nproc;
	mtype = FROM_MASTER;
	offset = rows;
	for (dest = 1; dest &lt; nproc; dest++) {
	    if (DEBUG)
		printf("   sending %d rows to task %d\n",rows,dest);
	    MPI_Send(&amp;offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
	    MPI_Send(&amp;rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
	    MPI_Send(&amp;a[offset][0], rows*SIZE, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
	    MPI_Send(&amp;b, SIZE*SIZE, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
	    offset += rows;
	}

	/* let master do its part of the work */
	for (i = 0; i &lt; rows; i++) {
	    for (j = 0; j &lt; SIZE; j++) {
		c[i][j] = 0.0;
		for (k = 0; k &lt; SIZE; k++)
		    c[i][j] = c[i][j] + a[i][k] * b[k][j];
	    }
	}

	/* collect the results from all the workers */
	mtype = FROM_WORKER;
	for (src = 1; src &lt; nproc; src++) {
	    MPI_Recv(&amp;offset, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &amp;status);
	    MPI_Recv(&amp;rows, 1, MPI_INT, src, mtype, MPI_COMM_WORLD, &amp;status);
	    MPI_Recv(&amp;c[offset][0], rows*SIZE, MPI_DOUBLE, src, mtype, MPI_COMM_WORLD, &amp;status);
	    if (DEBUG)
		printf("   recvd %d rows from task %d, offset = %d\n",
		       rows, src, offset);
	}

	end_time = MPI_Wtime();
	if (DEBUG)
	    /* Prints the resulting matrix c */
	    print_matrix();
	printf("Execution time on %2d nodes: %f\n", nproc, end_time-start_time);

    } else {
	/* Worker tasks */

	/* Receive data from master */
	mtype = FROM_MASTER;
	MPI_Recv(&amp;offset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &amp;status);
	MPI_Recv(&amp;rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD, &amp;status);
	MPI_Recv(&amp;a[offset][0], rows*SIZE, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &amp;status);
	MPI_Recv(&amp;b, SIZE*SIZE, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD, &amp;status);
	if (DEBUG)
	    printf ("Rank=%d, offset=%d, row =%d, a[offset][0]=%e, b[0][0]=%e\n",
		    myrank, offset, rows, a[offset][0], b[0][0]);

	/* do the workers part of the calculation */
	for (i=offset; i&lt;offset+rows; i++)
	    for (j=0; j&lt;SIZE; j++) {
		c[i][j] = 0.0;
		for (k=0; k&lt;SIZE; k++)
		    c[i][j] = c[i][j] + a[i][k] * b[k][j];
	    }
	if (DEBUG)
	    printf ("Rank=%d, offset=%d, row =%d, c[offset][0]=%e\n",
		    myrank, offset, rows, a[offset][0]);

	/* send the results to the master */
	mtype = FROM_WORKER;
	MPI_Send(&amp;offset, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
	MPI_Send(&amp;rows, 1, MPI_INT, 0, mtype, MPI_COMM_WORLD);
	MPI_Send(&amp;c[offset][0], rows*SIZE, MPI_DOUBLE, 0, mtype, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}

</pre></body></html>