
<!-- saved from url=(0536)https://files.itslearning.com/File/Download/GetFile.aspx?FileName=sor_seq.c&Path=S1BXclZZUCmNfJ6L5vHS9%2bb6XAaZxJeu9hf7pJVIUs8WlBmvpsojmj77p0uOxhLql6j8MVYlhMU%2ffCHrtbiHNf6Dr3gEwKwJZ1rLVuAUOUk6TT6s1EFFZn2NBC%2boVLQyE8HsnhlVcuEB12GkQQb6UIToVwyF7dNzKUaSYrW0R0E%3d&MimeType=c&Domain=www.itslearning.com&TimeStamp=635356895283726032&Unicode=True&Signature=V1DPx%2bmCiuvFKNelp5WjX0QVrreZQLLz4vv%2bIaU55FSofd0y6%2fjWiSvbqAEhLGT7d0EFt3JM3OxPs%2f1E0kjV8G5QnDy2Y4Ooxprz%2fbdn7Vzr6oROfWhHr7Y5jFvfhWkgIupTxa42wBFepehdy7R6ROJCOMlTXtc2JYCZNBrk7fc%3d -->
<html><head><meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></head><body><pre style="word-wrap: break-word; white-space: pre-wrap;">/*****************************************************
 *
 * S O R algorithm
 * ("Red-Black" solution to LaPlace approximation)
 *
 * sequential version
 *
 *****************************************************/

#include &lt;stdio.h&gt;
#include &lt;stdlib.h&gt;
#include &lt;math.h&gt;
#include &lt;malloc.h&gt;

#define MAX_SIZE 4096
#define EVEN_TURN 0 /* shall we calculate the 'red' or the 'black' elements */
#define ODD_TURN  1

typedef double matrix[MAX_SIZE+2][MAX_SIZE+2]; /* (+2) - boundary elements */

volatile struct globmem {
    int		N;		/* matrix size		*/
    int		maxnum;		/* max number of element*/
    char	*Init;		/* matrix init type	*/
    double	difflimit;	/* stop condition	*/
    double	w;		/* relaxation factor	*/
    int		PRINT;		/* print switch		*/
    matrix	A;		/* matrix A		*/
} *glob;

/* forward declarations */
int work();
void Init_Matrix();
void Print_Matrix();
void Init_Default();
int Read_Options(int, char **);

int 
main(int argc, char **argv)
{
    int i, timestart, timeend, iter;
 
    glob = (struct globmem *) malloc(sizeof(struct globmem));

    Init_Default();		/* Init default values	*/
    Read_Options(argc,argv);	/* Read arguments	*/
    Init_Matrix();		/* Init the matrix	*/
    iter = work();
    if (glob-&gt;PRINT == 1)
	Print_Matrix();
    printf("\nNumber of iterations = %d\n", iter);
}

int
work()
{
    double prevmax_even, prevmax_odd, maxi, sum, w;
    int	m, n, N, i;
    int finished = 0;
    int turn = EVEN_TURN;
    int iteration = 0;

    prevmax_even = 0.0;
    prevmax_odd = 0.0;
    N = glob-&gt;N;
    w = glob-&gt;w;
    
    while (!finished) {
	iteration++;
	if (turn == EVEN_TURN) {
	    /* CALCULATE part A - even elements */
	    for (m = 1; m &lt; N+1; m++) {
		for (n = 1; n &lt; N+1; n++) {
		    if (((m + n) % 2) == 0)
			glob-&gt;A[m][n] = (1 - w) * glob-&gt;A[m][n] 
			    + w * (glob-&gt;A[m-1][n] + glob-&gt;A[m+1][n] 
				   + glob-&gt;A[m][n-1] + glob-&gt;A[m][n+1]) / 4;
		}
	    }
	    /* Calculate the maximum sum of the elements */
	    maxi = -999999.0;
	    for (m = 1; m &lt; N+1; m++) {
		sum = 0.0;
		for (n = 1; n &lt; N+1; n++)
		    sum += glob-&gt;A[m][n];
		if (sum &gt; maxi)
		    maxi = sum;
	    }
	    /* Compare the sum with the prev sum, i.e., check wether 
	     * we are finished or not. */
	    if (fabs(maxi - prevmax_even) &lt;= glob-&gt;difflimit)
		finished = 1;
	    if ((iteration%100) == 0)
		printf("Iteration: %d, maxi = %f, prevmax_even = %f\n",
		       iteration, maxi, prevmax_even);
	    prevmax_even = maxi;
	    turn = ODD_TURN;

	} else if (turn == ODD_TURN) {
	    /* CALCULATE part B - odd elements*/
	    for (m = 1; m &lt; N+1; m++) {
		for (n = 1; n &lt; N+1; n++) {
		    if (((m + n) % 2) == 1)
			glob-&gt;A[m][n] = (1 - w) * glob-&gt;A[m][n] 
			    + w * (glob-&gt;A[m-1][n] + glob-&gt;A[m+1][n] 
				   + glob-&gt;A[m][n-1] + glob-&gt;A[m][n+1]) / 4;
		}
	    }
	    /* Calculate the maximum sum of the elements */
	    maxi = -999999.0;
	    for (m = 1; m &lt; N+1; m++) {
		sum = 0.0;
		for (n = 1; n &lt; N+1; n++)
		    sum += glob-&gt;A[m][n];	
		if (sum &gt; maxi)			
		    maxi = sum;
	    }
	    /* Compare the sum with the prev sum, i.e., check wether 
	     * we are finished or not. */
	    if (fabs(maxi - prevmax_odd) &lt;= glob-&gt;difflimit)
		finished = 1;
	    if ((iteration%100) == 0)
		printf("Iteration: %d, maxi = %f, prevmax_odd = %f\n",
		       iteration, maxi, prevmax_odd);
	    prevmax_odd = maxi;
	    turn = EVEN_TURN;
	} else {
	    /* something is very wrong... */
	    printf("PANIC: Something is really wrong!!!\n");
	    exit(-1);
	}
	if (iteration &gt; 100000) {
	    /* exit if we don't converge fast enough */
	    printf("Max number of iterations reached! Exit!\n");
	    finished = 1;
	}
    }
    return iteration;
}

/*--------------------------------------------------------------*/

void
Init_Matrix()
{
    int i, j, N, dmmy;
 
    N = glob-&gt;N;
    printf("\nsize      = %dx%d ",N,N);
    printf("\nmaxnum    = %d \n",glob-&gt;maxnum);
    printf("difflimit = %.7lf \n",glob-&gt;difflimit);
    printf("Init	  = %s \n",glob-&gt;Init);
    printf("w	  = %f \n\n",glob-&gt;w);
    printf("Initializing matrix...");
 
    /* Initialize all grid elements, including the boundary */
    for (i = 0; i &lt; glob-&gt;N+2; i++) {
	for (j = 0; j &lt; glob-&gt;N+2; j++) {
	    glob-&gt;A[i][j] = 0.0;
	}
    }
    if (strcmp(glob-&gt;Init,"count") == 0) {
	for (i = 1; i &lt; N+1; i++){
	    for (j = 1; j &lt; N+1; j++) {
		glob-&gt;A[i][j] = (double)i/2;
	    }
	}
    }
    if (strcmp(glob-&gt;Init,"rand") == 0) {
	for (i = 1; i &lt; N+1; i++){
	    for (j = 1; j &lt; N+1; j++) {
		glob-&gt;A[i][j] = (rand() % glob-&gt;maxnum) + 1.0;
	    }
	}
    }
    if (strcmp(glob-&gt;Init,"fast") == 0) {
	for (i = 1; i &lt; N+1; i++){
	    dmmy++;
	    for (j = 1; j &lt; N+1; j++) {
		dmmy++;
		if ((dmmy%2) == 0)
		    glob-&gt;A[i][j] = 1.0;
		else
		    glob-&gt;A[i][j] = 5.0;
	    }
	}
    }

    /* Set the border to the same values as the outermost rows/columns */
    /* fix the corners */
    glob-&gt;A[0][0] = glob-&gt;A[1][1];
    glob-&gt;A[0][N+1] = glob-&gt;A[1][N];
    glob-&gt;A[N+1][0] = glob-&gt;A[N][1];
    glob-&gt;A[N+1][N+1] = glob-&gt;A[N][N];
    /* fix the top and bottom rows */
    for (i = 1; i &lt; N+1; i++) {
	glob-&gt;A[0][i] = glob-&gt;A[1][i];
	glob-&gt;A[N+1][i] = glob-&gt;A[N][i];
    }
    /* fix the left and right columns */
    for (i = 1; i &lt; N+1; i++) {
	glob-&gt;A[i][0] = glob-&gt;A[i][1];
	glob-&gt;A[i][N+1] = glob-&gt;A[i][N];
    }

    printf("done \n\n");
    if (glob-&gt;PRINT == 1)
	Print_Matrix();
}

void
Print_Matrix()
{
    int i, j, N;
 
    N = glob-&gt;N;
    for (i=0; i&lt;N+2 ;i++){
	for (j=0; j&lt;N+2 ;j++){
	    printf(" %f",glob-&gt;A[i][j]);
	}
	printf("\n");
    }
    printf("\n\n");
}

void 
Init_Default()
{
    glob-&gt;N = 2048;
    glob-&gt;difflimit = 0.00001*glob-&gt;N;
    glob-&gt;Init = "rand";
    glob-&gt;maxnum = 15.0;
    glob-&gt;w = 0.5;
    glob-&gt;PRINT = 0;
}
 
int
Read_Options(int argc, char **argv)
{
    char    *prog;
 
    prog = *argv;
    while (++argv, --argc &gt; 0)
	if (**argv == '-')
	    switch ( *++*argv ) {
	    case 'n':
		--argc;
		glob-&gt;N = atoi(*++argv);
		glob-&gt;difflimit = 0.00001*glob-&gt;N;
		break;
	    case 'h':
		printf("\nHELP: try sor -u \n\n");
		exit(0);
		break;
	    case 'u':
		printf("\nUsage: sor [-n problemsize]\n");
		printf("           [-d difflimit] 0.1-0.000001 \n");
		printf("           [-D] show default values \n");
		printf("           [-h] help \n");
		printf("           [-I init_type] fast/rand/count \n");
		printf("           [-m maxnum] max random no \n");
		printf("           [-P print_switch] 0/1 \n");
		printf("           [-w relaxation_factor] 1.0-0.1 \n\n");
		exit(0);
		break;
	    case 'D':
		printf("\nDefault:  n         = %d ", glob-&gt;N);
		printf("\n          difflimit = 0.0001 ");
		printf("\n          Init      = rand" );
		printf("\n          maxnum    = 5 ");
		printf("\n          w         = 0.5 \n");
		printf("\n          P         = 0 \n\n");
		exit(0);
		break;
	    case 'I':
		--argc;
		glob-&gt;Init = *++argv;
		break;
	    case 'm':
		--argc;
		glob-&gt;maxnum = atoi(*++argv);
		break;
	    case 'd':
		--argc;
		glob-&gt;difflimit = atof(*++argv);
		break;
	    case 'w':
		--argc;
		glob-&gt;w = atof(*++argv);
		break;
	    case 'P':
		--argc;
		glob-&gt;PRINT = atoi(*++argv);
		break;
	    default:
		printf("%s: ignored option: -%s\n", prog, *argv);
		printf("HELP: try %s -u \n\n", prog);
		break;
	    } 
}
</pre></body></html>