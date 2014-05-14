
<!-- saved from url=(0539)https://files.itslearning.com/File/Download/GetFile.aspx?FileName=matmul_seq.c&Path=lVmIj28ENP9a0WWaWPMvaeJL2%2fvkFpZ4IEfp1Fn3Jzy%2bRkLJY8gSb4QKa5miL9t6ygn0Cnelap7sYjpUdnFU%2bB0m2cB1DmZCBIyne6Q99aoI9i1mczBYBDpmybzQz0Rav6ZGTP0CscIMgyp9G4dntHs6CXMfpTzENi8kgyZvcM4%3d&MimeType=c&Domain=www.itslearning.com&TimeStamp=635356894574636863&Unicode=True&Signature=O2RJASRyOqfKulUGuBGdUKtWvPqZeybjvW%2ba2zmn3G0Q0fMggBO1QiDvXFNqFNOfhvcTx5muEuT%2fAV6up0EHjKjJNCNgZdXUzUSKmlNI1qlELgtutj9l%2f5il17ob3K%2fUzPd%2fqYNPNllAuNBFLzC62WpcTkkRReRemexxEUSsST0%3d -->
<html><head><meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"></head><body><pre style="word-wrap: break-word; white-space: pre-wrap;">/***************************************************************************
 *
 * Sequential version of Matrix-Matrix multiplication
 *
 ***************************************************************************/

#include &lt;stdio.h&gt;
#include &lt;stdlib.h&gt;

#define SIZE 1024

static double a[SIZE][SIZE];
static double b[SIZE][SIZE];
static double c[SIZE][SIZE];

static void
init_matrix(void)
{
    int i, j;

    for (i = 0; i &lt; SIZE; i++)
        for (j = 0; j &lt; SIZE; j++) {
	    /* Simple initialization, which enables us to easy check
	     * the correct answer. Each element in c will have the same 
	     * value as SIZE after the matmul operation.
	     */
	    a[i][j] = 1.0;
	    b[i][j] = 1.0;
        }
}

static void
matmul_seq()
{
    int i, j, k;

    for (i = 0; i &lt; SIZE; i++) {
        for (j = 0; j &lt; SIZE; j++) {
            c[i][j] = 0.0;
            for (k = 0; k &lt; SIZE; k++)
                c[i][j] = c[i][j] + a[i][k] * b[k][j];
        }
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
    init_matrix();
    matmul_seq();
    //print_matrix();
}

</pre></body></html>