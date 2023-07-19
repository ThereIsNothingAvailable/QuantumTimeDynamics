#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<string.h>
#include <time.h>
#include "f2c.h"
#include "blaswrap.h"
#include "clapack.h"


/* Complex datatype */
struct _fcomplex { float re, im; };
typedef struct _fcomplex fcomplex;

fcomplex random_complex() {
    fcomplex z;
    z.re = (float)rand()/ RAND_MAX; // Random real part between 0 and 1
    z.im = (float)rand() /RAND_MAX; // Random imaginary part between 0 and 1
    return z;
}

/* Parameters */
#define N 100
#define LDA N
#define LDVL N
#define LDVR N

int main() {
    int m = 250;
    int k = 100; // Dimensions of the original matrix
    fcomplex h[(m+1)*m]; // Original matrix (flattened)
    //double complex vals_eff[100]; // Array to store eigenvalues
    //double complex vecs_eff[100 * 100]; // Array to store eigenvectors

    // Fill the original matrix 'h' with your data
    for (int i = 0; i < (m+1); i++) {
        for (int j = 0; j < m; j++) {
            h[i*m+j] = random_complex();
        }
    }
   
    
    // Printing the upper hessenburg matrix
    /*printf("Upper hessenburg matrix is : ");
    for( int i=0;i<(m+1);i++)
    {
    	for(int j =0;j<m;j++)
    	{
    		printf("for the index %d : %.2f + %.2f", (i*m+j),h[i*m+j].re,h[i*m+j].im);
    	}
    }*/
    
    // Allocate memory for sub-matrix
    fcomplex submatrix[k * k];

    // Copy the sub-matrix from the original matrix 'h'
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            submatrix[i * k + j] = h[i *m + j];
        }
    }
    

    // Variables to store eigenvalues and eigenvectors

    // LAPACK function to compute eigenvalues and eigenvectors
    char jobvl = 'N'; // Compute left eigenvectors (N for no, V for yes)
    char jobvr = 'V'; // Compute right eigenvectors (N for no, V for yes)
   
    int n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR, info, lwork;
    fcomplex wkopt;
    fcomplex* work;
    float rwork[2*N];
    fcomplex w[N], vl[LDVL*N], vr[LDVR*N];

    lwork = -1;
    cgeev_( "Vectors", "Vectors", &n, submatrix, &lda, w, vl, &ldvl, vr, &ldvr, &wkopt, &lwork, rwork, &info );

    lwork = (int)wkopt.re;
    work = (fcomplex*)malloc( lwork*sizeof(fcomplex) );
    /* Solve eigenproblem */
    cgeev_( "Vectors", "Vectors", &n, submatrix, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info );
    
    /* Check for convergence */
    if( info > 0 ) {
        printf( "The algorithm failed to compute eigenvalues.\n" );
        exit( 1 );
    }

    // Print the eigenvalues and eigenvectors
    printf("\nEigenvalues:\n");
    for (int i = 0; i < k; i++) {
        printf("%.14f + %.14fi\n", w[i].re, w[i].im);
    }

    // Print the left eigenvectors (if computed)
	if (jobvl == 'V') {
	    printf("\nLeft Eigenvectors:\n");
	    for (int i = 0; i < k; i++) {
		for (int j = 0; j <k; j++) {
		    printf("%.2f + %.2fi\t", vl[i*k + j].re, vl[i*k + j].im);
		}
		printf("\n");
	    }
	}

	// Print the right eigenvectors (if computed)
	if (jobvr == 'V') {
	    printf("\nRight Eigenvectors:\n");
	    for (int i = 0; i < k; i++) {
		for (int j = 0; j < k; j++) {
		    printf("%.2f + %.2fi\t", vr[i*k+j].re, vr[i*k+j].im);
		}
		printf("\n");
	    }
	}
   
    return 0;
}

