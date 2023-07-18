#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include<string.h>
#include <time.h>
#include "f2c.h"
#include "clapack.h"
#include "blaswrap.h"

#define MAX_SIZE 100

int main() {
    int n = 250;
    int k = 100; // Dimensions of the original matrix
    double complex h[251*250]; // Original matrix (flattened)
    double complex vals_eff[100]; // Array to store eigenvalues
    double complex vecs_eff[100 * 100]; // Array to store eigenvectors
    int sum =0;

    // Fill the original matrix 'h' with your data
    // Fill the original matrix 'h' with your data
    for(int i=0;i<(n+1);i++)
    {
    	for(int j=0;j<n;j++)
    	{
    		h[i*n + j] = (double)rand() / RAND_MAX;	
    	}
    }
    
    for(int i=0;i<(n+1);i++)
    {
    	for(int j=0;j<n;j++)
    	{
    		sum = sum + h[i*n + j];	
    	}
    }

    // Normalize the vector by dividing each element by the sum
    for(int i=0;i<(n+1);i++)
    {
    	for(int j=0;j<n;j++)
    	{
    		h[i*n + j] = sum/h[i*n + j];	
    	}
    }

    // Compute the effective eigenvalues and eigenvectors
    int info;
    int lda = k; // Leading dimension of the sub-matrix

    // Allocate memory for sub-matrix
    double complex submatrix[MAX_SIZE * MAX_SIZE];

    // Copy the sub-matrix from the original matrix 'h'
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            submatrix[i * k + j] = h[i * (n + 1) + j];
        }
    }

    char jobvl = 'N'; // Do not compute left eigenvectors
    char jobvr = 'V'; // Compute right eigenvectors

    // Temporary workspace and sizes
    double complex work_query;
    doublereal rwork;
    int lwork = -1;
    int ldvr =-1;
    int ldvl = -1;

    // Query the workspace size;
    zgeev_(&jobvl, &jobvr, &k, submatrix, &lda, vals_eff, vecs_eff, &ldvl, vecs_eff, &ldvr, &work_query, &lwork, &rwork, &info);


    if (info != 0) {
        printf("Error: zgeev_ failed with error code %d\n", info);
        return info;
    }

    // Allocate memory for workspace
    lwork = (int)creal(work_query);
    double complex* work = (double complex*)malloc(lwork * sizeof(double complex));

    // Compute eigenvalues and eigenvectors with allocated workspace
      zgeev_(&jobvl, &jobvr, &k, submatrix, &lda, vals_eff, vecs_eff, &ldvl, vecs_eff, &ldvr, &work_query, &lwork, &rwork, &info);


    if (info != 0) {
        printf("Error: zgeev_ failed with error code %d\n", info);
        free(work);
        return info;
    }

    // Print the effective eigenvalues
    printf("Effective Eigenvalues:\n");
    for (int i = 0; i < k; i++) {
        printf("%f + %fi\n", creal(vals_eff[i]), cimag(vals_eff[i]));
    }
    printf("\n");

    // Print the effective eigenvectors
    printf("Effective Eigenvectors:\n");
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            printf("%f + %fi\t", creal(vecs_eff[i * k + j]), cimag(vecs_eff[i * k + j]));
        }
        printf("\n");
    }

    // Free allocated memory
    free(work);

    return 0;
}
