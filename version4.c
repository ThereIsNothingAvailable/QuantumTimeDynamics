#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include<string.h>
#include <time.h>
#include "f2c.h"
#include "blaswrap.h"
#include "clapack.h"

struct _fcomplex { float re, im; };
typedef struct _fcomplex fcomplex;

int cutoff = 2;
float complex* destroy_op; 
float complex* identity_op;
float complex* tensor_result;
float complex* tensor_result1;
float complex* dagres;
float complex* matrixmul;
float complex* li_result;
float complex** c_ops;
void identity_operator(float complex* identity_op, int cutoff);
void destroy_operator(float complex* destroy_op, int cutoff);
float complex* dag(float complex* op, int size, float complex* dag_res);
float complex* matrix_mul(float complex* matrix1, int size1, float complex* matrix2, int size2, float complex* result);
float complex* tensor_product(float complex* op1, int size1, float complex* op2, float complex* result);
void arnoldi(float complex *A, int N, float *v, int m, int n,float complex *H, float complex *Q);
float* normalize_vector(float* q, int size);
void arnoldilindbard(float complex* op,float complex** dis_ops, float* rho,int size,int n,int T,int numsteps,char* condition,int tau, int min_check, int how_often, float complex* li_result);

int main() {
    int size = cutoff*cutoff; //size of the array
    int totalsize = size*size;
    float complex cop1[size*size]; //Represents the first collapse operator
    float complex cop2[size*size]; //Represents the second collapse operator
    float complex a[size*size]; //As per the python code (to build Hamiltonian H)
    float complex b[size*size]; //As per the python code (to build Hamiltonian H)
    float complex a_dag[size*size]; //As per the python code (to build Hamiltonian H)
    float complex b_dag[size*size]; //As per the python code (to build Hamiltonian H)
    float complex adaga[size*size];
    float complex bdagb[size*size];
    float complex adagb[size*size];
    float complex bdaga[size*size];
    float complex adagsq[size*size];
    float complex bdagsq[size*size];
    float complex asq[size*size];
    float complex bsq[size*size];
    float complex adagsqasq[size*size];
    float complex bdagsqbsq[size*size];
    
    
    
    float complex H[size*size]; //Hamiltonian H
    float vec_init[size*size];
    int k=0;
    int delta = 5;
    int gamma=1;
    float F=4.5;
    int U=20;
    int J=10;
   
    
    // Allocate memory for operators and result
    destroy_op = (float complex*)malloc(cutoff * cutoff * sizeof(float complex));
    identity_op = (float complex*)malloc(cutoff * cutoff * sizeof(float complex));
    tensor_result = (float complex*)malloc(cutoff * cutoff * cutoff * cutoff * sizeof(float complex));
    dagres = (float complex*)malloc(cutoff * cutoff * cutoff * cutoff * sizeof(float complex));
    matrixmul = (float complex*)malloc(cutoff * cutoff * cutoff * cutoff  * sizeof(float complex));
    tensor_result1 = (float complex*)malloc(size * size * size * size  * sizeof(float complex));
    li_result = (float complex*)malloc(size * size * size * size * sizeof(float complex));
    
    //To create a list of arrays in c_ops (where c_ops represents a list of the collapse operators)
    int sizelist = 2;
    int size_arraylist = size*size;
    c_ops = (float complex**)malloc(sizelist * sizeof(float complex*));
	for (int i = 0; i < sizelist; i++) {
	    c_ops[i] = (float complex*)malloc(size_arraylist * sizeof(float complex));
	}
    
    // Generate destroy operator
    destroy_operator(destroy_op, cutoff);
    
    // Generate identity operator
    identity_operator(identity_op, cutoff);
    
    // Compute the tensor product for a (where cutoff represents the rows, colums for both matrices)
    float complex* resa= tensor_product(destroy_op, cutoff, identity_op, tensor_result);
   
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        			a[i*size + j] = resa[i*size + j];	
        			
       
      	}
  }
    
    // Compute the tensor product for b
    float complex* resb= tensor_product(identity_op, cutoff,destroy_op, tensor_result);
    
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        			b[i*size + j] = resb[i*size + j];	
        			
       
      	}
  }
        
 
    free(tensor_result);
   // Compute the complex conjugate transpose for both a and b      
   float complex* a_dagres = dag(a,size,dagres);
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        			a_dag[i*size + j] = a_dagres[i*size + j];	
        			
       
      	}
  }
   
   float complex* b_dagres = dag(b, size, dagres);
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        			b_dag[i*size + j] = b_dagres[i*size + j];	
        			
       
      	}
  }
  
        
   //Multiply matrices as required in the hamiltonian
   float complex* adag_a = matrix_mul(a_dag,size,a,size,matrixmul);
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        			adaga[i*size + j] = adag_a[i*size + j];	
        			
       
      	}
  }
   
   float complex* bdag_b = matrix_mul(b_dag,size,b,size,matrixmul);
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        			bdagb[i*size + j] = bdag_b[i*size + j];	
        			
       
      	}
  }
   float complex* adag_b = matrix_mul(a_dag,size,b,size,matrixmul);
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        			adagb[i*size + j] = adag_b[i*size + j];	
        			
       
      	}
  }
   float complex* bdag_a = matrix_mul(b_dag,size,a,size,matrixmul);
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        			bdaga[i*size + j] = bdag_a[i*size + j];	
        			
       
      	}
  }
   
   //Here sq is used to represent square
   float complex* adag_sq = matrix_mul(a_dag,size,a_dag,size,matrixmul);
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        			adagsq[i*size + j] = adag_sq[i*size + j];	
        			
       
      	}
  }
   float complex* bdag_sq = matrix_mul(b_dag,size,b_dag,size,matrixmul);
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        			bdagsq[i*size + j] = bdag_sq[i*size + j];	
        			
       
      	}
  } 
   float complex* a_sq = matrix_mul(a,size,a,size,matrixmul);
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        			asq[i*size + j] = a_sq[i*size + j];	
        			
       
      	}
  }
   float complex* b_sq = matrix_mul(b,size,b,size,matrixmul);
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        			bsq[i*size + j] = b_sq[i*size + j];	
        			
       
      	}
  }
   float complex* adagsq_asq = matrix_mul(adag_sq,size,a_sq,size,matrixmul);
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        			adagsqasq[i*size + j] = adagsq_asq[i*size + j];	
        			
       
      	}
  }
   float complex* bdagsq_bsq = matrix_mul(bdag_sq,size,b_sq,size,matrixmul);
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        			bdagsqbsq[i*size + j] = bdagsq_bsq[i*size + j];	
        			
       
      	}
  }
   
   
   // Build the hamiltonian (sample as given in jupyter notebook)
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            H[i * size + j] = -delta * (adaga[i * size + j] + bdagb[i * size + j]) + F * (a[i * size + j] + a_dag[i * size + j] + b[i * size + j] + b_dag[i * size + j]) + U / 2.0 * (adagsqasq[i * size + j] + bdagsqbsq[i * size + j]);
    }
   }
   
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            H[i * size + j] = H[i * size + j] - J*(adagb[i * size + j]+bdaga[i * size + j]);
    }
   }
   
   // Build the c_ops list
   for (int j = 0; j < size; j++) {
       for (int k = 0; k < size; k++)
        	{
        		c_ops[0][j*size + k] = sqrt(gamma)*a[j*size + k];  	
        	}
    }
   
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
         	c_ops[1][i*size + j] = sqrt(gamma)*b[i*size + j]; 
    }
   }
   
  
   // Populate the cop1, cop2 seperately from the c_ops list
   	for(int i=0;i<size;i++)
	{
		for(int j=0;j<size;j++){
			cop1[i*size + j] = c_ops[0][i*size + j];
		}	
	}
	
	for(int i=0;i<size;i++)
	{
		for(int j=0;j<size;j++){
			cop2[i*size + j] = c_ops[1][i*size + j];
		}	
	}
	
	
	//Create an identity operator of the required size to be able to perform tensor product with the hamiltonian(H)
	identity_operator(identity_op,size);
	
	
	//Following code represents operations to create Lioulivillian operator, here id is identity operator, c1 is cop1, c2 is cop2
	float complex hid[totalsize*totalsize];
	float complex idh[totalsize*totalsize];
	float complex c1c1dag[totalsize*totalsize];
	float complex c2c2dag[totalsize*totalsize];
	float complex idc1dagc1[totalsize*totalsize];
	float complex c1dagc1id[totalsize*totalsize];
	float complex idc2dagc2[totalsize*totalsize];
	float complex c2dagc2id[totalsize*totalsize];
	float complex c1dag[size*size];
	float complex c2dag[size*size];
	float complex c1dagc1[size*size];
	float complex c2dagc2[size*size];

	// Define the tensor products
	float complex* h_id= tensor_product(H, size, identity_op, tensor_result1);
	
	for (int i = 0; i < totalsize; i++) {
        for (int j = 0; j < totalsize; j++) {
        			hid[i*totalsize + j] = h_id[i*totalsize + j];	
        			
       
      	}
  }

	
	float complex* id_h= tensor_product(identity_op, size, H, tensor_result1);
	
	for (int i = 0; i < totalsize; i++) {
        for (int j = 0; j < totalsize; j++) {
        			idh[i*totalsize + j] = id_h[i*totalsize + j];	
        			
       
      	}
  }
	
	// Define the conjugate transpose
	float complex* c1_dag = dag(cop1,size,dagres);
	for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        			c1dag[i*size + j] = c1_dag[i*size + j];	
        			
       
      	}
  }
	
	float complex* c2_dag = dag(cop2,size,dagres);
	for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        			c2dag[i*size + j] = c2_dag[i*size + j];	
        			
       
      	}
  }
	
	// Define the tensor products
	float complex* c1_c1dag = tensor_product(cop1, size, c1dag, tensor_result1);
	for (int i = 0; i < totalsize; i++) {
        for (int j = 0; j < totalsize; j++) {
        			c1c1dag[i*totalsize + j] = c1_c1dag[i*totalsize + j];	
        			
       
      	}
  }
	float complex* c2_c2dag = tensor_product(cop2, size, c2dag, tensor_result1);
	for (int i = 0; i < totalsize; i++) {
        for (int j = 0; j < totalsize; j++) {
        			c2c2dag[i*totalsize + j] = c2_c2dag[i*totalsize + j];	
        			
       
      	}
  }
  
	
	// Define the matrix multiplication
	float complex* c1dag_c1 = matrix_mul(c1dag,size,cop1,size,matrixmul);
	for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        			c1dagc1[i*size + j] = c1dag_c1[i*size + j];	
        			
       
      	}
  }
  
	float complex* c2dag_c2 = matrix_mul(c2dag,size,cop2,size,matrixmul);
	for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        			c2dagc2[i*size + j] = c2dag_c2[i*size + j];	
        			
       
      	}}
      	
  
	
	// Define the tensor products for 1_c2dagc1, c2dagc1_1
	float complex* id_c1dagc1 = tensor_product(identity_op, size, c1dagc1, tensor_result1);
	for (int i = 0; i < totalsize; i++) {
        for (int j = 0; j < totalsize; j++) {
        			idc1dagc1[i*totalsize + j] = id_c1dagc1[i*totalsize + j];	
        			
       
      	}
  }
	float complex* c1dagc1_id = tensor_product(c1dagc1, size, identity_op, tensor_result1);
	for (int i = 0; i < totalsize; i++) {
        for (int j = 0; j < totalsize; j++) {
        			c1dagc1id[i*totalsize + j] = c1dagc1_id[i*totalsize + j];	
        			
       
      	}
  }
	float complex* id_c2dagc2 = tensor_product(identity_op, size, c2dagc2, tensor_result1);
	for (int i = 0; i < totalsize; i++) {
        for (int j = 0; j < totalsize; j++) {
        			idc2dagc2[i*totalsize + j] = id_c2dagc2[i*totalsize + j];	
        			
       
      	}
  }
	float complex* c2dagc2_id = tensor_product(c2dagc2, size, identity_op, tensor_result1);
	for (int i = 0; i < totalsize; i++) {
        for (int j = 0; j < totalsize; j++) {
        			c2dagc2id[i*totalsize + j] = c2dagc2_id[i*totalsize + j];	
        			
       
      	}
  }
	
	
	//Using equation 21 of the paper we have created the L operator in the same format, where totalsize = size*size and used to represent the entire size of the operator

	for(int i=0;i<totalsize;i++){
		for(int j=0;j<totalsize;j++){
			li_result[i*totalsize + j] = 0.5*14*(2*c1c1dag[i*totalsize + j]-idc1dagc1[i*totalsize + j]-c1dagc1id[i*totalsize + j])+ (2*c2c2dag[i*totalsize + j]-idc2dagc2[i*totalsize + j]-c2dagc2id[i*totalsize + j]);	
		}
	}
	
	for(int i=0;i<totalsize;i++){
		for(int j=0;j<totalsize;j++){
			li_result[i*totalsize + j] = li_result[i*totalsize + j] + cimag(-1)*(hid[i*totalsize + j]-idh[i*totalsize + j]);
		}	
	}
	
	
   // Create the initial density matrix
   srand(time(NULL));
    float sum = 0.0f;

    // Generate random numbers for the vector
    for (int i = 0; i < totalsize; i++) {
        vec_init[i] = (float)rand() / RAND_MAX;
    }

    for (int i = 0; i < totalsize; i++) {
        sum += vec_init[i];
    }

    // Normalize the vector by dividing each element by the sum
    for (int i = 0; i < totalsize; i++) {
        vec_init[i] /= sum;
    }
    
    //Arnoldi lindbard time evolution function definition
  
    arnoldilindbard(H,c_ops,vec_init,size,250,10,100,"steady_state",10*(-3),80,10,li_result);

	
   // Printing the liouvillian array
    
  
  /*printf("Value of liouvillian ");
   for (int i = 0; i < totalsize; i++) {
        for (int j = 0; j < totalsize; j++) { 	
        			printf("for the index %d : %.2f + %.2fi", (i*totalsize+j),creal(li_result[i*totalsize+j]),cimag(li_result[i*totalsize+j]));	
        		
       
      	}}*/

    return 0;
    }
    

// Define the arnoldi lindbard time evolution
void arnoldilindbard(float complex* op,float complex** dis_ops, float* rho,int size,int n,int T,int numsteps,char* condition,int tau, int min_check, int how_often, float complex* li_result)
{
    int k =0;
    int totalsize = size*size;
    int m = totalsize;
    int sqrtm = sqrt(m);
    // Creating an initial density matrix
    float* initial_dm = (float*)malloc(size * size * sizeof(float)); //q matrix
    float complex* h = (float complex*)malloc((n+1) * n * sizeof(float complex));
    float complex* Q_1 = (float complex*)malloc(m * (n+1) * sizeof(float complex));
    double *tlist = (double *)malloc(numsteps * sizeof(double));
    
    // Creating the h matrix as a null matrix first
    for(int i=0;i<(n+1);i++)
    {
    	for(int j=0;j<n;j++)
    	{
    		h[i*(n+1) +j] =0;
    	}
    }
    
    for(int i = 0;i<size*size;i++)
    {
    	initial_dm[i] = rho[i];
    }
    normalize_vector(initial_dm,totalsize);
    //building the first krylov subspace
    for(int i=0;i<m;i++)
    {
    	for(int j=0;j<(n+1);j++)
    	{
    		if(j == 0){
    			Q_1[i*(n+1) + j] = initial_dm[k];
    			k++;
    		}
    		else{
    			Q_1[i*(n+1)+j]=creal(0)+cimag(0);
    		}
    	}
    }
    
    
    double dt = T / (numsteps - 1);
    for (int i = 0; i < numsteps; i++) {
        tlist[i] = i * dt;
    }
    arnoldi(li_result,totalsize,initial_dm,m,n,h,Q_1);
    
    /*printf("H:\n");
    for (int i = 0; i < n+1; i++) {
        for (int j = 0; j < n; j++) {
        	if(creal(h[i*(n) + j])!=0.00f)
        		printf("for index %d: %.2f + %.2fi", (i*(n) + j),creal(h[i*(n) + j]),cimag(h[i*(n) + j]));
           
        }
    }*/
    
    factorization(h,Q_1,n,min_check,how_often,m,initial_dm);
   
    
    free(initial_dm);
    free(h);
    free(Q_1);
    free(tlist);
}

// Define the arnoldi iteration
void arnoldi(float complex *A, int N, float *v, int m, int n, float complex *H, float complex *Q) {
    float complex *v_k = (float complex *)malloc(N * sizeof(float complex));
    float complex *q_k = (float complex *)malloc(m * sizeof(float complex));
    float complex *h_k = (float complex *)malloc((n + 1) * n * sizeof(float complex));
    //float complex ritz_values[];


         for (int k = 0; k < (n + 1); k++) {
            // Extract the k-th column of Q
             for (int i = 0; i < m; i++) {
               q_k[i] = Q[i * (n + 1) + k];
		//printf("%.2f + %.2fi",creal(q_k[i]),cimag(q_k[i]));
             }
        

        // Perform matrix-vector multiplication: v_k = A * q_k
        
        for (int i = 0; i < N; i++) {
            v_k[i] = 0.0;
            for (int j = 0; j < N; j++) {
                v_k[i] = v_k[i]+ q_k[j]*A[j * N + i];
            }
        }
        
        // Orthogonalization using modified Gram-Schmidt
        for (int j = 0; j <= k; j++) {
            // Extract the j-th column of Q
            for (int i = 0; i < m; i++) {
                q_k[i] = Q[i * (n + 1) + j];
            }

            // Compute the projection: h_k[j] = q_k^T * v_k
            h_k[j * n + k] = 0.0;
            for (int i = 0; i < m; i++) {
                h_k[j * n + k] += conj(q_k[i]) * (v_k[i] + cimag(0));
            }

            // Subtract the projection: v_k = v_k - h_k[j] * q_k
            for (int i = 0; i < N; i++) {
                v_k[i] = v_k[i]- h_k[j * n + k] * q_k[i];
            }}

        // Compute the new column of Q: q_k+1 = v_k / ||v_k||
        float complex norm_v_k = creal(0) + cimag(0);
        for (int i = 0; i < N; i++) {
            norm_v_k = norm_v_k+ creal(v_k[i]) * creal(v_k[i]) + cimag(v_k[i]) * cimag(v_k[i]);
        }
        norm_v_k = csqrt(norm_v_k);

	if (norm_v_k != 0.0) {
	    for (int i = 0; i < m; i++) {
		Q[i * (n + 1) + k + 1] = v_k[i] / norm_v_k;
	    }}

        // Set the upper Hessenberg matrix H_k to the projection coefficients
        for (int j = 0; j < k; j++) {
            H[j * n + k] = h_k[j * n + k];
            if(H[j * n + k]!=0)
            {
            	//printf("%.2f + %.2fi",creal(H[j * n + k]),cimag(H[j * n + k]));	
            }
        }
        }    
}

//Define the function factorization to get eigenvalues and eigenvectors
void factorization(float complex* upper_h, float complex *krylov, int n, int min_check,int how_often,int m, float* density_init )
{
   for(int k=0;k<(n+1);k++)
   {
	 float complex* vecs_right =(float complex*)malloc(k*k*sizeof(float complex));
	 float complex* inv = (float complex*)malloc(k*k * sizeof(float complex));
	 float complex* states = (float complex*)malloc(sqrt(m)*sqrt(m) * sizeof(float complex));
	 float complex vecs_left[k*k];
	 float complex vecs_eff[k*k]; 
	 float complex vals_eff[k];
	
	if(k%how_often==0 && k>min_check){
        
          int n = k;
          fcomplex hes[k*k]; // Submatrix
          for (int i = 0; i < k; i++) {
        	for (int j = 0; j < k; j++) {
            		hes[i*k+j].re = creal(upper_h[i*k+j]);
            		hes[i*k+j].im = cimag(upper_h[i*k+j]);
        		}
    		}
    		
    		
    	// LAPACK function to compute eigenvalues and eigenvectors
    	char jobvl = 'N'; // Compute left eigenvectors (N for no, V for yes)
    	char jobvr = 'V'; // Compute right eigenvectors (N for no, V for yes)
   
    	int lda = n, ldvl = n, ldvr = n, info, lwork;
    	fcomplex wkopt;
    	fcomplex* work;
    	float rwork[2*k];
    	fcomplex w[k], vl[ldvl*k], vr[ldvr*k];

    	lwork = -1;
        cgeev_( "Vectors", "Vectors", &n, hes, &lda, w, vl, &ldvl, vr, &ldvr, &wkopt, &lwork, rwork, &info );

	lwork = (int)(wkopt.re);
	work = (fcomplex*)malloc( lwork*sizeof(fcomplex) );
	/* Solve eigenproblem */
	cgeev_( "Vectors", "Vectors", &n, hes, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info ); //where w is eigenvalues, vl is left eigenvectors, vr is right eigenvectors
	    
	/* Check for convergence */
	if( info > 0 ) {
	   printf( "The algorithm failed to compute eigenvalues.\n" );
	   exit( 1 );
	 }
	  
	  // Store the left eigenvectors (if computed)
	if (jobvl == 'V') {
	    for (int i = 0; i < k; i++) {
		for (int j = 0; j <k; j++) {
			vecs_left[i*k+j] = vl[i*k+j].re;
			vecs_left[i*k+j] = vl[i*k+j].im;
		}
	    }
	}    

	 // Store the left eigenvectors (if computed)
	if (jobvr == 'V') {
	    for (int i = 0; i < k; i++) {
		for (int j = 0; j <k; j++) {
			vecs_right[i*k+j] = vr[i*k+j].re;
			vecs_right[i*k+j] = vr[i*k+j].im;
		}
	    }
	}
	
	// Store the eigenvalues
	for (int i = 0; i < k; i++) {
		vals_eff[k] = w[k].re;	
		vals_eff[k] = w[k].im;	
		}
	 
	// Take vecs_eff as right eigenvectors for now  (from left eigenvectors and right eigenvectors (how???))
	
	
	//For loop upto k terms
	for(int jj=0; jj<k;jj++){
	   float complex vec[k*jj];
	   float complex* result = (float complex*)malloc(m*jj * sizeof(float complex));
	   if(has_converged(vals_eff[jj],vecs_right,k)==1){
	   
	   // Create a eigenvector matrix with upto jj krylov vectors
	      for(int i=0;i<k;i++){
	      	for(int j=0;j<jj;j++){
	      	    vec[i*jj + j] = vecs_right[i*jj + j];
	      	}
	      }
	   //Dot product   
	    for (int i = 0; i < m; i++) {
              for (int j = 0; j < jj; j++) {
		    result[i*jj +j] = 0.0f;
		    for (int l = 0; l < k; l++) {
		        result[i*jj + j] = result[i*jj + j]+ krylov[i*k + l] * vec[l*jj + j];
		    }
                }
            }
            
            //Convert qq matrix to square form to pass to time_evolve (as done in python)
            int qq_size=m*jj;
            normalize_vector(result,qq_size);
          
            
           //Pass the created qq matrix to time_evolve to get eigenstates and eigenvalues by using cgeev and then get sigenstates in way you understood
            time_evolve(result,sqrt(m),sqrt(m),states, vec_init, density_init);
            
            
            free(result);
	      
	   }
	}
	   	
	
	}}
}

// Define the time evolution to get evolved eigenstates of the iterative function within the factorization() 
// Only function left to complete
float complex* time_evolve(float complex* matrix, int rows,int cols,float complex* eigenstates, float* initial_rho)
{
	int n = rows;
        fcomplex qq[rows*cols]; // Submatrix
        float complex left[rows*cols];
        float complex right[rows*cols];
        float complex* left_dag = (float complex*)malloc(n*n * sizeof(float complex));
	float complex* norm = (float complex*)malloc(n*n * sizeof(float complex));
	float complex* rho_t0 = (float complex*)malloc(n*n* sizeof(float complex));
	float complex* vals = (float complex*)malloc(n*n* sizeof(float complex));
	float complex* left_dag_res = (float complex*)malloc(n*n* sizeof(float complex));
	float complex* right_res = (float complex*)malloc(n*n* sizeof(float complex));
	
	for (int i = 0; i < rows; i++) {
        	for (int j = 0; j < cols; j++) {
            		rho_t0[i*cols+j] = initial_rho[i*cols+j];
            		rho_t0[i*cols+j] = cimag(0);
        		}
    		} 
	
        for (int i = 0; i < rows; i++) {
        	for (int j = 0; j < cols; j++) {
            		qq[i*cols+j].re = creal(matrix[i*cols+j]);
            		qq[i*cols+j].im = cimag(matrix[i*cols+j]);
        		}
    		}    
    		
    		
    	// LAPACK function to compute eigenvalues and eigenvectors
    	char jobvl = 'N'; // Compute left eigenvectors (N for no, V for yes)
    	char jobvr = 'V'; // Compute right eigenvectors (N for no, V for yes)
   
    	int lda = n, ldvl = n, ldvr = n, info, lwork;
    	fcomplex wkopt;
    	fcomplex* work;
    	float rwork[2*n];
    	fcomplex w[n], vl[ldvl*n], vr[ldvr*n];

    	lwork = -1;
        cgeev_( "Vectors", "Vectors", &n, qq, &lda, w, vl, &ldvl, vr, &ldvr, &wkopt, &lwork, rwork, &info );

	lwork = (int)(wkopt.re);
	work = (fcomplex*)malloc( lwork*sizeof(fcomplex) );
	/* Solve eigenproblem */
	cgeev_( "Vectors", "Vectors", &n, qq, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info ); //where w is eigenvalues, vl is left eigenvectors, vr is right eigenvectors
	    
	/* Check for convergence */
	if( info > 0 ) {
	   printf( "The algorithm failed to compute eigenvalues.\n" );
	   exit( 1 );
	 }
	 
	  // Store the left eigenvectors
	if (jobvl == 'V') {
	    for (int i = 0; i < n; i++) {
		for (int j = 0; j <n; j++) {
			left[i*n+j] = vl[i*n+j].re;
			left[i*n+j] = vl[i*n+j].im;
		}
	    }
	}    

	 // Store the right eigenvectors
	if (jobvr == 'V') {
	    for (int i = 0; i < n; i++) {
		for (int j = 0; j <n; j++) {
			right[i*n+j] = vr[i*n+j].re;
			right[i*n+j] = vr[i*n+j].im;
		}
	    }
	}
	
	// Store eigenvalues within a diagonal matrix
	for (int i = 0; i < n; i++) {
	    for(int j=0;j<n;j++){
	        if(i==j){
	                vals[i*n + j] = w[j].re;	
		        vals[i*n + j] = w[j].im;	
	        }
	        else
	               vals[i*n + j] = 0;
	    }
	}
	
	// Compute the conjugate transpose of left and right eigenvectors
	dag(left,n,left_dag);

	// Matrix mulitplication of left_dag and right
	matrix_mul(left_dag,n,right,n,norm);
	
	
	//Normalization of left_dag and right
	for(int i =0;i<n;i++){
	   for(int j=0;j<n;j++){
	       right[i*n + j] = right[i*n + j]/sqrt(norm[i*n + j]);
	       left_dag[i*n + j] = left_dag[i*n + j]/sqrt(norm[i*n + j]); 
	   }
	}
	
	// Can take t from user also, but for now take arbitrary value
	float t = 145.56f;
	
	//Create a matrix with e^lamda(t-t0), assuming t0 = 0
	for(int i =0;i<n;i++)
	{
	   for(int j =0;j<n;j++){
	       vals[i*n+j] = exp(vals[i*n + j]*t);
	   }
	}
	
	//To create the eigenstates 
	matrix_mul(left_dag,n,rho_t0,n,left_dag_res);
	matrix_mul(right,n,vals,n,right_res);
	matrix_mul(left_dag_res,n,right_res,n,eigenstates);
	
	return eigenstates;
	
	 
	     
}


//Define the convergence as required for arnoldi eigenvalues
int has_converged(float complex val, float complex* all_vals, int size) {
    // Check if the eigenvalue 'val' has converged based on relative differences
    for (int i = 0; i < size; i++) {
        double diff = fabs((val - all_vals[i]) / all_vals[i]);
        if (diff >= 0.95) {
            return 0; // The eigenvalue has not converged
        }
    }
    return 1; // The eigenvalue has converged
}


//Define the do_lioulivillian 
float complex* do_liouvillian(float complex* rho, float complex* H, int size, float complex** cc_ops, int num_ops) {
    float complex* Lrho = (float complex*)malloc(size * size * sizeof(float complex));
    float complex* temp1 = (float complex*)malloc(size * size * sizeof(float complex));
    float complex* temp2 = (float complex*)malloc(size * size * sizeof(float complex));
    float complex* temp3 = (float complex*)malloc(size * size * sizeof(float complex));
    float complex* temp4 = (float complex*)malloc(size * size * sizeof(float complex));

    // Calculate -1.j * (HH * rho - rho * HH)
    matrix_mul(H,size, rho, size, temp1);
    matrix_mul(rho,size, H, size, temp2);
    for (int i = 0; i < size * size; i++) {
        Lrho[i] = -1.0 * cimag(1) * (temp1[i] - temp2[i]);
    }

    // Calculate Liouvillian terms for each jump operator
    for (int op = 0; op < num_ops; op++) {
        // Calculate jump * rho * jump^H
        dag(cc_ops[op],size,temp4);
        matrix_mul(rho,size,temp4,size,temp1);
        matrix_mul(cc_ops[op],size,temp1, size, temp2);
        for (int i = 0; i < size * size; i++) {
            Lrho[i] += temp2[i];
        }

        // Calculate -0.5 * (rho * jump^H * jump + jump^H * jump * rho)
        matrix_mul(temp4,size, cc_ops[op], size, temp1);
        matrix_mul(rho,size, temp1, size, temp2);
        matrix_mul(temp1,size, rho, size, temp3);
        for (int i = 0; i < size * size; i++) {
            Lrho[i] -= 0.5 * (temp2[i] + temp3[i]);
        }
    }

    free(temp1);
    free(temp2);
    free(temp3);
    free(temp4);

    return Lrho;
}



    
// Define the normalize vector function
float* normalize_vector(float* q, int size) {
    float norm = 0.0 + 0.0 * I;

    // Compute the norm of the vector
    for (int i = 0; i < size; i++) {
        norm += q[i] * q[i];
    }
    norm = sqrt(norm);

    // Normalize the vector by dividing each element by the norm
    for (int i = 0; i < size; i++) {
        q[i] /= norm;
    }
    return q;
}


// Define the conjugate transpose function
float complex* dag(float complex* op, int size, float complex* dag_res){
  for(int i=0;i<size;i++)
  {
  	for(int j=0;j<size;j++)
  	{
  		dag_res[j*size + i]=conj(op[i*size + j]); 
  	}
  }
  return dag_res;
  	
}

// Define the matrix multiplication

float complex* matrix_mul(float complex* matrix1, int size1, float complex* matrix2, int size2, float complex* result){
     for (int i = 0; i < size1; i++) {
        for (int j = 0; j < size2; j++) {
            result[i*size2+j] = 0.0;
            for (int k = 0; k < size1; k++) {
                result[i*size2 + j] += matrix1[i*size1 + k] * matrix2[k*size2 + j];
            }
        }
    }
    
    return result;
}	

// Define the tensor product function
float complex* tensor_product(float complex* op1, int cutoff, float complex* op2, float complex* result) {
    //int total_size = size1 * size2;
    
    for (int i = 0; i < cutoff; i++) {
        for (int j = 0; j < cutoff; j++) {
            for (int k = 0; k < cutoff; k++) {
                for (int l = 0; l < cutoff; l++) {
                    int index = (i * cutoff + k) * (cutoff * cutoff) + (j * cutoff + l);
                    result[index] = op1[i * cutoff + j] * op2[k * cutoff + l];
                    //printf("for the index %d : %.2f + %.2fi",index,creal(result[index]),cimag(result[index]));
                }
            }
        }
    }
    return result;
}

// Define the destroy operator
void destroy_operator(float complex* destroy_op, int cutoff) {
      for(int i = 0; i < cutoff; i++)
      {
      	for(int j = 0; j<cutoff; j++)
      	{
      		if(i==(j-1))
      		{
      		  destroy_op[i*cutoff + j]=sqrt(j);
      		}
      	}
      }
}

// Define the identity operator
void identity_operator(float complex* identity_op, int cutoff) {
    for (int i = 0; i < cutoff; i++) {
        for (int j = 0; j < cutoff; j++) {
            identity_op[i*cutoff + j] = (i == j) ? 1 : 0;
        }
    }
}

