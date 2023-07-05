#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>


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
            result[i*size1+j] = 0.0;
            for (int k = 0; k < size1; k++) {
                result[i*size1 + j] += matrix1[i*size1 + k] * matrix2[k*size1 + j];
            }
        }
    }
    
    return result;
}	

// Define the tensor product function
float complex* tensor_product(float complex* op1, int size1, float complex* op2, int size2, float complex* result) {
    int total_size = size1 * size2;
    
    for (int i = 0; i < size1; i++) {
        for (int j = 0; j < size2; j++) {
            result[i*size2 + j] = op1[i] * op2[j];
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

int main() {
    int cutoff = 8;
    int size = cutoff*cutoff; // Specify the size of the array
    float complex a[size*size];
    float complex b[size*size];
    float complex a_dag[size*size];
    float complex b_dag[size*size];
    float complex H[size*size];
    float complex c_ops[size*size];
    int k=0;
    int delta = 5;
    int gamma=1;
    float F=4.5;
    int U=20;
    int J=10;
   
    
    // Allocate memory for operators and result
    float complex* destroy_op = (float complex*)malloc(cutoff * cutoff * sizeof(float complex));
    float complex* identity_op = (float complex*)malloc(cutoff * cutoff * sizeof(float complex));
    float complex* tensor_result = (float complex*)malloc(cutoff * cutoff * cutoff * cutoff * sizeof(float complex));
    float complex* dagres = (float complex*)malloc(cutoff * cutoff * cutoff * cutoff * sizeof(float complex));
    float complex* matrixmul = (float complex*)malloc(cutoff * cutoff * cutoff * cutoff * sizeof(float complex));
   
    
    // Generate destroy operator
    destroy_operator(destroy_op, cutoff);
    
    // Generate identity operator
    identity_operator(identity_op, cutoff);
    
    // Compute the tensor product for a
    float complex* resa= tensor_product(destroy_op, size, identity_op, size, tensor_result);
    
    // Compute the tensor product for b
    float complex* resb= tensor_product(identity_op, size,destroy_op, size, tensor_result);
    
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
           a[i*cutoff + j]=resa[i*cutoff + j];
        }}
        
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
           b[i*cutoff + j]=resb[i*cutoff + j];
        }}
        
        
   // Compute the complex conjugate transpose for both a and b      
   float complex* a_dagres = dag(a,size,dagres);
   float complex* b_dagres = dag(b, size, dagres);
   
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
           a_dag[i*cutoff + j]=a_dagres[i*cutoff + j];
        }}
        
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
           b_dag[i*cutoff + j]=b_dagres[i*cutoff + j];
        }}
        
   //Multiply matrices as required in the hamiltonian
   float complex* adag_a = matrix_mul(a_dag,size,a,size,matrixmul);
   float complex* bdag_b = matrix_mul(b_dag,size,b,size,matrixmul);
   float complex* adag_b = matrix_mul(a_dag,size,b,size,matrixmul);
   float complex* bdag_a = matrix_mul(b_dag,size,a,size,matrixmul);
   float complex* adag_sq = matrix_mul(a_dag,size,a_dag,size,matrixmul);
   float complex* bdag_sq = matrix_mul(b_dag,size,b_dag,size,matrixmul); 
   float complex* a_sq = matrix_mul(a,size,a,size,matrixmul);
   float complex* b_sq = matrix_mul(b,size,b,size,matrixmul);
   float complex* adagsq_asq = matrix_mul(adag_sq,size,a_sq,size,matrixmul);
   float complex* bdagsq_bsq = matrix_mul(bdag_sq,size,b_sq,size,matrixmul);
   
   
   // Build the hamiltonian (sample as given in jupyter notebook)
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            H[i * size + j] = -delta * (adag_a[i * size + j] + bdag_b[i * size + j]) + F * (a[i * size + j] + a_dag[i * size + j] + b[i * size + j] + b_dag[i * size + j]) + U / 2.0 * (adagsq_asq[i * size + j] + bdagsq_bsq[i * size + j]);
    }
   }
   
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            H[i * size + j] = H[i * size + j] - J*(adag_b[i * size + j]+bdag_a[i * size + j]);
    }
   }
   
   
   // Build the c_ops list
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
         	c_ops[i*size + j] = sqrt(gamma)*a[i*size + j];   
    }
   }
   
   int t=0,l=0;
   for (int i = size; i < (2*size); i++) {
        for (int j = size; j < (2*size); j++) {
         	c_ops[i*size + j] = sqrt(gamma)*b[t*size + l];
         	 l++;   
    }
    t++;
   }
   
   
   
   // Print whatever you want (timepass)😵️
    
    printf("Destroy op: ");
    for (int i = 0; i < cutoff; i++) {
        for (int j = 0; j < cutoff; j++) {
            printf("%.2f + %.2fi ", creal(destroy_op[i*cutoff + j]),cimag(destroy_op[i*cutoff + j]));
        }}
    
    
   printf("Value of tensor product: a");
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            printf("%.2f + %.2fi ", creal(a[i*cutoff + j]),cimag(a[i*cutoff + j]));
      	}
  }
   
   printf("Value of tensor product: b");
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            printf("%.2f + %.2fi ", creal(b[i*cutoff + j]),cimag(b[i*cutoff + j]));
      	}
  }
  
  printf("Value of Hamiltonian");
   for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
        	if(creal(H[i*cutoff + j])!=0)
        	{
        		printf("%.2f + %.2fi ", creal(H[i*cutoff + j]),cimag(H[i*cutoff + j]));	
        	}
      	}
  }
  
  
    
    // Free allocated memory
    free(destroy_op);
    free(identity_op);
    free(tensor_result);
    
    return 0;
}
