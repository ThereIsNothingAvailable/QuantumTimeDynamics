#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

float complex* destroy_op; 
float complex* identity_op;
float complex* tensor_result;
float complex* dagres;
float complex* matrixmul;
float complex** c_ops;
float complex* dag(float complex* op, int size, float complex* dag_res);
float complex* matrix_mul(float complex* matrix1, int size1, float complex* matrix2, int size2, float complex* result);
float complex* tensor_product(float complex* op1, int size1, float complex* op2, int size2, float complex* result);
float complex* liovillian(float complex* H, int size1, float complex** c_ops, int size2, float complex* result); 

// Define the liovillian operator
float complex* liovillian(float complex* H, int size1, float complex** c_ops, int size2, float complex* result){

	int totalsize = size1*size1;
	float complex* cop1 = (float complex*)malloc(size2 * size2* sizeof(float complex));
	float complex* cop2 = (float complex*)malloc(size2 * size2* sizeof(float complex));

	
	// Define the identity operator
	float complex identity[size1*size1];
	for(int i=0;i<size1;i++){
		for(int j=0;j<size1;j++){
			if(i==j)
			{
				identity[i*size1 + j] = 1;
			}
			else
				identity[i*size1 + j] =0;
		}
	}
	
	for(int i=0;i<size2;i++)
	{
		for(int j=0;j<size2;j++){
			cop1[i*size2 + j] = c_ops[0][i*size2 + j];
		}	
	}
	
	for(int i=0;i<size2;i++)
	{
		for(int j=0;j<size2;j++){
			cop2[i*size2 + j] = c_ops[1][i*size2 + j];
		}	
	}

	// Define the tensor products
	float complex* h_id= tensor_product(H, totalsize, identity, totalsize, tensor_result);
	float complex* id_h= tensor_product(identity, totalsize, H, totalsize, tensor_result);
	float complex* c1_c2 = tensor_product(cop1, totalsize, cop2, totalsize, tensor_result);
	
	// Define the conjugate transpose
	float complex* c2dag = dag(cop1,size2,dagres);
	
	// Define the matrix multiplication
	float complex* c2dagc1 = matrix_mul(c2dag,size2,cop1,size2,matrixmul);
	
	// Define the tensor products for 1_c2dagc1, c2dagc1_1
	float complex* id_c2dagc1 = tensor_product(identity, totalsize, c2dagc1, totalsize, tensor_result);
	float complex* c2dagc1_id = tensor_product(c2dagc1, totalsize, identity, totalsize, tensor_result);
	
	// Define the term1 for L
	float complex liouvillian_op[totalsize*totalsize];

	for(int i=0;i<totalsize;i++){
		for(int j=0;j<totalsize;j++){
			result[i*totalsize + j] = cimag(-1)*(h_id[i*totalsize + j]-id_h[i*totalsize + j]) + 0.5*5*(2*c2dagc1[i*totalsize + j]-id_c2dagc1[i*totalsize + j]-c2dagc1_id[i*totalsize + j]);	
		}
	}
	
	free(cop1);
	free(cop2);
	
	return result;
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
    float complex* li_result = (float complex*)malloc(size * size * size * size * sizeof(float complex));
    
    //To create a list of arrays in c_ops
    int sizelist = 2;
    int size_arraylist = size*size;
    float complex** c_ops = (float complex**)malloc(sizelist * sizeof(float complex*));
	for (int i = 0; i < sizelist; i++) {
	    c_ops[i] = (float complex*)malloc(size_arraylist * sizeof(float complex));
	}


   
    
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
   
   float complex* li_op = liovillian(H,size,c_ops,size,li_result);
   
   // Print whatever you want (timepass)😵️
    
    /*printf("Destroy op: ");
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
  }*/
  
  
  /*printf("Value of c_ops");
   for (int i = 0; i < sizelist; i++) {
        for (int j = 0; j < size_arraylist; j++) {
        	if(creal(c_ops[i][j])!=0)
        	{
        		printf("%.2f + %.2fi ", creal(c_ops[i][j]),cimag(c_ops[i][j]));	
        	}
      	}
  }*/
  
  printf("%.2f + %.2fi ", creal(c_ops[0][45]),cimag(c_ops[0][45]));
  
    
    // Free allocated memory
    free(destroy_op);
    free(identity_op);
    free(tensor_result);
    free(dagres);
    free(matrixmul);
    free(li_result);
    free(c_ops);
    
    return 0;
}
