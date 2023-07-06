#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include<string.h>

int cutoff = 8;
float complex* destroy_op; 
float complex* identity_op;
float complex* tensor_result;
float complex* tensor_result1;
float complex* dagres;
float complex* matrixmul;
float complex* li_result;
float complex** c_ops;
float complex* identity_operator(float complex* identity_op, int cutoff);
float complex* destroy_operator(float complex* destroy_op, int cutoff);
float complex* dag(float complex* op, int size, float complex* dag_res);
float complex* matrix_mul(float complex* matrix1, int size1, float complex* matrix2, int size2, float complex* result);
float complex* tensor_product(float complex* op1, int size1, float complex* op2, int size2, float complex* result);


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
float complex* destroy_operator(float complex* destroy_op, int cutoff) {
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
float complex* identity_operator(float complex* identity_op, int cutoff) {
    for (int i = 0; i < cutoff; i++) {
        for (int j = 0; j < cutoff; j++) {
            identity_op[i*cutoff + j] = (i == j) ? 1 : 0;
        }
    }
}

int main() {
    int size = cutoff*cutoff; //size of the array
    int totalsize = size*size;
    float complex cop1[size*size]; //Represents the first collapse operator
    float complex cop2[size*size]; //Represents the second collapse operator
    float complex a[size*size]; //As per the python code (to build Hamiltonian H)
    float complex b[size*size]; //As per the python code (to build Hamiltonian H)
    float complex a_dag[size*size]; //As per the python code (to build Hamiltonian H)
    float complex b_dag[size*size]; //As per the python code (to build Hamiltonian H)
    float complex H[size*size]; //Hamiltonian H
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
    
    // Compute the tensor product for a
    float complex* resa= tensor_product(destroy_op, size, identity_op, size, tensor_result);
    
    // Compute the tensor product for b
    float complex* resb= tensor_product(identity_op, size,destroy_op, size, tensor_result);
    
   //Assign the tensor_result values into corresponding a,b 
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
   
   //Here sq is used to represent square
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

	// Define the tensor products
	float complex* h_id= tensor_product(H, totalsize, identity_op, totalsize, tensor_result1);
	float complex* id_h= tensor_product(identity_op, totalsize, H, totalsize, tensor_result1);
	
	// Define the conjugate transpose
	float complex* c1dag = dag(cop1,size,dagres);
	float complex* c2dag = dag(cop2,size,dagres);
	
	// Define the tensor products
	float complex* c1_c1dag = tensor_product(cop1, totalsize, c1dag, totalsize, tensor_result1);
	float complex* c2_c2dag = tensor_product(cop2, totalsize, c2dag, totalsize, tensor_result1);
	
	// Define the matrix multiplication
	float complex* c1dagc1 = matrix_mul(c1dag,size,cop1,size,matrixmul);
	float complex* c2dagc2 = matrix_mul(c2dag,size,cop2,size,matrixmul);
	
	// Define the tensor products for 1_c2dagc1, c2dagc1_1
	float complex* id_c1dagc1 = tensor_product(identity_op, totalsize, c1dagc1, totalsize, tensor_result1);
	float complex* c1dagc1_id = tensor_product(c1dagc1, totalsize, identity_op, totalsize, tensor_result1);
	float complex* id_c2dagc2 = tensor_product(identity_op, totalsize, c2dagc2, totalsize, tensor_result1);
	float complex* c2dagc2_id = tensor_product(c2dagc2, totalsize, identity_op, totalsize, tensor_result1);
	
	
	//Using equation 21 of the paper we have created the L operator in the same format, where totalsize = size*size and used to represent the entire size of the operator

	for(int i=0;i<totalsize;i++){
		for(int j=0;j<totalsize;j++){
			li_result[i*totalsize + j] = 0.5*19*(2*c1_c1dag[i*totalsize + j]-id_c1dagc1[i*totalsize + j]-c1dagc1_id[i*totalsize + j])+ 2*c2_c2dag[i*totalsize + j]-id_c2dagc2[i*totalsize + j]-c2dagc2_id[i*totalsize + j]);	
		}
	}
	
	for(int i=0;i<totalsize;i++){
		for(int j=0;j<totalsize;j++){
			li_result[i*totalsize + j] = li_result[i*totalsize + j] + cimag(-1)*(h_id[i*totalsize + j]-id_h[i*totalsize + j]);
		}	
	}
	
   // Printing the liouvillian array
    
  
  printf("Value of liouvillian ");
   for (int i = 0; i < totalsize; i++) {
        for (int j = 0; j < totalsize; j++) {
        		if(creal(H[i*totalsize+j])!=0 || cimag(H[i*totalsize+j])!=0)
        		{
        			printf("%.2f + %.2fi ", creal(H[i*totalsize+j]),cimag(H[i*totalsize+j]));	
        		}	
       
      	}
  }
 
    return 0;
}
   
