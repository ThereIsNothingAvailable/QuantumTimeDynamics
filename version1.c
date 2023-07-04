#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Define the tensor product function
void tensor_product(int* op1, int size1, int* op2, int size2, int* result) {
    int total_size = size1 * size2;
    
    for (int i = 0; i < size1; i++) {
        for (int j = 0; j < size2; j++) {
            result[i*size2 + j] = op1[i] * op2[j];
        }
    }
}

// Define the destroy operator
void destroy_operator(int* destroy_op, int cutoff) {
      for(int i = 0; i < cutoff; i++)
      {
      	for(int j = 0; j<cutoff; j++)
      	{
      		if(i==(j-1))
      		{
      		  destroy_op[i][j]=sqrt(j);
      		}
      	}
      }
    /*for (int i = 0; i < cutoff; i++) {
    
        destroy_op[i] = i; 
    }*/
}

// Define the identity operator
void identity_operator(int* identity_op, int cutoff) {
    for (int i = 0; i < cutoff; i++) {
        for (int j = 0; j < cutoff; j++) {
            identity_op[i*cutoff + j] = (i == j) ? 1 : 0;
        }
    }
}

int main() {
    int cutoff = 8;
    int size = cutoff*cutoff; // Specify the size of the array
    int a[size];
    int k=0;
    
    // Allocate memory for the array a
    //int* a = (int*)malloc(size * sizeof(int));
    
    // Allocate memory for operators and result
    //int* destroy_op = (int*)malloc(cutoff * sizeof(int));
    int* destroy_op = (int*)malloc(cutoff * cutoff * sizeof(int));
    int* identity_op = (int*)malloc(cutoff * cutoff * sizeof(int));
    int* tensor_result = (int*)malloc(cutoff * cutoff * sizeof(int));
    
    // Generate destroy operator
    destroy_operator(destroy_op, cutoff);
    
    // Generate identity operator
    identity_operator(identity_op, cutoff);
    
    // Compute the tensor product
    tensor_product(destroy_op, cutoff, identity_op, cutoff, tensor_result);
    
    // Compute tne tensor product result in a 
    for (int i = 0; i < cutoff; i++) {
        for (int j = 0; j < cutoff; j++) {
            //a[k] = tensor_result[i*cutoff + j];
            printf("%d ", tensor_result[i*cutoff + j]);
            k++;
        }
        printf("\n");
        k++;
    }
    
    //Print the contents of a
    int arrsize = sizeof(a) / sizeof(a[0]); // Calculate the size of the array
    
    // Print the elements of the array
    /*for (int i = 0; i < arrsize; i++) {
        printf("%d ", a[i]);
    }
    printf("\n");*/
    
    // Free allocated memory
    free(destroy_op);
    free(identity_op);
    free(tensor_result);
    free(a);
    
    return 0;
}
