#include<stdio.h>
#include <stdlib.h>
#include <time.h>


int main()
{

	int n= 750, b = 50;
	double N = n*n;
	double *A, *B, *C;
	
	/* Allocate memory for each array */

  	A = (double *)malloc(N*sizeof(double));
  	B = (double *)malloc(N*sizeof(double));
  	C = (double *)malloc(N*sizeof(double));
	
	/* see memory allocated or not */	

	if(A == NULL)                     
    	{
       		printf("Error! memory not allocated for A.");
       		exit(0);
    	}
	if(B == NULL)                     
    	{
        	printf("Error! memory not allocated B.");
        	exit(0);
    	}
	if(C == NULL)                     
    	{
        	printf("Error! memory not allocated for C.");
        	exit(0);
    	}

	/* Initialize array1 */
  	
	for (int indx = 0; indx < N; indx++)
  	{
    		A[indx] = 2;
	
  	}
	printf("A initialized \n");
	/* Initialize array2 */
	for (int indx = 0; indx < N; indx++)
        {
                B[indx] = 3;
	
        }
	printf("B initialized \n");
	/* Initialize array2 */
        for (int indx = 0; indx < N; indx++)
        {
                C[indx] = 0;

        }
	printf("C initialized \n");
	
	/* Matrix multiplication*/
	
	clock_t begin = clock();

	for(int i = 0; i < n; i+=b)
        {
                for (int j = 0; j < n; j+=b)
                {
                        for (int k = 0;k < n; k+=b)
                        {
                                /* b x b mini matrix multiplications */
				for (int i1 = i; i1 < i+b; i1++)
				{
					for (int j1 = j; j1 < j+b; j1++)
					{
						double cij = C[i1+j1*n];
						for (int k1 = k; k1 < k+b; k1++)
						{
							cij += A[i1 + k1*n]*B[k1 + j1*n];
							C[i1+j1*n] = cij; 
						}
					}
				}
                        }
                }
        }
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	
	printf("multiplication done \n");
	printf("4th elem of c is :%lf ", C[3]);
	printf("time taken is: %lf \n", time_spent);
	free(A);
  	free(B);
  	free(C);

  	return 0;
}
