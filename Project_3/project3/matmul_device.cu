#include <cassert>
#include <cuda_runtime.h>
#include "matmul_device.cuh"

/*
 * Read TODO items below
 */




__global__
void naiveMatmul(float *a, float *b, float *c, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    float acc = 0;
    for (int k=0; k<n; k++) {
	acc += a[i*n+k] * b[k*n+j];
    }
    c[i*n+j] = acc;
}



__global__ void cacheMatmul(float *a, float *b, float *c, int n)
{

     int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    float acc = 0;
    for(int k1=0;k1<n;k1+=gridDim.x)
    {
    acc=c[i*n+j];
    for(int k=k1;k<k1+gridDim.x;k++)
     {
    	acc += a[i*n+k] * b[k*n+j];
     }
      c[i*n+j] = acc;
     }
}


__global__ void sharedMatmul(float *a, float *b, float *c, int n)
{

    __shared__ float A_tile[32][32];
    __shared__ float B_tile[32][32];
    int width = gridDim.x*blockDim.x;

    float acc = 0;   
    
    int i = blockIdx.x*32 + threadIdx.x;
    int j = blockIdx.y*32 + threadIdx.y;
    

    /* Accumulate C tile by tile. */

    for (int tileIdx = 0; tileIdx < gridDim.x ; tileIdx+=1)
    {

        /* Load one tile of A and one tile of B into shared mem */
    
	A_tile[threadIdx.y][ threadIdx.x] = a[j * width + tileIdx*32+threadIdx.x];  
        B_tile[threadIdx.y][threadIdx.x] = b[(tileIdx * 32 + threadIdx.y)* width+ i ]; 
    
        __syncthreads();                                        

        /* Accumulate one tile of C from tiles of A and B in shared mem */

        for (int k = 0 ;k < 32; k++)
        {   
            acc += A_tile[threadIdx.y][k] * B_tile[k][threadIdx.x];    
        }
    
        __syncthreads();                                                            

    }

    c[j * width + i ] = acc;                            
    
}


void cudaMatmul(float *a, float *b, float *c, int n, MatmulImplementation type)
{
    // TODO: play with the gridSize and blockSize to find the best one
    if (type == NAIVE) {
        dim3 blockSize(32, 32);
        dim3 gridSize(n / 32, n / 32);
        naiveMatmul<<<gridSize, blockSize>>>(a,b,c,n);
    }
    else if (type == CACHE) {
        dim3 blockSize(32, 32);
        dim3 gridSize(n / 32, n / 32);
        cacheMatmul<<<gridSize, blockSize>>>(a,b,c,n);
    }
    else if (type == SHARED) {
        dim3 blockSize(32, 32);
        dim3 gridSize(n / 32, n / 32);
        sharedMatmul<<<gridSize, blockSize>>>(a,b,c,n);
    }
    // Unknown type
    else
        assert(false);
}
