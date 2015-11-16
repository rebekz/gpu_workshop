#include <stdio.h>

__device__ int getGlobalIdx()
{
	return blockIdx.x * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;
}

__global__ void kernel(int * d_in)
{
	int global_idx = getGlobalIdx();
	printf("Hello world! I'm a thread %d in block %d, my global id is %d and my value is %d\n", threadIdx.x, blockIdx.x, global_idx, d_in[global_idx]);
}

int main(int argc, char** argv)
{
	int array_size = 8;
	int array_bytes = array_size * sizeof(float);
        
	// generate the input array on the host
	int h_in[array_size];
        for(int i = 0; i < array_size; i++)
	{
	  h_in[i] = i;
        }

        //declare GPU memory pointers
        int * d_in;
        
        //allocate GPU memory
	cudaMalloc((void**) &d_in, array_bytes);

        //transfer the array to the GPU
        cudaMemcpy(d_in, h_in, array_bytes, cudaMemcpyHostToDevice);
	
        //launch the kernel
        kernel<<<1, array_size>>>(d_in);
        
        //force the printf()s to flush
        cudaDeviceSynchronize();
        
        cudaFree(d_in);
       
        return 0;        
}

