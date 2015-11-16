#include <stdio.h>

__device__ int getGlobalIdx()
{
	return blockIdx.x * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;
}

__global__ void kernel(int * d_in, int * d_out)
{
	int global_idx = getGlobalIdx();
	printf("Hello world! I'm a thread %d in block %d, my global id is %d and my value is %d\n", threadIdx.x, blockIdx.x, global_idx, d_in[global_idx]);
	d_out[global_idx] = d_in[global_idx] * d_in[global_idx];
}

int main(int argc, char** argv)
{
	int array_size = 8;
	int array_bytes = array_size * sizeof(float);
        
	// generate the input array on the host
	int h_in[array_size];
	printf("input: \n");
        for(int i = 0; i < array_size; i++)
	{
	  h_in[i] = i;
	  printf("%d ", h_in[i]);
        }
        int h_out[array_size];
	printf("\n");        

        //declare GPU memory pointers
        int * d_in;
        int * d_out;
        
        //allocate GPU memory
	cudaMalloc((void**) &d_in, array_bytes);
	cudaMalloc((void**) &d_out, array_bytes);

        //transfer the array to the GPU
        cudaMemcpy(d_in, h_in, array_bytes, cudaMemcpyHostToDevice);
	
        //launch the kernel
        kernel<<<1, array_size>>>(d_in, d_out);
        
        //force the printf()s to flush
        cudaDeviceSynchronize();
        
	//copy back the result array to the CPU
	cudaMemcpy(h_out, d_out, array_bytes, cudaMemcpyDeviceToHost);

	//print out the resulting array
        printf("Result:\n");
	for(int i = 0; i < array_size; i++)
	{
		printf("%d ", h_out[i]);
	}
	
	printf("\n");

        cudaFree(d_in);
	cudaFree(d_out);
       
        return 0;        
}

