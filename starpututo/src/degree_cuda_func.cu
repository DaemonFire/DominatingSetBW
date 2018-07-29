#include <starpu.h>

static void degree_compute_cuda(int n, int i, int *mat){
	unsigned j= blockIdx.x*blockDim.x + threadIdx.x;
	if ((j<n)&&(j>0))
		mat[i*n]+=mat[i*n+j];
}

extern "C" void degree_cuda_func (void *buffers[], void *_args){
	int n = STARPU_MATRIX_GET_NX(buffers[0]);
	int* i=(int*)_args;
	int *mat = (int*)STARPU_MATRIX_GET_PTR(buffers[0]);
	unsigned threads_per_block = 64;
	unsigned nblocks = (n + threads_per_block-1) / threads_per_block;

	degree_compute_cuda<<<nblocks, threads_per_block, 0, starpu_cuda_get_local_stream()>>>(n,*i,mat);
	cudaStreamSynchronize(starpu_cuda_get_local_stream());
}
