#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <fcntl.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>
#include <starpu.h>



void degree_compute (void *buffers[], void *cl_arg){
	int* i=cl_arg;
	int n = STARPU_MATRIX_GET_NX(buffers[0]);
	int *mat = (int*)STARPU_MATRIX_GET_PTR(buffers[0]);
	int degree=0;
	for (int a=0; a<n; a++){
		degree+=mat[*i*n+a];
	}
	mat[*i*n]=degree;
	printf("for %d, degree is %d\n",*i, degree);
}

struct starpu_codelet cl = {
	.cpu_funcs = {degree_compute},
#ifdef STARPU_USE_CUDA
	.cuda_funcs = {degree_cuda_func},
#endif
	.nbuffers=1,
	.modes = {STARPU_RW}
};

typedef struct graph {
	int* matrix;
	int* pos;
	int size;
} graph;


void generateEdges (graph g, int threshold){
	for (int i=0; i<g.size*g.size;i++)
		g.matrix[i]=0;
	for (int i = 0; i<g.size; i++) {
		for (int j= 0; j<i; j++) {
			if ((g.pos[3*i+1]-g.pos[3*j+1])*(g.pos[3*i+1]-g.pos[3*j+1])+(g.pos[3*i+2]-g.pos[3*j+2])*(g.pos[3*i+2]-g.pos[3*j+2])<=(threshold*threshold)){
				g.matrix[g.size*i+j]=1;
				g.matrix[g.size*j+i]=1;
			}
			else {
				g.matrix[g.size*i+j]=0;		
				g.matrix[g.size*i+j]=0;			 
			}		
		}
		g.matrix[g.size*i+i]=0;	 
	}
}


int computeSizeFromFile(char* path){
	FILE *f;
	char buffer[80];

	if ((f=fopen(path,"r"))==NULL) {
		fprintf(stderr,"Error fopen");
		perror("fopen");
	}

	int i = 0;

	while (fgets(buffer,80,f)!=NULL){
		i++;
	}

	fclose(f);
	return i;
}


graph loadgraph (char* path, int threshold) {
	FILE *f;
	char buffer[21];
	char* subtoken;

	graph g;
	g.size = computeSizeFromFile(path);	
	if ((f=fopen(path,"r")) == NULL) {
		fprintf(stderr,"Error fopen");
		perror("fopen");
	}
	g.pos = (int*)malloc(3*g.size*sizeof(int));
	g.matrix = (int*)malloc(g.size*g.size*sizeof(int));


	for (int i=0; i<g.size; i++){
		fgets(buffer,80,f);
		int a = atoi(subtoken=strtok(buffer," "));
		int b = atoi(subtoken=strtok(NULL," "));
		g.pos[3*i]=i;	
		g.pos[3*i+1]=a;
		g.pos[3*i+2]=b;
	}

	fclose(f);
	generateEdges(g,threshold);	

	return g;
}

int main (int argc, char** argv){
	if (argc!=2)
		return EXIT_FAILURE;
	graph g = loadgraph(argv[1], 180);
	printf("g.size=%d\n", g.size);
	for (int i=0; i<g.size; i++){
		for (int j=0; j<g.size;j++){
			printf ("%d ",g.matrix[i*g.size+j]);
		}
		printf("\n");
	}


	starpu_init(NULL);
	starpu_data_handle_t matrix_handle;
	starpu_matrix_data_register(&matrix_handle, 0, (uintptr_t)g.matrix, g.size, g.size, g.size, sizeof(int));
	int handler[g.size];
	for (int i=0;i<g.size;i++){
		handler[i]=i;
		struct starpu_task *task = starpu_task_create();

		task->cl=&cl;
		task->handles[0] = matrix_handle;
		task->cl_arg = &handler[i];
		task->cl_arg_size=sizeof(handler[i]);
	
		starpu_task_submit(task);	
	}
	starpu_data_unregister(matrix_handle);
	starpu_shutdown();
	int* result = (int*)malloc(g.size*g.size*sizeof(int));
	for (int i=0;i<g.size*g.size;i++)
		result[i]=0;
	for (int i=0;i<g.size;i++){
		result[g.matrix[i*g.size]*g.size+i]=1;
	}
	for (int i=0;i<g.size;i++){
		for (int j=0;j<g.size;j++){
			if (result[i*g.size+j]==1)		
				printf("%d, ", j);
		}
	}
	printf("\n");

	return EXIT_SUCCESS;
}
