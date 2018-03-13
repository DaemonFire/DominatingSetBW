#ifndef _DATAREPRESENTATION_
#define _DATAREPRESENTATION_

#include <stddef.h>

/* 
Data structure for decomposition tree representation. Both void* will store pointers on sons if they exist or none if they don't exist. The label will be 0 if the node is not a terminal singleton. The use of size_t ensures that we can deal with large graphs where int would cause some issues with big graphs, as it would loop back to negative values when it has reached the maximum 
*/

typedef struct dectree {
	struct dectree *left;
	struct dectree *right;
	int label;

	int *matrixrevisited;
	int *tc;
	int nrep;
	int *complementtc;
	int nrepincomp;
	int *pointtorep;
	int *pointtorepincomp;

	int **lra;
	int **lrcompa;
	int **assoc;

} dectree; 


/*
This data structure is a representation of the graph. It is basically its adjacency matrix and x and y positions of each node. In order to be easier to use, we also add the size, for it to be more accessible
*/

typedef struct graph {
	int* matrix;
	int* pos;
	int size;
} graph;

#endif


