#ifndef _DATAREPRESENTATION_
#define _DATAREPRESENTATION_
#include <stddef.h>

typedef struct pointset{
	int size;
	int *members;
} pointset;

typedef struct pasta {
	pointset where;
	int howmany;
} pasta;

typedef struct pastabox{
	pasta inleft;
	pasta inright;
} pastabox;

typedef struct cutdata {
	int *matrixrevisited;
	int na;
	int *a;
	int nacomp;
	int *acomp;

	int *tc;
	int nrep;
	int *complementtc;
	int nrepincomp;
	int *pointtorep;
	int *pointtorepincomp;

	pointset *lra;
	int lracard;
	pointset *lnra;
	int lnracard;
	pointset *assoc;

	pointset *lracomp;
	int lracompcard;
	pointset *lnracomp;
	int lnracompcard;
	pointset *assoccomp;

	pointset *m;
	pointset *mcomp;

	int* tab;

	pastabox* box;
} cutdata;

typedef struct dectree {
	struct dectree *left;
	struct dectree *right;
	int label;
	cutdata c;
	int computed;
} dectree; 

typedef struct graph {
	int* matrix;
	int* pos;
	int size;
	int sizeset;
	int* domset;
} graph;

typedef struct setwithinsets {
	pointset *set;
	int size;
} setwithinsets;

#endif


