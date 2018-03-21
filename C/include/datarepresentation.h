#ifndef _DATAREPRESENTATION_
#define _DATAREPRESENTATION_

#include <stddef.h>

/* 
Data structure for decomposition tree representation. Both pointers will store pointers on sons if they exist or be NULL if they don't exist. The label will be -1 if the node is not a terminal singleton, and any other int if it's a singleton, the value being the point of the graph represented by the leaf 
*/

typedef struct dectree {
	struct dectree *left;
	struct dectree *right;
	int label;
} dectree; 

// This structure allows to store a set of point in a single structure, while having its size available at all time
typedef struct pointset{
	int size;
	int *members;
} pointset;



/* 
This structure is linked to a cut of the decomposition tree t (cut of the left or right son, depending on choiceofson value, 0 or 1). It will be used by the algorithms and contains, the matrix of adjacency limited to the concerned points (giving infos about edges between points of the two sub-graphs thus formed), the numbers and lists of points in each of the subgraphs, the number and list of representants of equivalency classes in both sub-graphs, a list of every points of each subgraph associated to the representant of its equivalency class, the list of representative sets (see second preprocessing) in the primary sub-graph, the list of neighboorhoods of each representative set of the primary sub-graph in the second sub-graph and a list stating the associations between each representative set and his neighboorhood
*/
typedef struct cutdata {
	dectree t;
	int choiceofson;
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

	pointset *m;
	pointset *mcomp;

	int* tab;
} cutdata;



/*
This data structure is a representation of the graph. It is basically its adjacency matrix and a matrix giving the name of each point and its x and y positions. In order to be easier to use, we also add the size, for it to be more accessible
*/

typedef struct graph {
	int* matrix;
	int* pos;
	int size;
} graph;

#endif


