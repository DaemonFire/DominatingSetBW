#include "../include/datarepresentation.h"
#include "../include/graphplot.h"
#include "../include/generatedata.h"
#include "../include/treeprimitives.h"
#include "../include/algorithms.h"

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <fcntl.h>
#include <sys/time.h>
#include <time.h>

int main (int argc, char** argv){
	if (argc!=2)
		return EXIT_FAILURE;
	graph* g = (graph*)malloc(sizeof(graph));
	*g =loadgraph(argv[1], 180);
	graph *h = (graph*)malloc(sizeof(graph));
	*h = *g;
	int* sol = (int*)malloc(g->size*sizeof(int));
	int size=0;


	size = preprocessingsolopoints (g, sol, 180);
	printf ("size=%d\n", size);
	//graph g = loadgraphformat2(argv[1]);

	printf("g.size=%d\n", g->size);
	printf("number of edges=%d\n", getEdgeNumber(*g));
	//dectree t = loadtree("tiefighter.tree");
	//dectree *t = generateTree(p,g,0);

	struct timeval stop, start;
	gettimeofday(&start, NULL);

	dectree *t=generateTreeBW (*g);
	gettimeofday(&stop, NULL);
	int timeToTree=stop.tv_usec - start.tv_usec;
	printf("Time elapsed for tree=%d\n",timeToTree);
/*
	FILE *f;
	if ((f=fopen("arbre.tree","r"))==NULL) {
		fprintf(stderr,"Error fopen");
		perror("fopen");
	}


	storetree (*t, f, "0");*/
/*	int bwmax = getBW(*t,g);
	printf("Boolean-Width=%d\n",bwmax);*/

//	printTree (*t);
	pointset x;
	if (g->size>0)
		x = toplevelalgorithm (t, g);
	else
		x.size=0;

	if (x.size<0)
		x.size=0;
	gettimeofday(&stop, NULL);
	for (int i=0; i<x.size; i++)
		sol[size+i]=x.members[i];
	size+=x.size;

	printf("Minimum Dominating Set is of size %d\n",size);
/*	for (int i=0;i<x.size;i++){
		printf("%d %d\n",g.pos[3*x.members[i]+1], g.pos[3*x.members[i]+2]);
	}*/

	for (int i=0; i<size-x.size; i++){
		printf("(%d, %d),  ", h->pos[2*sol[i]], h->pos[2*sol[i]+1]);
	}
	for (int i=0; i<x.size; i++){
		printf("(%d, %d), ", g->pos[2*sol[size-x.size+i]], g->pos[2*sol[size-x.size+i]+1]);
	}
	printf("\n");
	
	//printf("\n");
	//generatePlotFile (*t, g);


	int timeToSet=stop.tv_usec - start.tv_usec;

	printf("Time elapsed for set=%d\n",timeToSet);


	return EXIT_SUCCESS;
}
