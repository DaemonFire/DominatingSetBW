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
	printf(" ===== %s =====\n", argv[1]);
	if (argc!=4)
		return EXIT_FAILURE;
	graph* g = (graph*)malloc(sizeof(graph));
	int threshold = atoi(argv[2]);
	*g =loadgraph(argv[1], threshold);
	graph *h = (graph*)malloc(sizeof(graph));
	*h = *g;
	int* sol = (int*)malloc(g->size*sizeof(int));
	int size=0;
	size = preprocessingsolopoints (g, sol, threshold);
	int sizeinit=size;
	//graph g = loadgraphformat2(argv[1]);

	printf("g.size=%d\n", g->size);
	printf("number of edges=%d\n", getEdgeNumber(*g));
	//dectree t = loadtree("tiefighter.tree");
	//dectree *t = generateTree(p,g,0);
	graph** components = (graph**)malloc(g->size*sizeof(graph*));
	int ncomp = computeconnexcomposants (g, components, threshold);
	struct timeval stop, start;
	long timeToSet = 0;
	for (int i=0; i<ncomp; i++){
		//printf("composant of size %d\n",components[i]->size);
		dectree *t;
		if (argv[3][0]=='n')
			t=generateTreeBW (*components[i]);
		else
			t=generateTreeBWcrude (*components[i], atoi(argv[3]));
		pointset x;
		if (components[i]->size>0){
			gettimeofday(&start, NULL);
			x = toplevelalgorithm (t, components[i]);
			gettimeofday(&stop, NULL);
			timeToSet+=1000000*(stop.tv_sec-start.tv_sec)+stop.tv_usec - start.tv_usec;
			printf("Time for this composant is %ld\n", 1000000*(stop.tv_sec-start.tv_sec)+stop.tv_usec - start.tv_usec);
		}	
		else
			x.size=0;
	
		if (x.size<0)
			x.size=0;
		
		for (int j=0; j<x.size; j++)
			sol[size+j]=x.members[j];
		size+=x.size;
		components[i]->sizeset=x.size;
		components[i]->domset=x.members;
	}
	for (int i=0; i<sizeinit; i++){
		printf("(%d, %d),  ", h->pos[2*sol[i]], h->pos[2*sol[i]+1]);
	}
	printf("|");
	int cursor =sizeinit;
	for (int i=0; i<ncomp; i++){
		for (int j=0; j<components[i]->sizeset; j++){
			printf("(%d, %d), ", components[i]->pos[2*components[i]->domset[j]], components[i]->pos[2*components[i]->domset[j]+1]);
		}
		cursor+=components[i]->sizeset;
		printf("|");
	}
	printf("\n");
	printf("There are %d composants\n", ncomp);
	printf("Minimum Dominating Set is of size %d\n",size);
	//printf("\n");
	//generatePlotFile (*t, g);
	printf("Time elapsed for set=%ld\n",timeToSet);
	printf("\n\n\n");

	return EXIT_SUCCESS;
}
