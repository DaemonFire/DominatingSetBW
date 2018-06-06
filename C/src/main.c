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

int main (int argc, char** argv){
	if (argc!=2)
		return EXIT_FAILURE;
	graph g = loadgraph(argv[1], 180);

	//graph g = loadgraphformat2(argv[1]);

	printf("g.size=%d\n", g.size);
	//dectree t = loadtree("tiefighter.tree");
	//dectree *t = generateTree(p,g,0);
	dectree *t=generateTreeBW (g);

/*
	FILE *f;
	if ((f=fopen("arbre.tree","r"))==NULL) {
		fprintf(stderr,"Error fopen");
		perror("fopen");
	}


	storetree (*t, f, "0");*/

	printTree (*t);

	pointset x = toplevelalgorithm (*t, g);
	printf("Minimum Dominating Set is of size %d\n",x.size);
	for (int i=0;i<x.size;i++){
		printf("%d %d\n",g.pos[3*x.members[i]+1], g.pos[3*x.members[i]+2]);
	}
	printf("\n");

	//generatePlotFile (*t, g);

	return EXIT_SUCCESS;
}
