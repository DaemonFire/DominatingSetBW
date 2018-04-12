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

	graph g = loadgraph("antoine2.points", 180);
	pointset p;
	p.size=g.size;
	p.members=(int*)malloc(g.size*sizeof(int));
	printf("g.size=%d\n", g.size);
	for (int i=0;i<g.size;i++){
		p.members[i]=i;
	}
	//dectree t = loadtree("tiefighter.tree");
	dectree *t = generateTree(p,g,0);

	FILE *f;
	if ((f=fopen("arbre.tree","r"))==NULL) {
		fprintf(stderr,"Error fopen");
		perror("fopen");
	}


	storetree (*t, f, "0");
	int x = toplevelalgorithm (*t, g);

	printf("Minimum Dominating Set is of size %d\n",x);
	//generatePlotFile (*t, g);

	return EXIT_SUCCESS;
}
