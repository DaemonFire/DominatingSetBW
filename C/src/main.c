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
	//graph g = generategraph (10,1000,1000,450);
	//storegraph(g, "straight.points");
	graph g=loadgraph("stardestroyer.points", 180);
	dectree t=loadtree("stardestroyer.tree");

	int *l;

	int n = getnumberofleaves(t);
	l=(int*)malloc(n*sizeof(int));
	
	getallleaves(*(t.right),l);
	
	firstpreprocess(g, t, 0);

	//generatePlotFile (t, g);

	return EXIT_SUCCESS;
}
