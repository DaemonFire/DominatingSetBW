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

	graph g = loadgraph("tiefighter.points", 180);
	dectree t = loadtree("tiefighter.tree");

	int x = toplevelalgorithm (t, g);
	printf("Minimum Dominating Set is of size %d\n",x);


	return EXIT_SUCCESS;
}
