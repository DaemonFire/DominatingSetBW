#include "../include/datarepresentation.h"
#include "../include/graphplot.h"
#include "../include/generatedata.h"

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
	graph g=loadgraph("straight.points", 180);
	dectree t=loadtree("straight.tree");

	generatePlotFile (t, g);

	return EXIT_SUCCESS;
}
