#include "../include/datarepresentation.h"
#include "../include/graphplot.h"

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <fcntl.h>
#include <string.h>

#define GNUPLOT_PATH "/usr/bin/gnuplot"
#define PLOTFILE_PATH "res/config.txt"

int generatePlotFile (dectree t, graph g){

//open file in which we will write down our data
	FILE * temp = fopen ("res/config.txt", "w");

//open the GnuPlot connector
	FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");

	fprintf(gnuplotPipe, "plot ");
	fprintf(temp, "# x y\n\n");

	int i;	
	int j;
	for (i=0; i<g.size;i++) {
		int friend=0;
		for (j=0; j<i;j++){
			if (g.matrix[g.size*i+j]==1){
				fprintf(temp, "%d %d %d\n", g.pos[3*i+1], g.pos[3*i+2], g.pos[3*i]);   // we use adjacency matrix to get edges to plot and 
				fprintf(temp, "%d %d %d\n", g.pos[3*j+1], g.pos[3*j+2], g.pos[3*j]);   // write thoose down in our config file
				fprintf(temp, "\n");
				friend = 1;
			}
		}
		if (friend ==0){      														// if a node has no neighbor (i.e he's isolated), he
			fprintf(temp, "%d %d %d\n", g.pos[3*i+1], g.pos[3*i+2], g.pos[3*i]);	// wouldn't be written down in our config file
			fprintf(temp, "%d %d %d\n", g.pos[3*i+1], g.pos[3*i+2], g.pos[3*i]);	// we know this thanks to our "friend" value and if
			fprintf(temp, "\n");													// this value is null, we write the node in our config
		}																			// file as a neighbor of himself (that won't alter the
	}																				// adjacency matrix as it's only gnuplot commands


// this long line is just a sequency of gnuplot instruction to plot, first the edges as lines, then the nodes as cycles, then the labels of the nodes
	fprintf(gnuplotPipe, "\"res/config.txt\" using 1:2 with lines lc rgb \"black\" lw 2 notitle, \\\n \"res/config.txt\" using 1:2:(0.2) with circles fill solid lc rgb \"black\" notitle, \\\n \"res/config.txt\" using 1:2:3 with labels tc rgb \"white\" offset (0,0) font 'Arial Bold' notitle\n");

// we flush the gnuplotPipe buffer in order for the display to appear upon termination of this function. 
	fflush(gnuplotPipe);
	return EXIT_SUCCESS;
}


