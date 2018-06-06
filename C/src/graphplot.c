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

	FILE * temp = fopen ("res/config.txt", "w");
	FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
	fprintf(gnuplotPipe, "plot ");
	fprintf(temp, "# x y\n\n");
	int i;	
	int j;
	for (i=0; i<g.size;i++) {
		int friend=0;
		for (j=0; j<i;j++){
			if (g.matrix[g.size*i+j]==1){
				fprintf(temp, "%d %d %d\n", g.pos[3*i+1], g.pos[3*i+2], g.pos[3*i]);
				fprintf(temp, "%d %d %d\n", g.pos[3*j+1], g.pos[3*j+2], g.pos[3*j]);
				fprintf(temp, "\n");
				friend = 1;
			}
		}
		if (friend ==0){
			fprintf(temp, "%d %d %d\n", g.pos[3*i+1], g.pos[3*i+2], g.pos[3*i]);
			fprintf(temp, "%d %d %d\n", g.pos[3*i+1], g.pos[3*i+2], g.pos[3*i]);
			fprintf(temp, "\n");
		}
	}

	fprintf(gnuplotPipe, "\"res/config.txt\" using 1:2 with lines lc rgb \"black\" lw 2 notitle, \\\n \"res/config.txt\" using 1:2:(0.2) with circles fill solid lc rgb \"black\" notitle, \\\n \"res/config.txt\" using 1:2:3 with labels tc rgb \"white\" offset (0,0) font 'Arial Bold' notitle\n");

	fflush(gnuplotPipe);
	return EXIT_SUCCESS;
}


