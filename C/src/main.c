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
	printf("DO YOU BELIEVE IT? DOMINATING SET IS OF SIZE %d\n",x);
	/*cutdata c = cutThatTree (g, *(t.right), 0);

	for (int i=0;i<g.size;i++){
		for (int j=0;j<g.size;j++){
			printf("%d ", g.matrix[i*g.size+j]);
		}
		printf("\n");
	}

	c=firstpreprocess(g, t, c);
	printf("choice of son: %d na = %d ,nacomp = %d nrep = %d nrepincomp= %d\n", c.choiceofson, c.na, c.nacomp, c.nrep, c.nrepincomp);


	for (int i=0;i<c.nrep;i++){
		for (int j=0;j<c.nrepincomp;j++){
			printf("%d ", c.matrixrevisited[i*c.nrepincomp+j]);
		}
		printf("\n");
	}

	c=secondpreprocess (t, c, g);
	printf("choice of son: %d na = %d ,nacomp = %d nrep = %d nrepincomp= %d, c.lra= %d, c.lnra= %d\n", c.choiceofson, c.na, c.nacomp, c.nrep, c.nrepincomp, c.lracard, c.lnracard);

	c=thirdpreprocess (t,c,g);

	printf("lnra = ");
	for (int i=0;i<c.lnracard;i++){
		for (int j=0;j<c.lnra[i].size;j++)
			printf("%d ", c.lnra[i].members[j]);
		printf("\n");
	}*/
	//generatePlotFile (t, g);

	return EXIT_SUCCESS;
}
