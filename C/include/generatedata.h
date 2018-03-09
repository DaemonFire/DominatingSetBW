#ifndef _GENERATEDATA_
#define _GENERATEDATA_

#include "../include/datarepresentation.h"


	graph generategraph (int nbpoint, int xrange, int yrange, int threshold);
	
	int storegraph (graph g, char* name);

	graph loadgraph (char* path, int threshold);

	void generateEdges(graph g, int threshold);

	int computeSizeFromFile (char* path);

	//dectree loadtree (char* path);

	//dectree *lookTree (char *path, char *node)

#endif
