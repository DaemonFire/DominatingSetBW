#ifndef _GENERATEDATA_
#define _GENERATEDATA_

#include "../include/datarepresentation.h"
#include <stdio.h>


// Generating a geometric graph with nbpoint points ranging from 0 to xrange in x and from 0 to yrange in y. Points are connected by an edge if their distance is lower than the threshold
	graph generategraph (int nbpoint, int xrange, int yrange, int threshold);
	

// Once we've generated a graph g, this function allows you to save this graph in a file, the path being name, for future usage
	int storegraph (graph g, char* name);


// We're not going to construct the graph in static code, or we'll never be able to run tests on large data sets. This function allows to load a file at path and build the graphs containing all points given in the file and using the threshold to compute edges
	graph loadgraph (char* path, int threshold);


// Whether you're generating a graph or loading one, this function will allow you to compute edges for points of g, with respect to the given threshold
	void generateEdges(graph g, int threshold);


// Basically the same as the previous one without the initializations and specific to one node
	dectree loadtree (char* path);


// This function computes the number of points of the graph by counting lines in the data file until it reaches the end of the document
	dectree lookTree (FILE *f, char *node);


// Once your graph is loaded, it's wise to load up the decomposition tree for this graph, which will be used by the algorithm we implement here.
	int computeSizeFromFile (char* path);

	int storetree (dectree t, FILE* f, char* node);

#endif
