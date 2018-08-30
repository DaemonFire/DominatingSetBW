#ifndef _ALGORITHMS_
#define _ALGORITHMS_
#include "../include/datarepresentation.h"

int getnumberofnodes(dectree* t);

int fillThevoid (dectree* t, graph* g);

cutdata cutThatTree (graph* g, dectree* t);

int firstpreprocess(graph* g,  cutdata* c);

int matchTwins (cutdata c, int* twinsleft, int* twinsright);

int findTwins(int n,int* mat,int* twins);

int secondpreprocess (cutdata* c, graph* g);

int thirdpreprocess (cutdata* c, graph* g);

pointset toplevelalgorithm (dectree* t, graph* g);

int stepalgorithm (dectree* t, graph* g);

pointset computeDS (dectree* t, int much, int a, int ac);

int getBW (dectree* t, graph* g);

#endif
