#ifndef _ALGORITHMS_
#define _ALGORITHMS_
#include "../include/datarepresentation.h"


cutdata cutThatTree (graph g, dectree t, int choiceofson);

cutdata firstpreprocess(graph g,  cutdata c);

int matchTwins (cutdata c, int* twinsleft, int* twinsright);

int findTwins(int n,int* mat,int* twins);

cutdata secondpreprocess (cutdata c, graph g);

cutdata thirdpreprocess (cutdata c, graph g);

pointset toplevelalgorithm (dectree t, graph g);

dectree stepalgorithm (dectree t, graph g);

pointset computeDS (dectree t, int muchleft, int aleft, int acleft, int muchright, int bright, int bcright);

#endif
