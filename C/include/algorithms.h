#ifndef _ALGORITHMS_
#define _ALGORITHMS_
#include "../include/datarepresentation.h"

int firstpreprocess(graph g, dectree t, int choiceofson);

int findTwins(int n,int* mat,int* twins);

int matchTwins (int nleft, int nright, int* mat, int* twinsleft, int* twinsright);

#endif
