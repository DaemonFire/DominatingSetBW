#ifndef _TREEPRIMITIVES_
#define _TREEPRIMITIVES_
#include "../include/datarepresentation.h"

// Compute all leaves of the tree. This function will be useful to compute leaves of subtrees, allowing for cuts in the graph by dividing the decomposition tree
int getallleaves(dectree t, int *list);


// This function allows for computation of the number of leaves inside a tree (and therefore, a subtree), which will be useful if we want to allocate memory for lists that will be filled by leaves's labels in getallleaves function
int getnumberofleaves (dectree t);

dectree *generateTree (pointset p, graph g, int verthor);
#endif
