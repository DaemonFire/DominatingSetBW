#ifndef _ALGORITHMS_
#define _ALGORITHMS_
#include "../include/datarepresentation.h"

// This function takes a graph, a decomposition tree and an int stating if we're to cut the left or right son of the tree t to create the
// cut in the decomposition graph
cutdata cutThatTree (graph g, dectree t, int choiceofson);


/*
This first pre-process takes a graph, its decomposition tree and a cut of this decomposition tree as input and computes its equivalency 
classes, the output being the lists of representants of the equivalency classes of both primary and secondary sub-graphs, their numbers
and two lists which associate each point of a subgraph to the representant of its equivalency class 
*/
cutdata firstpreprocess(graph g,  cutdata c);



// This if the function of wonders! It computes twins in each subgraph by checking if their induced vectors in the adjacency matrix match
int matchTwins (cutdata c, int* twinsleft, int* twinsright);


/*
This function is another way of computing twins, using the partition refinment procedure but it has not yet been tailored to fit this
project. At least it is statically ok and doesn't cause compilation problem. TODO: Tailor this shit to be used instead of the previous function because it is more badass.
*/
int findTwins(int n,int* mat,int* twins);


/* 
This is the second pre-processing function. We give it a tree and a cutdata on which the first preprocess has been carried on 
and it uses informations on representants to compute the list of representant sets, their respectives neighboorhoods and associations
between those two sets of data, to be stored in the cutdata
*/
cutdata secondpreprocess (cutdata c, graph g);

cutdata thirdpreprocess (cutdata c, graph g);

int toplevelalgorithm (dectree t, graph g);

cutdata *stepalgorithm (dectree t, dectree tparticular, graph g);


#endif
