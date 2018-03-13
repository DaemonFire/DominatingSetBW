#include "../include/treeprimitives.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <errno.h>
#include <sys/stat.h>


// Compute all leaves of the tree. This function will be useful to compute leaves of subtrees, allowing for cuts in the graph by dividing the decomposition tree
int getallleaves(dectree t, int *list){
	int n=0;
	int *lleft;
	int *lright;
	int nleft=0;
	int nright=0;

	if (t.left!=NULL){						// we get the list of leaves in the first son if it exists, starting by getting the number of leaves in this subtree in order to allocate memory for the list we're going to fill
		nleft= getnumberofleaves (*(t.left));
		lleft=(int*)malloc(nleft*sizeof(int));
		getallleaves(*(t.left), lleft);
	}

	if (t.right!=NULL){						// we do the same for the second son
		nright=getnumberofleaves(*(t.right));
		lright=(int*)malloc(nright*sizeof(int));
		getallleaves(*(t.right),lright);
	}

	if ((t.left==NULL)&&(t.right==NULL)){				// this is the case of the leaf. If the node is a leaf, the list will be a singleton containing only the label of the leaf
		list[0]=t.label;
		n+=1;
	}
	else {								// else, we concatene both list filled before inside the list provided by the caller, minding indexes
		n=nleft+nright;
		for (int i=0; i<nleft;i++){
			list[i]=lleft[i];
		}	

		for (int i=0;i<nright;i++){
			list[nleft+i]=lright[i];
		}
	}
	return n;
}


// This function allows for computation of the number of leaves inside a tree (and therefore, a subtree), which will be useful if we want to allocate memory for lists that will be filled by leaves's labels in getallleaves function
int getnumberofleaves (dectree t){
	int n=0;
	int nleft=0;
	int nright=0;

	if (t.left!=NULL)						// it's a simple enough recursive call on sons and addition of those 2 numbers
		nleft=getnumberofleaves(*(t.left));
	

	if (t.right!=NULL)
		nright=getnumberofleaves(*(t.right));
	

	if ((t.left==NULL)&&(t.right==NULL))				// or, in the case of the leaf, just returning 1
		n+=1;
	
	else 
		n=nleft+nright;
	

	return n;
}
