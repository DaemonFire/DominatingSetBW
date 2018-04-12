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

dectree *generateTree (pointset p, graph g, int verthor){
	dectree *t;
	t = (dectree*)malloc(sizeof(dectree));
	pointset p1;
	pointset p2;
	p1.size=0;
	p2.size=0;
	p1.members=(int*)malloc(p.size*sizeof(int));
	p2.members=(int*)malloc(p.size*sizeof(int));

	int vnext=0;
	if (p.size==0)
		return NULL;
	else if (p.size==1){

		(*t).label=p.members[0];
		(*t).left=NULL;
		(*t).right=NULL;

	}
	else if (p.size==2){

		(*t).label=-1;
		(*t).left=(dectree*)malloc(sizeof(dectree));
		(*t).right=(dectree*)malloc(sizeof(dectree));

		dectree *t1;
		dectree *t2;
		t1=(dectree*)malloc(sizeof(dectree));
		t2=(dectree*)malloc(sizeof(dectree));
		(*t).left=t1;
		(*t).right=t2;

		(*t1).label=p.members[0];

		(*t1).left=NULL;
		(*t1).right=NULL;

		(*t2).label=p.members[1];
		(*t2).left=NULL;
		(*t2).right=NULL;
	}
	else {
		if (verthor==0){
			int xmoy=0;
			for (int i=0;i<p.size;i++){
				xmoy+=g.pos[2*p.members[i]];
			}
			xmoy=xmoy/p.size;
			
			for (int i=0;i<p.size;i++){
				if (g.pos[2*p.members[i]]<xmoy){
					p1.size++;
					p1.members[p1.size-1]=p.members[i];
				}
				else {
					p2.size++;
					p2.members[p2.size-1]=p.members[i];
				}
			}
			if (p1.size==0){
				p1.size++;
				p1.members[0]=p2.members[p2.size-1];
				p2.size--;
			}
			if (p2.size==0){
				p2.size++;
				p2.members[0]=p1.members[p1.size-1];
				p1.size--;
			}

			vnext=1;
		}

		if (verthor==1){
			int ymoy=0;
			for (int i=0;i<p.size;i++){
				ymoy+=g.pos[2*p.members[i]+1];
			}
			ymoy=ymoy/p.size;
			
			for (int i=0;i<p.size;i++){
				if (g.pos[2*p.members[i]+1]<ymoy){
					p1.size++;
					p1.members[p1.size-1]=p.members[i];
				}
				else {
					p2.size++;
					p2.members[p2.size-1]=p.members[i];
				}
			}
			if (p1.size==0){
				p1.size++;
				p1.members[0]=p2.members[p2.size-1];
				p2.size--;
			}
			if (p2.size==0){
				p2.size++;
				p2.members[0]=p1.members[p1.size-1];
				p1.size--;
			}

			vnext=0;
		}


		(*t).label=-1;

		(*t).left=(dectree*)malloc(sizeof(dectree));
		(*t).left=generateTree(p1,g,vnext);
		(*t).right=(dectree*)malloc(sizeof(dectree));
		(*t).right=generateTree(p2,g,vnext);

	}
	
	return t;
}
