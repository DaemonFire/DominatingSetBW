#include "../include/treeprimitives.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <errno.h>
#include <sys/stat.h>

int getallleafs(dectree t, int *list){
	int n=0;
	int *lleft;
	int *lright;
	int nleft=0;
	int nright=0;

	if (t.left!=NULL){	
		nleft= getnumberofleafs (*(t.left));
		lleft=(int*)malloc(nleft*sizeof(int));
		getallleafs(*(t.left), lleft);
	}

	if (t.right!=NULL){
		nright=getnumberofleafs(*(t.right));
		lright=(int*)malloc(nright*sizeof(int));
		getallleafs(*(t.right),lright);
	}

	if ((t.left==NULL)&&(t.right==NULL)){
		list[0]=t.label;
		n+=1;
	}
	else {
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

int getnumberofleafs (dectree t){
	int n=0;
	int nleft=0;
	int nright=0;

	if (t.left!=NULL)	
		nleft=getnumberofleafs(*(t.left));
	

	if (t.right!=NULL)
		nright=getnumberofleafs(*(t.right));
	

	if ((t.left==NULL)&&(t.right==NULL))
		n+=1;
	
	else 
		n=nleft+nright;
	

	return n;
}
