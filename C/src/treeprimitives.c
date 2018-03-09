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
	if (t.left!=NULL){	
		n+=getallleafs(*((dectree*)(t.left)), list);
	}
	int *l=list+n;
	if (t.right!=NULL){
		n+=getallleafs(*((dectree*)(t.right)),l);
	}
	if ((t.left==NULL)&&(t.right==NULL)){
		n+=1;
		list=&(t.label);
	}
	return n;
}
