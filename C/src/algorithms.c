#include "../include/algorithms.h"
#include "../include/treeprimitives.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <errno.h>
#include <sys/stat.h>
#define MAT(mat,x,y,n) (mat[x*n+y])

int firstpreprocess(graph g, dectree t, int choiceofson){

	if((t.left==NULL)||(t.right==NULL){
		return EXIT_FAILURE;
	}

	int nleft=0;
	int nright=0;
	int *verticesleft=(int*)malloc(g.size*sizeof(int));
	int *verticesright=(int*)malloc(g.size*sizeof(int));
	
	if (choiceofson==0){
		nleft=int getallleafs(*(t.left), *verticesleft);
		nright=g.size-nleft;
		
		int i=0;
		int j=0;
		int k=0;
		for (i=0;i<g.size;i++){
			int inleft=0;
			for (j=0;j<nleft;j++){
				if (verticesleft[j]==i){
					inleft=1;
					break;
				}
			}
			if (inleft==0){
				verticesright[k]=i;
				k++;
				if (k==nright-1)
					break;
			}
		}



		
	}

	else if (choiceofson==1){
		nright=int getallleafs(*(t.right), *verticesright);
		nleft=g.size-nright;
		
		int i=0;
		int j=0;
		int k=0;
		for (i=0;i<g.size;i++){
			int inright=0;
			for (j=0;j<nright;j++){
				if (verticesright[j]==i){
					inright=1;
					break;
				}
			}
			if (inright==0){
				verticesleft[k]=i;
				k++;
				if (k==nkeft-1)
					break;
			}
		}


	}


	int *matrixwrtcut = (int*)malloc(nleft*nright*sizeof(int));
		
	for (i=0;i<nleft;i++){
		for (j=0;j<nright;j++){
			matrixwrtcut[i*nleft+j]=g.matrix[verticesleft[i]*g.size+verticesright[j]];
		}
	}
	
	int *twinsleft = (int*)malloc(nleft*nleft*sizeof(int));
	findTwins(g.size,matrixwrtcut,twinsleft);
	int *twinsright = (int*)malloc(nright*nright*sizeof(int));
	findTwins(g.size,matrixwrtcut,twinsright);

}


int findTwins(int n,int* mat,int* twins){
	int i,j,k;	
	int cursor=0;
	int cursortwins=0;
	
	for (i=0;i<n*n;i++)
		twins[i]=0;
	
	int* partition;
	partition=(int*)malloc(n*n*sizeof(int));

	for (i=0;i<n*n;i++){
		if (i<n)
			partition[i]=1;
		else
			partition[i]=0;
	}

	for(i=0;i<n;i++){
		int* tmp;
		tmp=(int*)malloc(n*n*sizeof(int));

		for (j=0;j<n*n;j++)
			tmp[j]=0;
		cursor=0;
		int tw=0;
		
		for (j=0;j<n;j++){
			for(k=0;k<n;k++){
				if ((MAT(mat,i,k,n)==1)&&(MAT(partition,j,k,n)==1)){
					MAT(tmp,cursor,k,n)=1;
					tw++;
				}
			}	
			if (tw!=0){
				cursor++;
				tw=0;
			}

			for(k=0;k<n;k++){
				if ((MAT(mat,i,k,n)==0)&&(MAT(partition,j,k,n)==1)){
					MAT(tmp,cursor,k,n)=1;
					tw++;
				}
			}

			if (tw!=0){
				cursor++;
				tw=0;
			}
		}
		
		for (j=0;j<n*n;j++)
			partition[j]=tmp[j];
		free(tmp);
	}


	
	for (i=0;i<n;i++){
		int tw3=0;
		for (j=0;j<n;j++){
		
			if (MAT(partition,i,j,n)==1){
				if(tw3<2){
					MAT(twins,cursortwins,tw3,2)=j+1;
					tw3++;
				}
				else {
					printf("Module de cardinal n>2\n");
					tw3++;
				}
			}
		}
		
		if(tw3==2)
			cursortwins++;
		else{
			MAT(twins,cursortwins,0,2)=0;
			MAT(twins,cursortwins,1,2)=0;
		}
	}
	


	for (i=0;i<n;i++)
		MAT(mat,i,i,n)=1;
	

	for (i=0;i<n*n;i++){
		if (i<n)
			partition[i]=1;
		else
			partition[i]=0;
	}

	for(i=0;i<n;i++){
		int* tmp;
		tmp=(int*)malloc(n*n*sizeof(int));

		for (j=0;j<n*n;j++)
			tmp[j]=0;

		cursor=0;
		int tw=0;
		
		for (j=0;j<n;j++){
			for(k=0;k<n;k++){
				if ((MAT(mat,i,k,n)==1)&&(MAT(partition,j,k,n)==1)){
					MAT(tmp,cursor,k,n)=1;
					tw++;
				}
			}
	
			if (tw!=0){
				cursor++;
				tw=0;
			}

			for(k=0;k<n;k++){
				if ((MAT(mat,i,k,n)==0)&&(MAT(partition,j,k,n)==1)){
					MAT(tmp,cursor,k,n)=1;
					tw++;
				}
			}

			if (tw!=0){
				cursor++;
				tw=0;
			}
		}
		
		for (j=0;j<n*n;j++)
			partition[j]=tmp[j];
		free(tmp);
	}


	
	for (i=0;i<n;i++){
		int tw3=0;
		for (j=0;j<n;j++){
		
			if (MAT(partition,i,j,n)==1){
				if(tw3<2){
					MAT(twins,cursortwins,tw3,2)=j+1;
					tw3++;
				}
				else {
					printf("Module de cardinal n>2\n");
					tw3++;
				}
			}
		}
		
		if(tw3==2){
			cursortwins++;
		}
		else{
			MAT(twins,cursortwins,0,2)=0;
			MAT(twins,cursortwins,1,2)=0;
		}
	}

	for (i=0;i<n;i++)
		MAT(mat,i,i,n)=0;
	

	return EXIT_SUCCESS;
}
