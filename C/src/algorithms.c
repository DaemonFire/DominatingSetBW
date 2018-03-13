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

	if((t.left==NULL)||(t.right==NULL)){
		return EXIT_FAILURE;
	}

	int nleft=0;
	int nright=0;
	int *verticesleft=NULL;
	int *verticesright=NULL;

	if (choiceofson==0){
		nleft=getnumberofleaves (*(t.left));
		verticesleft=(int*)malloc(nleft*sizeof(int));
		getallleaves(*(t.left), verticesleft);
		nright=g.size-nleft;
		verticesright=(int*)malloc(nright*sizeof(int));

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
				if (k==nright)
					break;
			}
		}


		
	}

	else if (choiceofson==1){
		nright=getnumberofleaves (*(t.right));
		verticesright=(int*)malloc(nright*sizeof(int));
		getallleaves(*(t.right), verticesright);
		nleft=g.size-nright;
		verticesleft=(int*)malloc(nleft*sizeof(int));
		
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
				if (k==nleft)
					break;
			}
		}


	}


	int *matrixwrtcut = (int*)malloc(nleft*nright*sizeof(int));

	for (int i=0;i<nleft;i++){
		for (int j=0;j<nright;j++){
			matrixwrtcut[i*nright+j]=g.matrix[verticesleft[i]*g.size+verticesright[j]];
		}
	printf("what happen %d\n",i);
	}
	for (int i=0;i<nleft;i++){
		for (int j=0;j<nright;j++){
			printf("%d|",matrixwrtcut[i*nright+j]);
		}
		printf("\n");
	}

	int *twinsleft = (int*)malloc(2*nleft*nleft*sizeof(int));
	int *twinsright = (int*)malloc(2*nright*nright*sizeof(int));
	//printf("what happen\n");
	matchTwins(nleft,nright,matrixwrtcut,twinsleft,twinsright);
	//printf("what happen\n");
	for (int i=0;i<nleft*nleft;i++){
		printf("Voici une classe à gauche: ");
		for (int j=0;j<2;j++){
			printf("%d ",twinsleft[i*2+j]);
		}
		printf("\n");
	}

	for (int i=0;i<nright*nright;i++){
		printf("Voici une classe à droite: ");
		for (int j=0;j<2;j++){
			printf("%d ",twinsright[i*2+j]);
		}
		printf("\n");
	}

	t.pointtorep=(int*)malloc(2*nleft*sizeof(int));
	t.pointtorepincomp=(int*)malloc(2*nright*sizeof(int));
	for (int i=0;i<nleft;i++){
		printf("Sommet à gauche : %d\n",verticesleft[i]);
	}
	for (int i=0;i<nright;i++){
		printf("Sommet à droite : %d\n",verticesright[i]);
	}
	for (int i=0; i<2*nleft;i++){
		t.pointtorep[i]=-1;
	}

	for (int i=0;i<2*nright;i++){
		t.pointtorepincomp[i]=-1;
	}

	for (int i=0;i<nleft;i++){
		t.pointtorep[i*2]=verticesleft[i];
		for (int j=0;j<nleft*nleft;j++){
			int a = twinsleft[j*2];
			int b = twinsleft[j*2+1];
		
			if (a==i){
				if (verticesleft[b]>t.pointtorep[i*2+1]){
					t.pointtorep[i*2+1]=verticesleft[b];
				}
			}
			if (b==i){
				if (verticesleft[a]>t.pointtorep[i*2+1]){
					t.pointtorep[i*2+1]=verticesleft[a];
				}
			}
		}
		if ((t.pointtorep[i*2+1]==-1)||(t.pointtorep[i*2+1]<verticesleft[i])){
			t.pointtorep[i*2+1]=verticesleft[i];
		}
	}
	for (int i=0;i<nleft;i++){
		printf("Voici une classe à gauche: ");
		for (int j=0;j<2;j++){
			printf("%d ",t.pointtorep[i*2+j]);
		}
		printf("\n");
	}

	for (int i=0;i<nright;i++){
		t.pointtorepincomp[i*2]=verticesright[i];
		for (int j=0;j<nright*nright;j++){
			int a = twinsright[j*2];
			int b = twinsright[j*2+1];
		
			if (a==i){
				if (verticesright[b]>t.pointtorepincomp[i*2+1]){
					t.pointtorepincomp[i*2+1]=verticesright[b];
				}
			}
			if (b==i){
				if (verticesright[a]>t.pointtorepincomp[i*2+1]){
					t.pointtorepincomp[i*2+1]=verticesright[a];
				}
			}
		}
		if ((t.pointtorepincomp[i*2+1]==-1)||(t.pointtorepincomp[i*2+1]<verticesright[i])){
			t.pointtorepincomp[i*2+1]=verticesright[i];
		}
	}
	for (int i=0;i<nright;i++){
		printf("Voici une classe à droite: ");
		for (int j=0;j<2;j++){
			printf("%d ",t.pointtorepincomp[i*2+j]);
		}
		printf("\n");
	}
	
	t.tc=(int*)malloc(nleft*sizeof(int));
	t.complementtc=(int*)malloc(nright*sizeof(int));
	int cursor = 0;
	for (int i=0;i<nleft;i++){
		int here=0;
		for (int j=0;j<nleft;j++){
			if (t.pointtorep[2*i+1]==t.tc[j]){
				here=1;
				break;
			}
		}
		if (here==0){
			t.tc[cursor]=t.pointtorep[2*i+1];
			cursor++;
		}
	}

	t.nrep=cursor;
	realloc(t.tc,cursor*sizeof(int));
	cursor=0;
	for (int i=0;i<nright;i++){
		int here=0;
		for (int j=0;j<nright;j++){
			if (t.pointtorepincomp[2*i+1]==t.complementtc[j]){
				here=1;
				break;
			}
		}
		if (here==0){
			t.complementtc[cursor]=t.pointtorepincomp[2*i+1];
			cursor++;
		}
	}
	t.nrepincomp=cursor;
	realloc(t.complementtc,cursor*sizeof(int));

	t.matrixrevisited=matrixwrtcut;

	return EXIT_SUCCESS;

}

int matchTwins (int nleft, int nright, int* mat, int* twinsleft, int* twinsright){
	int i,j,k;
	int cursor=0;


	for(i=0;i<2*nleft*nleft;i++)
		twinsleft[i]=-1;

	for(i=0;i<2*nright*nright;i++)
		twinsright[i]=-1;

	for (i=0;i<nleft;i++){
		for (j=i+1;j<nleft;j++){	
			int tw=1;	

			for (k=0;k<nright;k++){
			
				if(mat[i*nright+k]!=mat[j*nright+k]){
				
					tw=0;
					break;
				}
			}	
			
			if (tw==1){	
				printf("%d, %d\n",i,j);		
				twinsleft[cursor]=i;
				cursor++;
				twinsleft[cursor]=j;
				cursor++;
			}
			
		} 		
	}
	cursor = 0;
	for (i=0;i<nright;i++){
		for (j=i+1;j<nright;j++){	
			int tw=1;	

			for (k=0;k<nleft;k++){
			
				if(mat[k*nright+i]!=mat[k*nright+j]){
				
					tw=0;
					break;
				}
			}	
			
			if (tw==1){		
				printf("%d, %d\n",i,j);			
				twinsright[cursor]=i;
				cursor++;
				twinsright[cursor]=j;
				cursor++;
			}
			
		} 		
	}

	return EXIT_SUCCESS;
}


int findTwins(int n,int* mat,int* twins){
	int i,j,k;	
	int cursor=0;
	int cursortwins=0;
	
	for (i=0;i<n*n;i++)
		twins[i]=0;
	
	int* partition;
	printf("sup\n");
	partition=(int*)malloc(n*n*sizeof(int));
	printf("sup\n");

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

int secondepreprocess (dectree t){

	t.lra=NULL;
	t.lrcompa=NULL;
	int **nextLevel=NULL;
	int **lastLevel=(int*)malloc(sizeof(int*));
	lastLevel[0]=(int*)malloc(sizeof(int));
	lastLevel[0][0]=-1;
	int sizeoflast=1;
	
	while (lastLevel!=NULL) {
		for (int i=0; i<sizeoflast;i++){
			int j=0;
			int *r;
			while (lastLevel[i][j]!=-1){
				realloc(r,(j+1)*sizeof(int));
				r[j]=lastLevel[i][j];
			}

			int k=0;

			while (t.tc[k]!=NULL){

				int *rprime = (int*)malloc((j+2)*sizeof(int));

				for (int l=0;l<j+1;l++){
					rprime[l]=r[l];
				}

				rprime[j+1]=t.tc[k];
				int *n=NULL;
				int numn=0;

				for (int l=0;l<j+2;l++){

					for (int m=0; m<t.nrepincomp; m++){

						if(t.matrixrevisited[rprime[l]*t.nrepincomp+m]==1){

							int alreadyin=0;
							for (int o=0; o<numn;o++){

								if (n[o]==t.complementtc[m]){
									alreadyin=1;
									break;
								}

							}

							if (alreadyin==0){
								numn++;
								realloc(n,numn*sizeof(int));
								n[numn-1]=t.complementtc[m];
							}
						}
					}
				}

				
			}
		}
	}

	return EXIT_SUCCESS;
}
