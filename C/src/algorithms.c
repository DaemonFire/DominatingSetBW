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

int tocompute;
dectree** nodestocompute;

int getnumberofnodes(dectree* t){
	int n=1;
	if (t->left!=NULL)
		n+=getnumberofnodes(t->left);
	if (t->right!=NULL)
		n+=getnumberofnodes(t->right);
	return n;
}

int fillThevoid (dectree* t, graph* g){
	t->computed=0;
	nodestocompute[tocompute]=t;
	tocompute++;
	if (t->left!=NULL)
		fillThevoid(t->left, g);
	if (t->right!=NULL)
		fillThevoid(t->right, g);
	return EXIT_SUCCESS;
}

cutdata cutThatTree (graph* g, dectree* t){
	cutdata c;
	c.na=0;
	c.nacomp=0;
	c.a=NULL;
	c.acomp=NULL;

	c.na=getnumberofleaves (*(t));
	c.a=(int*)malloc(c.na*sizeof(int));
	getallleaves(*(t), c.a);
	
	c.nacomp=g->size-c.na;

	c.acomp=(int*)malloc(c.nacomp*sizeof(int));

	int i=0;
	int j=0;
	int k=0;
	for (i=0;i<g->size;i++){
		int ina=0;
		for (j=0;j<c.na;j++){
			if (c.a[j]==i){
				ina=1;
				break;
			}
		}
		if (ina==0){
			c.acomp[k]=i;
			k++;
			if (k==c.nacomp)
				break;
		}
	}

	c.matrixrevisited = (int*)malloc(c.na*c.nacomp*sizeof(int));


	for (int i=0;i<c.na;i++){
		for (int j=0;j<c.nacomp;j++)
			c.matrixrevisited[i*c.nacomp+j]=g->matrix[c.a[i]*g->size+c.acomp[j]];
	}
	return c;
}


int firstpreprocess(graph* g,  cutdata* c){

	int *twinsleft = (int*)malloc(2*c->na*c->na*sizeof(int));
	int *twinsright = (int*)malloc(2*c->nacomp*c->nacomp*sizeof(int));

	matchTwins(*c ,twinsleft,twinsright);

	c->pointtorep=(int*)malloc(2*c->na*sizeof(int));
	c->pointtorepincomp=(int*)malloc(2*c->nacomp*sizeof(int));
	
	for (int i=0; i<2*c->na;i++){
		c->pointtorep[i]=-1;
	}

	for (int i=0;i<2*c->nacomp;i++){
		c->pointtorepincomp[i]=-1;
	}

	for (int i=0;i<c->na;i++){
		c->pointtorep[i*2]=c->a[i];

		for (int j=0;j<c->na*c->na;j++){
			int a = twinsleft[j*2];
			int b = twinsleft[j*2+1];

			if (a==i){
				if (c->a[b]>c->pointtorep[i*2+1])
					c->pointtorep[i*2+1]=c->a[b];
			}
			if (b==i){
				if (c->a[a]>c->pointtorep[i*2+1])
					c->pointtorep[i*2+1]=c->a[a];
			}
		}

		if (c->pointtorep[i*2+1]<c->a[i])
			c->pointtorep[i*2+1]=c->a[i];
	}


	for (int i=0;i<c->nacomp;i++){
		c->pointtorepincomp[i*2]=c->acomp[i];
		for (int j=0;j<c->nacomp*c->nacomp;j++){
			int a = twinsright[j*2];
			int b = twinsright[j*2+1];
		
			if (a==i){
				if (c->acomp[b]>c->pointtorepincomp[i*2+1])
					c->pointtorepincomp[i*2+1]=c->acomp[b];
			}
			if (b==i){
				if (c->acomp[a]>c->pointtorepincomp[i*2+1])
					c->pointtorepincomp[i*2+1]=c->acomp[a];
			}
		}

		if (c->pointtorepincomp[i*2+1]<c->acomp[i])
			c->pointtorepincomp[i*2+1]=c->acomp[i];
	}

	c->tc=(int*)malloc(c->na*sizeof(int));
	for (int i=0;i<c->na;i++)
		c->tc[i]=-1;
	c->complementtc=(int*)malloc(c->nacomp*sizeof(int));
	int cursor = 0;
	for (int i=0;i<c->na;i++){
		int here=0;
		for (int j=0;j<c->na;j++){
			if (c->pointtorep[2*i+1]==c->tc[j]){
				here=1;
				break;
			}
		}
		if (here==0){
			c->tc[cursor]=c->pointtorep[2*i+1];			
			cursor++;
		}
	}

	c->nrep=cursor;

	for (int i=0;i<c->nacomp;i++)
		c->complementtc[i]=-1;

	cursor=0;
	for (int i=0;i<c->nacomp;i++){
		int here=0;
		for (int j=0;j<c->nacomp;j++){
			if (c->pointtorepincomp[2*i+1]==c->complementtc[j]){
				here=1;
				break;
			}
		}
		if (here==0){
			c->complementtc[cursor]=c->pointtorepincomp[2*i+1];
			cursor++;
		}
	}
	c->nrepincomp=cursor;

	return EXIT_SUCCESS;

}


int matchTwins (cutdata c, int* twinsleft, int* twinsright){
	int i,j,k;
	int cursor=0;


	for(i=0;i<2*c.na*c.na;i++)
		twinsleft[i]=-1;

	for(i=0;i<2*c.nacomp*c.nacomp;i++)
		twinsright[i]=-1;

	for (i=0;i<c.na;i++){
		for (j=i+1;j<c.na;j++){
			int tw=1;	

			for (k=0;k<c.nacomp;k++){
				if(c.matrixrevisited[i*c.nacomp+k]!=c.matrixrevisited[j*c.nacomp+k]){
					tw=0;
					break;
				}
			}	
			
			if (tw==1){
				twinsleft[cursor]=i;
				cursor++;
				twinsleft[cursor]=j;
				cursor++;
			}
			
		} 		
	}

	cursor = 0;

	for (i=0;i<c.nacomp;i++){
		for (j=i+1;j<c.nacomp;j++){	
			int tw=1;	

			for (k=0;k<c.na;k++){
				if(c.matrixrevisited[k*c.nacomp+i]!=c.matrixrevisited[k*c.nacomp+j]){
					tw=0;
					break;
				}
			}	
			
			if (tw==1){		
				twinsright[cursor]=i;
				cursor++;
				twinsright[cursor]=j;
				cursor++;
			}
		} 		
	}

	return EXIT_SUCCESS;
}



/*
This function is another way of computing twins, using the partition refinment procedure but it has not yet been tailored to fit this
project. At least it is statically ok and doesn't cause compilation problem. TODO: Tailor this shit to be used instead of the previous function because it is more badass.
*/
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


int secondpreprocess (cutdata* c, graph* g){
/*
	c->lra=NULL;
	c->lnra=NULL;
	c->lracomp=NULL;
	c->lnracomp=NULL;
	c->assoc=NULL;
	c->assoccomp=NULL;
*/
	c->lra=(pointset*)malloc(4*c->nrep*c->nrep*c->nrep*c->nrep*sizeof(pointset));	//TODO: Find more accurate measures of memory needed for 
	c->lnra=(pointset*)malloc(4*c->nrep*c->nrep*c->nrep*c->nrep*sizeof(pointset)); //those allocations
	c->lracard=1;
	c->lnracard=1;
	pointset *nextLevel=(pointset*)malloc(4*c->nrep*c->nrep*c->nrep*c->nrep*sizeof(pointset));
	pointset s;
	s.size=0;						
	pointset *lastLevel=(pointset*)malloc(4*c->nrep*c->nrep*c->nrepincomp*c->nrepincomp*sizeof(pointset));
	lastLevel[0]=s;
	int sizeoflast=1;
	int sizeofnext=0;
	c->assoc=(pointset*)malloc(2*c->nrep*c->nrep*c->nrepincomp*c->nrepincomp*sizeof(pointset));

	c->lra[0]=s;
	c->lnra[0]=s;
	c->assoc[0]=s;
	c->assoc[1]=s;

	while (sizeoflast!=0){
		for (int i=0; i<sizeoflast;i++){
			pointset r;
			r.size = lastLevel[i].size;

			if(r.size!=0)
				r.members=(int*)malloc((r.size)*sizeof(int));

			for (int j=0; j< r.size;j++)
				r.members[j]=lastLevel[i].members[j];
			for (int j=0; j<c->nrep; j++){
				pointset rprime;
				rprime.size=0;
				int alreadyin = 0;

				rprime.members=(int*)malloc((r.size+1)*sizeof(int));

				for (int k=0;k<r.size;k++){
					rprime.size++;

					rprime.members[k]=r.members[k];
					if (rprime.members[k]==c->tc[j]){
						alreadyin = 1;
						break;
					}
				}

				if (alreadyin==0) {
					rprime.size++;
					rprime.members[rprime.size-1]=c->tc[j];		
					pointset n;
					n.members=(int*)malloc(c->nrepincomp*sizeof(int));								
					n.size=0;
					
					for (int l=0; l<c->nrepincomp; l++){
						for (int k=0;k<rprime.size;k++){
							if(g->matrix[rprime.members[k]*g->size + c->complementtc[l]]==1){
									n.size++;
									n.members[n.size-1]=c->complementtc[l];
									break;
							}
						}
					}

					int alreadyin = 0;
					
					if (n.size==0)
						alreadyin = 1;
					for (int k=0; k<c->lnracard; k++) {
						int common = 0;
						if (c->lnra[k].size==n.size){
							for (int l = 0; l<n.size; l++){
								for (int m =0; m<c->lnra[k].size;m++){
									if (c->lnra[k].members[m]==n.members[l]){		
										common++;
										break;
									}
								}
							}

							if (common == n.size){
								alreadyin = 1;
								break;
							}
						}
					}

					if (alreadyin == 0) {				
						c->lracard++;
						c->lra[c->lracard-1].size=rprime.size;
						c->lra[c->lracard-1].members=(int*)malloc(rprime.size*sizeof(int));
						for (int k=0; k<rprime.size; k++)
							c->lra[c->lracard-1].members[k]=rprime.members[k];
						sizeofnext++;
						nextLevel[sizeofnext-1].size=rprime.size;
						nextLevel[sizeofnext-1].members=(int*)malloc(rprime.size*sizeof(int));
						for (int k=0; k<rprime.size; k++)
							nextLevel[sizeofnext-1].members[k]=rprime.members[k];			
						c->lnracard++;
						c->lnra[c->lnracard-1].size=n.size;
						c->lnra[c->lnracard-1].members=(int*)malloc(n.size*sizeof(int));
						for (int k=0; k<n.size; k++)
							c->lnra[c->lnracard-1].members[k]=n.members[k];
						c->assoc[c->lracard*2-2]=rprime;
						c->assoc[c->lracard*2-1]=n;
					}
				}
			}			
		}
		lastLevel=nextLevel;
		nextLevel=NULL;
		nextLevel=(pointset*)malloc(c->nrep*c->nrep*c->nrep*c->nrep*sizeof(pointset));
		sizeoflast=sizeofnext;
		sizeofnext=0;
	}

	c->lracomp=(pointset*)malloc(4*c->nrep*c->nrep*c->nrep*c->nrep*sizeof(pointset));
	c->lnracomp=(pointset*)malloc(4*c->nrep*c->nrep*c->nrep*c->nrep*sizeof(pointset));
	c->lracompcard=1;
	c->lnracompcard=1;
	free(nextLevel);
	nextLevel=(pointset*)malloc(4*c->nrep*c->nrep*c->nrep*c->nrep*sizeof(pointset));
	s.size=0;
	free(lastLevel);						
	lastLevel=(pointset*)malloc(4*c->nrep*c->nrep*c->nrepincomp*c->nrepincomp*sizeof(pointset));
	lastLevel[0]=s;
	sizeoflast=1;
	sizeofnext=0;
	c->assoccomp=(pointset*)malloc(2*c->nrep*c->nrep*c->nrepincomp*c->nrepincomp*sizeof(pointset));
	c->lracomp[0]=s;
	c->lnracomp[0]=s;
	c->assoccomp[0]=s;
	c->assoccomp[1]=s;
	
	while (sizeoflast!=0){								
		for (int i=0; i<sizeoflast;i++){					
			pointset r;
			r.size = lastLevel[i].size;

			if(r.size!=0)
				r.members=(int*)malloc(r.size*sizeof(int));

			for (int j=0; j< r.size;j++)
				r.members[j]=lastLevel[i].members[j];						
			for (int j=0; j<c->nrepincomp; j++){					
				pointset rprime;
				rprime.size=0;
				int alreadyin = 0;
				rprime.members=(int*)malloc((r.size+1)*sizeof(int));

				for (int k=0;k<r.size;k++){
					rprime.size++;
					rprime.members[k]=r.members[k];
					if (rprime.members[k]==c->complementtc[j]){
						alreadyin = 1;
						break;
					}
				}

				if (alreadyin==0) {							
					rprime.size++;							
					rprime.members[rprime.size-1]=c->complementtc[j];		
					pointset n;
					n.members=(int*)malloc(c->nrep*sizeof(int));								
					n.size=0;
					for (int l=0; l<c->nrep; l++){	
						for (int k=0;k<rprime.size;k++){

							if(g->matrix[rprime.members[k]+g->size*c->tc[l]]==1){
								int alreadyin=0;			
								for (int m=0; m<n.size;m++){
									if (n.members[m]==c->tc[l])
										alreadyin=1;
								}
								if (alreadyin==0){			
									n.size++;
									n.members[n.size-1]=c->tc[l];
									break;
								}
							}
						}
					}

					int alreadyin = 0;
					
					if (n.size==0)
						alreadyin = 1;
					for (int k=0; k<c->lnracompcard; k++) {
						if (c->lnracomp[k].size==n.size){
							int common = 0;
							for (int l = 0; l<n.size; l++){
								for (int m =0; m<c->lnracomp[k].size;m++){				
									if (c->lnracomp[k].members[m]==n.members[l]){		
										common++;
										break;
									}
								}
							}

							if (common == n.size){							
								alreadyin = 1;
								break;
							}
						}
					}
					
					if (alreadyin == 0) {
						c->lracompcard++;
						c->lracomp[c->lracompcard-1]=rprime;
						/*c->lracomp[c->lracompcard-1].size=rprime.size;
						c->lracomp[c->lracompcard-1].members=(int*)malloc(rprime.size*sizeof(int));
						for (int k=0; k<rprime.size; k++)
							c->lracomp[c->lracompcard-1].members[k]=rprime.members[k];*/
						sizeofnext++;
						nextLevel[sizeofnext-1].size=rprime.size;
						nextLevel[sizeofnext-1].members=(int*)malloc(rprime.size*sizeof(int));
						for (int k=0; k<rprime.size; k++)
							nextLevel[sizeofnext-1].members[k]=rprime.members[k];
						c->lnracompcard++;
						c->lnracomp[c->lnracompcard-1].size=n.size;
						c->lnracomp[c->lnracompcard-1].members=(int*)malloc(n.size*sizeof(int));
						for (int k=0; k<n.size; k++)
							c->lnracomp[c->lnracompcard-1].members[k]=n.members[k];
						c->assoccomp[c->lracompcard*2-2]=rprime;
						c->assoccomp[c->lracompcard*2-1]=n;
						
					}
				}
			}			
		}
		lastLevel=nextLevel;											
		nextLevel=NULL;			
		nextLevel=(pointset*)malloc(c->nrep*c->nrep*c->nrep*c->nrep*sizeof(pointset));
		sizeoflast=sizeofnext;
		sizeofnext=0;
	}

	return EXIT_SUCCESS;
}


int thirdpreprocess (cutdata* c, graph* g){

	c->m=(pointset*)malloc(c->lracard*c->nrep*sizeof(pointset));

	for (int i=0;i<c->nrep;i++){
		int v = c->tc[i];
		for (int j=0;j<c->lracard;j++){
			pointset r = c->lra[j];
			pointset rprime;
			rprime.size=r.size;
			rprime.members=(int*)malloc((r.size+1)*sizeof(int));
			int alreadyin = 0;
			for (int k=0;k<r.size;k++){
				rprime.members[k]=r.members[k];
				if (r.members[k]==v)
					alreadyin=1;
			}

			if (alreadyin==0){
				rprime.size++;
				rprime.members[r.size]=v;
			}

			pointset n;
			n.size=0;
			n.members = (int*)malloc(c->nrepincomp*sizeof(pointset));

			for (int k=0;k<c->nrepincomp;k++){
				int neighboor=0;
				for (int l=0;l<rprime.size;l++){
					if (g->matrix[rprime.members[l]*g->size+c->complementtc[k]]==1){
						neighboor=1;
						break;
					}
				}

				if (neighboor==1){
					n.size++;
					n.members[n.size-1]=c->complementtc[k];
				}
			}

			for (int k=0; k<c->lnracard;k++){
				if (c->lnra[k].size==n.size){
					int common = 0;
					for (int l=0; l<n.size;l++){
						for (int m=0; m<c->lnra[k].size; m++){
							if (n.members[l]==c->lnra[k].members[m])
								common++;
						}
					}

					if (common==n.size){
						c->m[j*c->nrep+i]=c->assoc[2*k];
						break;
					}
					else {
						pointset s;
						s.size=0;
						c->m[j*c->nrep+i]=s;
					}
				}
			}
		}
	}

	c->mcomp=(pointset*)malloc(c->lracompcard*c->nrepincomp*sizeof(pointset));

	for (int i=0;i<c->nrepincomp;i++){
		int v = c->complementtc[i];
		for (int j=0;j<c->lracompcard;j++){
			pointset r = c->lracomp[j];
			pointset rprime;
			rprime.size=r.size;
			rprime.members=(int*)malloc((r.size+1)*sizeof(int));
			int alreadyin = 0;
			for (int k=0;k<r.size;k++){
				rprime.members[k]=r.members[k];
				if (r.members[k]==v)
					alreadyin=1;
			}

			if (alreadyin==0){
				rprime.size++;
				rprime.members[r.size]=v;
			}

			pointset n;
			n.size=0;
			n.members = (int*)malloc(c->nrep*sizeof(pointset));

			for (int k=0;k<c->nrep;k++){
				int neighboor=0;
				for (int l=0;l<rprime.size;l++){
					if (g->matrix[rprime.members[l]+g->size*c->tc[k]]==1){
						neighboor=1;
						break;
					}
				}

				if (neighboor==1){
					n.size++;
					n.members[n.size-1]=c->tc[k];
				}
			}

			for (int k=0; k<c->lnracompcard;k++){
				if (c->lnracomp[k].size==n.size){
					int common = 0;
					for (int l=0; l<n.size;l++){
						for (int m=0; m<c->lnracomp[k].size; m++){
							if (n.members[l]==c->lnracomp[k].members[m])
								common++;
						}
					}

					if (common==n.size){
						c->mcomp[j*c->nrepincomp+i]=c->assoccomp[2*k];
						break;
					}
					else {
						pointset s;
						s.size=0;
						c->mcomp[j*c->nrepincomp+i]=s;
					}
				}
			}
		}
	}

	return EXIT_SUCCESS;
}

pointset toplevelalgorithm (dectree* t, graph* g){

	if ((t->right==NULL)||(t->left==NULL)){
		pointset s;
		s.size=-1;
		return s;
	}

	int nnodes = getnumberofnodes(t)-1;
	tocompute=0;
	nodestocompute=(dectree**)malloc(nnodes*sizeof(dectree*));
	if (t->right!=NULL)
		fillThevoid(t->right,g);

	if (t->left!=NULL)
		fillThevoid(t->left, g);
	int i=0;
	while (tocompute!=0){
		if (nodestocompute[i]->computed==0)
			stepalgorithm(nodestocompute[i],g);
		i=(i+1)%nnodes;
	}
	cutdata c1 = t->left->c;
	cutdata c2= t->right->c;
	int size=-1;	
	int amax=-1;
	int acmax=-1;
	int bmax=-1;
	int bcmax=-1;

	for (int i=0;i<c1.lracard;i++){
		for (int j=0;j<c2.lracard;j++){
			pointset p = c1.lra[i];
			pointset q = c2.lra[j];
			pointset po;
			po.size=0;
			po.members= (int*)malloc(c2.lracompcard*sizeof(int));

			for (int l=0; l<p.size; l++){
				int z=0;
				int y=0;
				int x=0;
				for (int m=0; m<c2.nacomp;m++){
					if (c2.pointtorepincomp[2*m]==p.members[l]){
						z=c2.pointtorepincomp[2*m+1];
						break;
					}
				} 
				for (int m=0; m<c2.nrepincomp; m++){
					if (c2.complementtc[m]==z){
						x=m;
						break;
					}
				}
				for (int m=0;m<c2.lracompcard;m++){
					int common=0;
					if (c2.lracomp[m].size==po.size){
						for (int n=0; n<po.size; n++){
							for (int o=0; o<c2.lracomp[m].size; o++){
								if (c2.lracomp[m].members[o]==po.members[n]){
									common++;
									break;
								}
							}
						}
						if (common==po.size){
							y = m;
							break;
						}
					}
				}

				po = c2.mcomp[y*c2.nrepincomp+x];
			}

			int pindex = 0;
			for (int l=0; l<c2.lracompcard; l++){
				int common=0;
				if (c2.lracomp[l].size==po.size){
					for (int m=0; m<po.size; m++){
						for (int o=0; o<c2.lracomp[l].size; o++){
							if (c2.lracomp[l].members[o]==po.members[m]){
								common++;
								break;
							}
						}
					}

					if (common ==po.size){
						pindex= l;
						break;
					}
				}
			}

			pointset qo;
			qo.size=0;
			qo.members= (int*)malloc(c1.lracompcard*sizeof(int));

			for (int l=0; l<q.size; l++){
				int z=0;
				int y=0;
				int x=0;
				for (int m=0; m<c1.nacomp;m++){
					if (c1.pointtorepincomp[2*m]==q.members[l]){
						z=c1.pointtorepincomp[2*m+1];
						break;
					}
				} 
				for (int m=0; m<c1.nrepincomp; m++){
					if (c1.complementtc[m]==z){
						x=m;
						break;
					}
				}
				for (int m=0;m<c1.lracompcard;m++){
					int common=0;
					if (c1.lracomp[m].size==qo.size){
						for (int n=0; n<qo.size; n++){
							for (int o=0; o<c1.lracomp[m].size; o++){
								if (c1.lracomp[m].members[o]==qo.members[n]){
									common++;
									break;
								}
							}
						}
						if (common==qo.size){
							y = m;
							break;
						}
					}
				}

				qo = c1.mcomp[y*c1.nrepincomp+x];
			}

			int qindex = 0;
			for (int l=0; l<c1.lracompcard; l++){
				int common=0;
				if (c1.lracomp[l].size==qo.size){
					for (int m=0; m<qo.size; m++){
						for (int o=0; o<c1.lracomp[l].size; o++){
							if (c1.lracomp[l].members[o]==qo.members[m]){
								common++;
								break;
							}
						}
					}

					if (common ==qo.size){
						qindex= l;
						break;
					}
				}
			}

			if (size==-1){
				if ((c1.tab[i*c1.lracompcard+qindex]!=-1)&&(c2.tab[j*c2.lracompcard+pindex]!=-1)){
					size = c1.tab[i*c1.lracompcard+qindex]+c2.tab[j*c2.lracompcard+pindex];
					amax = i;
					acmax = qindex;
					bmax= j;
					bcmax = pindex;
				}
			}
			else {
				if ((c1.tab[i*c1.lracompcard+qindex]!=-1)&&(c2.tab[j*c2.lracompcard+pindex]!=-1)){
					if (size > c1.tab[i*c1.lracompcard+qindex]+c2.tab[j*c2.lracompcard+pindex]){
						size = c1.tab[i*c1.lracompcard+qindex]+c2.tab[j*c2.lracompcard+pindex];
						amax = i;
						acmax = qindex;
						bmax= j;
						bcmax = pindex;
					}
				}
			}			
		}
	}
/*	printf("composant=");
	for (int i=0; i<g->size; i++)
		printf("(%d, %d), ",g->pos[2*i], g->pos[2*i+1]);
	printf("\n");
	printf("C1=\n");
	for (int i=0; i<c1.lracard; i++){
		for (int j=0; j<c1.lracompcard; j++)
			printf("%d ",c1.tab[i*c1.lracompcard+j]);
		printf("\n");
	}
	printf("C2=\n");
	for (int i=0; i<c2.lracard; i++){
		for (int j=0; j<c2.lracompcard; j++)
			printf("%d ", c2.tab[i*c2.lracompcard+j]);
		printf("\n");
	}*/
	pointset sol;
	sol.size=c1.tab[amax*c1.lracompcard+acmax];
	sol.members= (int*)malloc(size*sizeof(int));
	sol = computeDS (t->left, c1.tab[amax*c1.lracompcard+acmax], amax, acmax);
	pointset tmp;
	tmp.size=c2.tab[bmax*c2.lracompcard+bcmax];
	tmp.members= (int*)malloc(tmp.size*sizeof(int));
	tmp = computeDS (t->right, c2.tab[bmax*c2.lracompcard+bcmax], bmax, bcmax);

	for (int i=0; i<tmp.size; i++){
		sol.members[sol.size]=tmp.members[i];
		sol.size++;
	}
	return sol;
}


int stepalgorithm (dectree* t, graph* g){

	if ((t->left!=0)&&(t->right!=NULL)){
		if ((t->right->computed==0)||(t->left->computed==0))
			return 0;
	}

	pointset s;
	s.size=0;
	pasta sp;
	sp.where=s;
	sp.howmany=0;
	pastabox spb;
	spb.inleft=sp;
	spb.inright=sp;
	t->c = cutThatTree (g, t);
	firstpreprocess (g,&(t->c));
	secondpreprocess (&(t->c), g);
	thirdpreprocess (&(t->c), g);

	t->c.tab = (int*)malloc(t->c.lracard*t->c.lracompcard*sizeof(int));
	t->c.box = (pastabox*)malloc(t->c.lracard*t->c.lracompcard*sizeof(pastabox));
	for (int i= 0; i<t->c.lracard*t->c.lracompcard; i++){
		t->c.tab[i]=-1;
		t->c.box[i]=spb;
	}

	if ((t->left!=NULL)&&(t->right!=NULL)){
		for (int i=0; i< t->left->c.lracard; i++){
			for (int j= 0;j < t->right->c.lracard; j++){
				for (int k=0; k< t->c.lracompcard; k++){
					pointset ra = t->left->c.lra[i];
					pointset rb = t->right->c.lra[j];
					pointset rwc = t->c.lracomp[k];
					pointset ua;
					ua.size = rb.size;
					ua.members = (int*) malloc ((rb.size+rwc.size)*sizeof(int));
					for (int l = 0; l<rb.size;l++)
						ua.members[l]=rb.members[l];
					for (int l=0; l<rwc.size;l++){
						int alreadyin = 0;
						for (int m = 0; m<ua.size; m++){
							if (ua.members[m]==rwc.members[l]){
								alreadyin = 1;
								break;
							}
						}
						if (alreadyin==0){
							ua.size++;
							ua.members[ua.size-1]=rwc.members[l];
						}
					}
					pointset rac;
					rac.size=0;
					rac.members= (int*)malloc(t->left->c.lracompcard*sizeof(int));

					for (int l=0; l<ua.size; l++){
						int z=0;
						int y=0;
						int x=0;
						for (int m=0; m<t->left->c.nacomp;m++){
							if (t->left->c.pointtorepincomp[2*m]==ua.members[l]){
								z=t->left->c.pointtorepincomp[2*m+1];
								break;
							}
						} 
						for (int m=0; m<t->left->c.nrepincomp; m++){
							if (t->left->c.complementtc[m]==z){
								x=m;
								break;
							}
						}
						for (int m=0;m<t->left->c.lracompcard;m++){
							int common=0;
							if (t->left->c.lracomp[m].size==rac.size){
								for (int n=0; n<rac.size; n++){
									for (int o=0; o<t->left->c.lracomp[m].size; o++){
										if (t->left->c.lracomp[m].members[o]==rac.members[n]){
											common++;
											break;
										}
									}
								}
								if (common==rac.size){
									y = m;
									break;
								}
							}
						}
						rac = t->left->c.mcomp[y*t->left->c.nrepincomp+x];
					}
					
					pointset ub;
					ub.size = ra.size;
					ub.members = (int*) malloc ((ra.size+rwc.size)*sizeof(int));
					for (int l = 0; l<ra.size;l++)
						ub.members[l]=ra.members[l];
					for (int l=0; l<rwc.size;l++){
						int alreadyin = 0;
						for (int m = 0; m<ub.size; m++){
							if (ub.members[m]==rwc.members[l]){
								alreadyin = 1;
								break;
							}
						}
						if (alreadyin==0){
							ub.size++;
							ub.members[ub.size-1]=rwc.members[l];
						}
					}
					pointset rbc;
					rbc.size=0;
					rbc.members= (int*)malloc(t->right->c.lracompcard*sizeof(int));

					for (int l=0; l<ub.size; l++){
						int z=0;
						int y=0;
						int x=0;
						for (int m=0; m<t->right->c.nacomp;m++){
							if (t->right->c.pointtorepincomp[2*m]==ub.members[l]){
								z=t->right->c.pointtorepincomp[2*m+1];
								break;
							}
						} 
						for (int m=0;m<t->right->c.nrepincomp;m++){
							if (t->right->c.complementtc[m]==z){
								x=m;
								break;
							}
						}
						for (int m=0;m<t->right->c.lracompcard;m++){
							int common=0;
							if (t->right->c.lracomp[m].size==rbc.size){
								for (int n=0; n<rbc.size; n++){
									for (int o=0; o<t->right->c.lracomp[m].size; o++){
										if (t->right->c.lracomp[m].members[o]==rbc.members[n]){
											common++;
											break;
										}
									}
								}
								if (common==rbc.size){
									y = m;
									break;
								}
							}
						}
						rbc = t->right->c.mcomp[y*t->right->c.nrepincomp+x];
					}

					pointset uw;
					uw.size = ra.size;
					uw.members = (int*) malloc ((ra.size+rb.size)*sizeof(int));
					for (int l = 0; l<ra.size;l++)
						uw.members[l]=ra.members[l];
					for (int l=0; l<rb.size;l++){
						int alreadyin = 0;
						for (int m = 0; m<uw.size; m++){
							if (uw.members[m]==rb.members[l]){
								alreadyin = 1;
								break;
							}
						}
						if (alreadyin==0){
							uw.size++;
							uw.members[uw.size-1]=rb.members[l];
						}
					}
				
					pointset rw;
					rw.size=0;
					rw.members= (int*)malloc(t->c.lracard*sizeof(int));

					for (int l=0; l<uw.size; l++){
						int z=0;
						int y=0;
						int x=0;
						for (int m=0; m<t->c.na;m++){
							if (t->c.pointtorep[2*m]==uw.members[l]){
								z=t->c.pointtorep[2*m+1];
								break;
							}
						} 
						for (int m=0;m<t->c.nrep;m++){
							if (t->c.tc[m]==z){
								x=m;
								break;
							}
						}
						for (int m=0;m<t->c.lracard;m++){
							int common=0;
							if (t->c.lra[m].size==rw.size){
								for (int n=0; n<rw.size; n++){
									for (int o=0; o<t->c.lra[m].size; o++){
										if (t->c.lra[m].members[o]==rw.members[n]){
											common++;
											break;
										}
									}
								}
								if (common==rw.size){
									y = m;
									break;
								}
							}
						}
						rw = t->c.m[y*t->c.nrep+x];
						
					}

					int ain=0;
					int bin=0;
					int win=0;
					int acin=0;
					int bcin=0;
					int wcin=0;
							
					for (int l=0;l<t->left->c.lracard;l++){
						if (t->left->c.lra[l].size==ra.size){
							int common=0;
							for (int m=0;m<t->left->c.lra[l].size;m++){
								for (int n=0;n<ra.size; n++){
									if (t->left->c.lra[l].members[m]==ra.members[n]){
										common ++;
										break;
									}
								}
							}
							if (common==ra.size){
								ain=l;
								break;
							}
						}
					}

					for (int l=0;l<t->right->c.lracard;l++){
						if (t->right->c.lra[l].size==rb.size){
							int common=0;
							for (int m=0;m<t->right->c.lra[l].size;m++){
								for (int n=0;n<rb.size; n++){
									if (t->right->c.lra[l].members[m]==rb.members[n]){
										common ++;
										break;
									}
								}
							}
							if (common==rb.size){
								bin=l;
								break;
							}
						}
					}


					for (int l=0;l<t->c.lracard;l++){
						if (t->c.lra[l].size==rw.size){
							int common=0;
							for (int m=0;m<t->c.lra[l].size;m++){
								for (int n=0;n<rw.size; n++){
									if (t->c.lra[l].members[m]==rw.members[n]){
										common ++;
										break;
									}
								}
							}
							if (common==rw.size){
								win=l;
								break;
							}
						}
					}

					for (int l=0;l<t->left->c.lracompcard;l++){
						if (t->left->c.lracomp[l].size==rac.size){
							int common=0;
							for (int m=0;m<t->left->c.lracomp[l].size;m++){
								for (int n=0;n<rac.size; n++){
									if (t->left->c.lracomp[l].members[m]==rac.members[n]){
										common ++;
										break;
									}
								}
							}
							if (common==rac.size){
								acin=l;
								break;
							}
						}
					}


					for (int l=0;l<t->right->c.lracompcard;l++){
						if (t->right->c.lracomp[l].size==rbc.size){
							int common=0;
							for (int m=0;m<t->right->c.lracomp[l].size;m++){
								for (int n=0;n<rbc.size; n++){
									if (t->right->c.lracomp[l].members[m]==rbc.members[n]){
										common ++;
										break;
									}
								}
							}
							if (common==rbc.size){
								bcin=l;
								break;
							}
						}
					}


					for (int l=0;l<t->c.lracompcard;l++){
						if (t->c.lracomp[l].size==rwc.size){
							int common=0;
							for (int m=0;m<t->c.lracomp[l].size;m++){
								for (int n=0;n<rwc.size; n++){
									if (t->c.lracomp[l].members[m]==rwc.members[n]){
										common ++;
										break;
									}
								}
							}
							if (common==rwc.size){
								wcin=l;
								break;
							}
						}
					}
					//printf("ain=%d, acin=%d, bin=%d, bcin=%d, win=%d, wcin=%d, sizea=%d, sizeb=%d, sizew=%d\n", ain, acin, bin, bcin, win, wcin, t->left->c.tab[ain*t->left->c.lracompcard+acin], t->right->c.tab[bin*t->right->c.lracompcard+bcin], t->c.tab[win*t->c.lracompcard+wcin]);								
					if ((t->left->c.tab[ain*t->left->c.lracompcard+acin]!=-1)&&(t->right->c.tab[bin*t->right->c.lracompcard+bcin]!=-1)){
						if ((t->c.tab[win*t->c.lracompcard+wcin]==-1)||(t->c.tab[win*t->c.lracompcard+wcin]>t->left->c.tab[ain*t->left->c.lracompcard+acin]+t->right->c.tab[bin*t->right->c.lracompcard+bcin])){
			
							t->c.tab[win*t->c.lracompcard+wcin]=t->left->c.tab[ain*t->left->c.lracompcard+acin]+t->right->c.tab[bin*t->right->c.lracompcard+bcin];

							pasta pasta1;
							pasta1.where=t->left->c.lra[ain];
							pasta1.howmany=t->left->c.tab[ain*t->left->c.lracompcard+acin];

							pasta pasta2;
							pasta2.where=t->right->c.lra[bin];
							pasta2.howmany=t->right->c.tab[bin*t->right->c.lracompcard+bcin];

							pastabox pastab;
							pastab.inleft=pasta1;
							pastab.inright=pasta2;
							t->c.box[win*t->c.lracompcard+wcin]=pastab;
						}
					}					

	
				}
		
			}
			
		}
	}
	else {
		t->c.tab[0]=-1;
		t->c.tab[1]=0;
		t->c.tab[2]=1;
		t->c.tab[3]=1;
	}
	t->computed=1;
	tocompute--;
/*	printf("a=");
	for (int i=0; i<t->c.na; i++)
		printf("%d, ",t->c.a[i]);
	printf("\n");
	printf("matrix=\n");
	for (int i=0; i<t->c.na; i++){
		for (int j=0; j<t->c.nacomp; j++)
			printf("%d ",t->c.matrixrevisited[i*t->c.nacomp+j]);
		printf("\n");
	}
	printf("pointtorepincomp=\n");
	for (int i=0; i<t->c.nacomp; i++)
		printf("%d -> %d\n",t->c.pointtorepincomp[2*i], t->c.pointtorepincomp[2*i+1]);
	printf("lra=\n");
	for (int i=0; i<t->c.lracard; i++){
		for (int j=0; j<t->c.lra[i].size; j++)
			printf("%d ",t->c.lra[i].members[j]);
		printf("\n");
	}
	printf("lracomp=\n");
	for (int i=0; i<t->c.lracompcard; i++){
		for (int j=0; j<t->c.lracomp[i].size; j++)
			printf("%d ",t->c.lracomp[i].members[j]);
		printf("\n");
	}
	printf("mcomp=\n");
	for (int i=0; i<t->c.lracompcard; i++){
		for (int j=0; j<t->c.nrepincomp; j++){
			printf("(");
			for (int k=0; k<t->c.mcomp[i*t->c.nrepincomp+j].size; k++)
				printf("%d, ",t->c.mcomp[i*t->c.nrepincomp+j].members[k]);
			printf(") ");
		}
		printf("\n");
	}
	printf("tab=\n");
	for (int i=0; i<t->c.lracard; i++){
		for (int j=0; j<t->c.lracompcard; j++)
			printf("%d ", t->c.tab[i*t->c.lracompcard+j]);
		printf("\n");
	}*/
	return EXIT_SUCCESS;

}


pointset computeDS (dectree* t, int much, int a, int ac){
	pointset p;
	p.size=0;
	p.members=(int*)malloc((t->c.na)*sizeof(int));

	if ((t->left==NULL)&&(t->left==NULL)&&(much==1)){
		p.size++;
		p.members[p.size-1]=t->c.tc[0];
	}

	else if (much!=0){

		pastabox pbleft = t->c.box[a*t->c.lracompcard+ac];
		pasta p1 = pbleft.inleft;
		pasta p2 = pbleft.inright;
		pointset pointleft = p1.where;
		pointset pointright = p2.where;
		pointset complem = t->c.lracomp[ac];

		int anext=0;
		int acnext=0;
		for (int i=0;i<t->left->c.lracard;i++){
			if (t->left->c.lra[i].size==pointleft.size){
				int common = 0;
				for (int j=0;j<t->left->c.lra[i].size;j++){
					for (int k=0;k<pointleft.size;k++){
						if (t->left->c.lra[i].members[j]==pointleft.members[k]){
							common++;
							break;
						}
					}
				}
				if (common ==pointleft.size){
					anext=i;
					break;
				}
			}
		}

		
		pointset ua;
		ua.size = pointright.size;
		ua.members = (int*) malloc ((pointright.size+complem.size)*sizeof(int));
		for (int l = 0; l<pointright.size;l++)
			ua.members[l]=pointright.members[l];
		for (int l=0; l<complem.size;l++){
			int alreadyin = 0;
			for (int m = 0; m<ua.size; m++){
				if (ua.members[m]==complem.members[l]){
					alreadyin = 1;
					break;
				}
			}
			if (alreadyin==0){
				ua.size++;
				ua.members[ua.size-1]=complem.members[l];
			}
		}
		pointset rac;
		rac.size=0;
		rac.members= (int*)malloc(t->left->c.lracompcard*sizeof(int));

		for (int l=0; l<ua.size; l++){
			int z=0;
			int y=0;
			int x=0;
			for (int m=0; m<t->left->c.nacomp;m++){
				if (t->left->c.pointtorepincomp[2*m]==ua.members[l]){
					z=t->left->c.pointtorepincomp[2*m+1];
					break;
				}
			} 
			for (int m=0;m<t->left->c.nrepincomp;m++){
				if (t->left->c.complementtc[m]==z){
					x=m;
					break;
				}
			}
			for (int m=0;m<t->left->c.lracompcard;m++){
				int common=0;
				if (t->left->c.lracomp[m].size==rac.size){
					for (int n=0; n<rac.size; n++){
						for (int o=0; o<t->left->c.lracomp[m].size; o++){
							if (t->left->c.lracomp[m].members[o]==rac.members[n]){
								common++;
								break;
							}
						}
					}
					if (common==rac.size){
						y = m;
						break;
					}
				}
			}
			rac = t->left->c.mcomp[y*t->left->c.nrepincomp+x];
		}

		for (int i=0;i<t->left->c.lracompcard;i++){
			if (t->left->c.lracomp[i].size==rac.size){
				int common = 0;
				for (int j=0;j<t->left->c.lracomp[i].size;j++){
					for (int k=0;k<rac.size;k++){
						if (t->left->c.lracomp[i].members[j]==rac.members[k]){
							common++;
							break;
						}
					}
				}
				if (common ==rac.size){
					acnext=i;
					break;
				}
			}
		}

		int bnext=0;
		int bcnext=0;
		for (int i=0;i<t->right->c.lracard;i++){
			if (t->right->c.lra[i].size==pointright.size){
				int common = 0;
				for (int j=0;j<t->right->c.lra[i].size;j++){
					for (int k=0;k<pointright.size;k++){
						if (t->right->c.lra[i].members[j]==pointright.members[k]){
							common++;
							break;
						}
					}
				}
				if (common ==pointright.size){
					bnext=i;
					break;
				}
			}
		}
		
		pointset ub;
		ub.size = pointleft.size;
		ub.members = (int*) malloc ((pointleft.size+complem.size)*sizeof(int));
		for (int l = 0; l<pointleft.size;l++)
			ub.members[l]=pointleft.members[l];
		for (int l=0; l<complem.size;l++){
			int alreadyin = 0;
			for (int m = 0; m<ub.size; m++){
				if (ub.members[m]==complem.members[l]){
					alreadyin = 1;
					break;
				}
			}
			if (alreadyin==0){
				ub.size++;
				ub.members[ub.size-1]=complem.members[l];
			}
		}
		pointset rbc;
		rbc.size=0;
		rbc.members= (int*)malloc(t->right->c.lracompcard*sizeof(int));

		for (int l=0; l<ub.size; l++){
			int z=0;
			int y=0;
			int x=0;
			for (int m=0; m<t->right->c.nacomp;m++){
				if (t->right->c.pointtorepincomp[2*m]==ub.members[l]){
					z=t->right->c.pointtorepincomp[2*m+1];
					break;
				}
			} 
			for (int m=0;m<t->right->c.nrepincomp;m++){
				if (t->right->c.complementtc[m]==z){
					x=m;
					break;
				}
			}
			for (int m=0;m<t->right->c.lracompcard;m++){
				int common=0;
				if (t->right->c.lracomp[m].size==rbc.size){
					for (int n=0; n<rbc.size; n++){
						for (int o=0; o<t->right->c.lracomp[m].size; o++){
							if (t->right->c.lracomp[m].members[o]==rbc.members[n]){
								common++;
								break;
							}
						}
					}
					if (common==rbc.size){
						y = m;
						break;
					}
				}
			}
			rbc = t->right->c.mcomp[y*t->right->c.nrepincomp+x];
		}

		for (int i=0;i<t->right->c.lracompcard;i++){
			if (t->right->c.lracomp[i].size==rbc.size){
				int common = 0;
				for (int j=0;j<t->right->c.lracomp[i].size;j++){
					for (int k=0;k<rbc.size;k++){
						if (t->right->c.lracomp[i].members[j]==rbc.members[k]){
							common++;
							break;
						}
					}
				}
				if (common ==rbc.size){
					bcnext=i;
					break;
				}
			}
		}
		pointset ptmp;
		ptmp.size=p1.howmany;
		ptmp.members=(int*)malloc(p1.howmany*sizeof(int));
		ptmp = computeDS (t->left, p1.howmany, anext, acnext);
		for (int i=0; i<ptmp.size; i++){
			p.size++;
			p.members[p.size-1]=ptmp.members[i];
		}
		ptmp.size= p2.howmany;
		ptmp.members=(int*)malloc(p2.howmany*sizeof(int));
		ptmp= computeDS(t->right, p2.howmany, bnext, bcnext);
		for (int i=0; i<ptmp.size; i++){
			p.size++;
			p.members[p.size-1]=ptmp.members[i];
		}
	}
	return p;
}

int getBW (dectree* t, graph* g){
	int bwmax=-1;
	if ((t==NULL)||(t->right==NULL)||(t->left==NULL))
		bwmax=2;
	else {
		t->c = cutThatTree (g, t);

		firstpreprocess (g,&(t->c));
		secondpreprocess (&(t->c), g);

		
		int p  = getBW(t->left, g);
		int q  = getBW(t->right, g);
		

		if (p>bwmax)
			bwmax=p;
		if (q>bwmax)
			bwmax=q;
		if (t->c.lracard>bwmax)
			bwmax=t->c.lracard;

	}
	return bwmax;
}
