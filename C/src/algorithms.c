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


// This function takes a graph, a decomposition tree and and int stating if we're to cut the left or right son of the tree t to create the
// cut in the decomposition graph
cutdata cutThatTree (graph g, dectree t, int choiceofson){
	cutdata c;
	c.t=t;
	c.choiceofson=choiceofson;
	
	c.na=0;
	c.nacomp=0;
	c.a=NULL;
	c.acomp=NULL;

	if (choiceofson==0){						// a choiceofson of 0 means that we're cutting the decomposition tree along the left exiting
		c.na=getnumberofleaves (*(t.left));		// edge of the current node. So primary sub-graph is going to be the set of points
		c.a=(int*)malloc(c.na*sizeof(int));		// represented by leaves of the left son node
		getallleaves(*(t.left), c.a);

	}
	
	else {
		c.na=getnumberofleaves (*(t.right));	// a choice of son 1 means that we're cutting the tree along the right exiting edge
		c.a=(int*)malloc(c.na*sizeof(int));
		getallleaves (*(t.right),c.a);
	}

	c.nacomp=g.size-c.na;

	c.acomp=(int*)malloc(c.nacomp*sizeof(int));

	int i=0;
	int j=0;
	int k=0;
	for (i=0;i<g.size;i++){						// all other points will be in the complementary, the second sub-graph. That's why we		
		int ina=0;								// iterate on graph's size, to get all points of the graph and not just points that
		for (j=0;j<c.na;j++){					// are contained in sons of the current node
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



	// Now that we have determined the vertices we're interested in, we compute the adjacency matrix reduced to those vertices
	for (int i=0;i<c.na;i++){
		for (int j=0;j<c.nacomp;j++){

			c.matrixrevisited[i*c.nacomp+j]=g.matrix[c.a[i]*g.size+c.acomp[j]];

		}

	}
	return c;
}


/*
This first pre-process takes a graph, its decomposition tree and a cut of this decomposition tree as input and computes its equivalency 
classes, the output being the lists of representants of the equivalency classes of both primary and secondary sub-graphs, their numbers
and two lists which associate each point of a subgraph to the representant of its equivalency class 
*/
cutdata firstpreprocess(graph g,  cutdata c){

	int *twinsleft = (int*)malloc(2*c.na*c.na*sizeof(int));			// we have at most na squared twin pairs, in the case of a complete graph
	int *twinsright = (int*)malloc(2*c.nacomp*c.nacomp*sizeof(int));

	matchTwins(c ,twinsleft,twinsright);							// we compute the list of twins in each sub-graph. This will be a 
																	// sequency with each eventh vertice being twin of the following oddth
																	// vertice, relative to the induced partition of the graph

	c.pointtorep=(int*)malloc(2*c.na*sizeof(int));					// we're going to associate the representant of the equivalency class
	c.pointtorepincomp=(int*)malloc(2*c.nacomp*sizeof(int));		// to each vertice of each subgraph
	
	for (int i=0; i<2*c.na;i++){
		c.pointtorep[i]=-1;
	}

	for (int i=0;i<2*c.nacomp;i++){
		c.pointtorepincomp[i]=-1;
	}

	for (int i=0;i<c.na;i++){
		c.pointtorep[i*2]=c.a[i];

		for (int j=0;j<c.na*c.na;j++){
			int a = twinsleft[j*2];
			int b = twinsleft[j*2+1];								// for each vertice in primary sub-graph A, we try to find every
																	// occurence of it in the twins matrix
			if (a==i){												// in the twins matrix, vertices are designed by indexes in the sub-graph
				if (c.a[b]>c.pointtorep[i*2+1]){					// once we find an occurence of the vertice, we check if the vertice to
					c.pointtorep[i*2+1]=c.a[b];						// which he is a twin has an index (in the whole graph) greater than
				}													// the one registered in the current pointerToRep matrix, in which case
			}														// we save this one instead
			if (b==i){
				if (c.a[a]>c.pointtorep[i*2+1]){
					c.pointtorep[i*2+1]=c.a[a];
				}
			}
		}

		// at the end of those loops, the pointerToRep matrix associates each vertex to the representant of his equivalency class
		// which is the vertex which is a twin of it with the greatest index possible
		if (c.pointtorep[i*2+1]<c.a[i]){
			c.pointtorep[i*2+1]=c.a[i];								// but we've not treated the case of the vertex itself, which could	
		}															// be the one with the greatest index of its equivalency class
	}																// if that's the case, this vertex becomes its own representant

	// we do the exact same thing for the complementary of A, the secondary sub-graph
	for (int i=0;i<c.nacomp;i++){
		c.pointtorepincomp[i*2]=c.acomp[i];
		for (int j=0;j<c.nacomp*c.nacomp;j++){
			int a = twinsright[j*2];
			int b = twinsright[j*2+1];
		
			if (a==i){
				if (c.acomp[b]>c.pointtorepincomp[i*2+1]){
					c.pointtorepincomp[i*2+1]=c.acomp[b];
				}
			}
			if (b==i){
				if (c.acomp[a]>c.pointtorepincomp[i*2+1]){
					c.pointtorepincomp[i*2+1]=c.acomp[a];
				}
			}
		}

		if (c.pointtorepincomp[i*2+1]<c.acomp[i]){
			c.pointtorepincomp[i*2+1]=c.acomp[i];
		}
	}

	// now that we have computed equivalency class and associated each verted to the representant of its class, we list those representants
	c.tc=(int*)malloc(c.na*sizeof(int));
	for (int i=0;i<c.na;i++)
		c.tc[i]=-1;
	c.complementtc=(int*)malloc(c.nacomp*sizeof(int));
	int cursor = 0;
	for (int i=0;i<c.na;i++){
		int here=0;
		for (int j=0;j<c.na;j++){
			if (c.pointtorep[2*i+1]==c.tc[j]){
				here=1;										// for each representant in the pointerToRep matrix, we check if it has already
				break;										// been registered in the list of representants. If it has, no need to register
			}
		}
		if (here==0){										// if it has not, we register it and increments our cursor
			c.tc[cursor]=c.pointtorep[2*i+1];			
			cursor++;
		}
	}

	c.nrep=cursor;											// this cursor has counted all representatives and we store that



	// we do the same for secondary sub-graph
	cursor=0;
	for (int i=0;i<c.nacomp;i++){
		int here=0;
		for (int j=0;j<c.nacomp;j++){
			if (c.pointtorepincomp[2*i+1]==c.complementtc[j]){
				here=1;
				break;
			}
		}
		if (here==0){
			c.complementtc[cursor]=c.pointtorepincomp[2*i+1];
			cursor++;
		}
	}
	c.nrepincomp=cursor;



	/*
	 now that we have reduced the graph to equivalency classes we'd like to reduce the matrix to the representants of those equivalency
	classes in order to speed up future computations ever further
	*/


	int *tmp = (int*) malloc(c.nrep*c.nrepincomp*sizeof(int));
	for (int i=0;i<c.nrep;i++){
		for (int j=0;j<c.nrepincomp;j++){
			tmp[i*c.nrepincomp+j]=g.matrix[c.tc[i]*g.size+c.complementtc[j]];
		}
	}

	c.matrixrevisited=tmp;

	return c;

}



// This if the function of wonders! It computes twins in each subgraph by checking if their induced vectors in the adjacency matrix match
int matchTwins (cutdata c, int* twinsleft, int* twinsright){
	int i,j,k;
	int cursor=0;


	for(i=0;i<2*c.na*c.na;i++)
		twinsleft[i]=-1;

	for(i=0;i<2*c.nacomp*c.nacomp;i++)
		twinsright[i]=-1;

	for (i=0;i<c.na;i++){
		for (j=i+1;j<c.na;j++){					// we take each pair of vertices of one subgraph	
			int tw=1;	

			for (k=0;k<c.nacomp;k++){			// and we check coordinates by coordinates that they have the same induced vector
				if(c.matrixrevisited[i*c.nacomp+k]!=c.matrixrevisited[j*c.nacomp+k]){
				
					tw=0;						// upon finding a coordinate that doesn't match, we flag this pair as non-twins
					break;
				}
			}	
			
			if (tw==1){							// if no conflit has been detected, the pair is a twins pair and we write it down in the
				twinsleft[cursor]=i;			// twin matrix
				cursor++;
				twinsleft[cursor]=j;
				cursor++;
			}
			
		} 		
	}

	cursor = 0;

	// we do the same on secondary sub-graph
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



/* 
This is the second pre-processing function. We give it a tree and a cutdata on which the first preprocess has been carried on 
and it uses informations on representants to compute the list of representant sets, their respectives neighboorhoods and associations
between those two sets of data, to be stored in the cutdata
*/
cutdata secondpreprocess (cutdata c, graph g){

	c.lra=(pointset*)malloc(c.nrep*c.nrep*c.nrep*c.nrep*sizeof(pointset));
	c.lnra=(pointset*)malloc(c.nrep*c.nrep*c.nrep*c.nrep*sizeof(pointset));
	c.lracard=1;
	c.lnracard=1;
	pointset *nextLevel=(pointset*)malloc(c.nrep*c.nrep*c.nrep*c.nrep*sizeof(pointset));
	pointset s;
	s.size=0;						
	pointset *lastLevel=(pointset*)malloc(c.nrep*c.nrep*c.nrepincomp*c.nrepincomp*sizeof(pointset));
	lastLevel[0]=s;
	int sizeoflast=1;
	int sizeofnext=0;
	c.assoc=(pointset*)malloc(2*c.nrep*c.nrep*c.nrepincomp*c.nrepincomp*sizeof(pointset));

	c.lra[0]=s;
	c.lnra[0]=s;
	c.assoc[0]=s;
	c.assoc[1]=s;
	
	while (sizeoflast!=0){								// we initialize lastLevel as a set containing only an empty set of points
		for (int i=0; i<sizeoflast;i++){					// and we loop until lastLevel is an empty set itself
			pointset r;
			r.size = lastLevel[i].size;

			if(r.size!=0)
				r.members=(int*)malloc(r.size*sizeof(int));

			for (int j=0; j< r.size;j++)
				r.members[j]=lastLevel[i].members[j];						// for each set of points in lastLevel
			for (int j=0; j<c.nrep; j++){					// we try adding each representants of equivalency classes of primary subgraph	
				pointset rprime;


				rprime.size=0;
				int alreadyin = 0;

				rprime.members=(int*)malloc((r.size+1)*sizeof(int));

				for (int k=0;k<r.size;k++){
					rprime.size++;

					rprime.members[k]=r.members[k];
					if (rprime.members[k]==c.tc[j])
						alreadyin = 1;
				}

				if (alreadyin==0) {							// if the vertex was already in the pointset, we have nothing to add, so we pass
					rprime.size++;							// on to the next representant. Else, we add it to a pointset and compute the
															// neighboorhood of this new set of points


					rprime.members[rprime.size-1]=c.tc[j];		
					pointset n;
					n.members=(int*)malloc(c.nrepincomp*sizeof(int));								
					n.size=0;

					for (int k=0;k<rprime.size;k++){

						for (int l=0; l<c.nrepincomp; l++){	// for each point of this new set we get neighboorhoods

							if(g.matrix[rprime.members[k]*g.size + c.complementtc[l]]==1){
															// for each representant of equivalence classes in the secondary subgraph, we
															// check if they are neighboors of the point computed

								int alreadyin=0;			// we check that this point has not yet been added to the neighboorhood
								for (int m=0; m<n.size;m++){

									if (n.members[m]==c.complementtc[l]){
										alreadyin=1;
									}

								}

								if (alreadyin==0){			// else, we add it
									n.size++;

									n.members[n.size-1]=c.complementtc[l];

								}
							}
						}
					}


					// we then check if this neighboorhood is already in lnra
					int alreadyin = 0;
					
					if (n.size==0)
						alreadyin = 1;
					for (int k=0; k<c.lnracard; k++) {

						pointset tmp = c.lnra[k];
						int common = 0;
						for (int l = 0; l<n.size; l++){
							for (int m =0; m<tmp.size;m++){				// we count the number of vertices of n in each set of lnra
								if (tmp.members[m]==n.members[l]){		
									common++;
									break;
								}
							}
						}

						if (common == n.size){							// if this number is equal to the size of n, n is already in lnra
							alreadyin = 1;
							break;
						}
					}

					if (alreadyin == 0) {								// if n isn't in lnra, we going to add it and the rprime associated

				
						c.lracard++;
						c.lra[c.lracard-1]=rprime;

						sizeofnext++;

						nextLevel[sizeofnext-1].size=rprime.size;
						nextLevel[sizeofnext-1].members=(int*)malloc(rprime.size*sizeof(int));
						for (int y=0; y<rprime.size; y++)
							nextLevel[sizeofnext-1].members[y]=rprime.members[y];
				
						c.lnracard++;

						c.lnra[c.lnracard-1]=n;


						c.assoc[c.lracard*2-2]=rprime;
						c.assoc[c.lracard*2-1]=n;

					}
					

				}


			}
			
		}
		lastLevel=nextLevel;											// we loop back until no representating set and no neighboorhoods 

		nextLevel=NULL;			// are interesting enough to be added
		nextLevel=(pointset*)malloc(c.nrep*c.nrep*c.nrep*c.nrep*sizeof(pointset));
		sizeoflast=sizeofnext;

		sizeofnext=0;
	}

	c.lracomp=(pointset*)malloc(c.nrep*c.nrep*c.nrep*c.nrep*sizeof(pointset));
	c.lnracomp=(pointset*)malloc(c.nrep*c.nrep*c.nrep*c.nrep*sizeof(pointset));
	c.lracompcard=1;
	c.lnracompcard=1;
	free(nextLevel);
	nextLevel=(pointset*)malloc(c.nrep*c.nrep*c.nrep*c.nrep*sizeof(pointset));
	s.size=0;
	free(lastLevel);						
	lastLevel=(pointset*)malloc(c.nrep*c.nrep*c.nrepincomp*c.nrepincomp*sizeof(pointset));
	lastLevel[0]=s;
	sizeoflast=1;
	sizeofnext=0;
	c.assoccomp=(pointset*)malloc(2*c.nrep*c.nrep*c.nrepincomp*c.nrepincomp*sizeof(pointset));

	c.lracomp[0]=s;
	c.lnracomp[0]=s;
	c.assoccomp[0]=s;
	c.assoccomp[1]=s;
	
	while (sizeoflast!=0){								
		for (int i=0; i<sizeoflast;i++){					
			pointset r;
			r.size = lastLevel[i].size;

			if(r.size!=0)
				r.members=(int*)malloc(r.size*sizeof(int));

			for (int j=0; j< r.size;j++)
				r.members[j]=lastLevel[i].members[j];						
			for (int j=0; j<c.nrepincomp; j++){					
				pointset rprime;


				rprime.size=0;
				int alreadyin = 0;

				rprime.members=(int*)malloc((r.size+1)*sizeof(int));

				for (int k=0;k<r.size;k++){
					rprime.size++;

					rprime.members[k]=r.members[k];
					if (rprime.members[k]==c.complementtc[j])
						alreadyin = 1;
				}

				if (alreadyin==0) {							
					rprime.size++;							

					rprime.members[rprime.size-1]=c.complementtc[j];		
					pointset n;
					n.members=(int*)malloc(c.nrep*sizeof(int));								
					n.size=0;

					for (int k=0;k<rprime.size;k++){

						for (int l=0; l<c.nrep; l++){	
							if(g.matrix[rprime.members[k]+g.size*c.tc[l]]==1){
													

								int alreadyin=0;			
								for (int m=0; m<n.size;m++){

									if (n.members[m]==c.tc[l]){
										alreadyin=1;
									}

								}

								if (alreadyin==0){			
									n.size++;
								
									n.members[n.size-1]=c.tc[l];

								}
							}
						}
					}

					int alreadyin = 0;
					
					if (n.size==0)
						alreadyin = 1;
					for (int k=0; k<c.lnracompcard; k++) {

						pointset tmp = c.lnracomp[k];
						int common = 0;
						for (int l = 0; l<n.size; l++){
							for (int m =0; m<tmp.size;m++){				
								if (tmp.members[m]==n.members[l]){		
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
				
					if (alreadyin == 0) {								

										
						c.lracompcard++;
			
						c.lracomp[c.lracompcard-1]=rprime;

						sizeofnext++;

						nextLevel[sizeofnext-1].size=rprime.size;
						nextLevel[sizeofnext-1].members=(int*)malloc(rprime.size*sizeof(int));
						for (int y=0; y<rprime.size; y++)
							nextLevel[sizeofnext-1].members[y]=rprime.members[y];
	
						c.lnracompcard++;

						c.lnracomp[c.lnracompcard-1]=n;


						c.assoccomp[c.lracompcard*2-2]=rprime;
						c.assoccomp[c.lracompcard*2-1]=n;

					}
					

				}


			}
			
		}
		lastLevel=nextLevel;											
		nextLevel=NULL;			
		nextLevel=(pointset*)malloc(c.nrep*c.nrep*c.nrep*c.nrep*sizeof(pointset));
		sizeoflast=sizeofnext;

		sizeofnext=0;
	}

	return c;
}

cutdata thirdpreprocess (cutdata c, graph g){

	c.m=(pointset*)malloc(c.lracard*c.nrep*sizeof(pointset));

	for (int i=0;i<c.nrep;i++){
		int v = c.tc[i];
		for (int j=0;j<c.lracard;j++){
			pointset r = c.lra[j];

			pointset rprime;
			rprime.size=r.size;
			rprime.members=(int*)malloc((r.size+1)*sizeof(int));
			for (int k=0;k<r.size;k++)
				rprime.members[k]=r.members[k];
			int alreadyin = 0;
			for (int k=0;k<r.size;k++){
				if (r.members[k]==v)
					alreadyin=1;
			}

			if (alreadyin==0){
				rprime.size++;
				rprime.members[r.size]=v;
			}

			pointset n;
			n.size=0;
			n.members = (int*)malloc(c.nrepincomp*sizeof(pointset));

			for (int k=0;k<c.nrepincomp;k++){
				int neighboor=0;
				for (int l=0;l<rprime.size;l++){
					if (g.matrix[rprime.members[l]*g.size+c.complementtc[k]]==1){
						neighboor=1;
						break;
					}
				}

				if (neighboor==1){
					n.size++;
					n.members[n.size-1]=c.complementtc[k];
				}
			}

			for (int k=0; k<c.lnracard;k++){
				int common = 0;
				for (int l=0; l<n.size;l++){
					for (int m=0; m<c.lnra[k].size; m++){
						if (n.members[l]==c.lnra[k].members[m])
							common++;
					}
				}

				if (common==n.size){

					c.m[j*c.nrep+i]=c.assoc[2*k];
					break;
				}
				else {
					pointset s;
					s.size=0;
					c.m[j*c.nrep+i]=s;
				}
			}
		}
	}


	c.mcomp=(pointset*)malloc(c.lracompcard*c.nrepincomp*sizeof(pointset));

	for (int i=0;i<c.nrepincomp;i++){
		int v = c.complementtc[i];
		for (int j=0;j<c.lracompcard;j++){
			pointset r = c.lracomp[j];

			pointset rprime;
			rprime.size=r.size;
			rprime.members=(int*)malloc((r.size+1)*sizeof(int));
			for (int k=0;k<r.size;k++)
				rprime.members[k]=r.members[k];
			int alreadyin = 0;
			for (int k=0;k<r.size;k++){
				if (r.members[k]==v)
					alreadyin=1;
			}

			if (alreadyin==0){
				rprime.size++;
				rprime.members[r.size]=v;
			}

			pointset n;
			n.size=0;
			n.members = (int*)malloc(c.nrep*sizeof(pointset));

			for (int k=0;k<c.nrep;k++){
				int neighboor=0;
				for (int l=0;l<rprime.size;l++){
					if (g.matrix[rprime.members[l]+g.size*c.tc[k]]==1){
						neighboor=1;
						break;
					}
				}

				if (neighboor==1){
					n.size++;
					n.members[n.size-1]=c.tc[k];
				}
			}

			for (int k=0; k<c.lnracompcard;k++){
				int common = 0;
				for (int l=0; l<n.size;l++){
					for (int m=0; m<c.lnracomp[k].size; m++){
						if (n.members[l]==c.lnracomp[k].members[m])
							common++;
					}
				}

				if (common==n.size){

					c.mcomp[j*c.nrepincomp+i]=c.assoccomp[2*k];
					break;
				}
				else {
					pointset s;
					s.size=0;
					c.mcomp[j*c.nrepincomp+i]=s;
				}
			}
		}
	}


	return c;
}

int toplevelalgorithm (dectree t, graph g){

	if ((t.right==NULL)||(t.left==NULL))
		return EXIT_FAILURE;

	cutdata c1 = cutThatTree (g, t, 0);

	c1 = firstpreprocess (g,c1);

	c1 = secondpreprocess (c1, g);



	c1 = thirdpreprocess (c1, g);


	c1.tab = (int*)malloc(c1.lracard*c1.lnracard*sizeof(int));
	for (int i= 0; i<c1.lracard*c1.lnracard; i++)
		c1.tab[i]=-1;


	cutdata c2 = cutThatTree (g, t, 1);

	c2 = firstpreprocess (g,c2);

	c2 = secondpreprocess (c2, g);

	c2 = thirdpreprocess (c2, g);

	c2.tab = (int*)malloc(c2.lracard*c2.lnracard*sizeof(int));
	for (int i= 0; i<c2.lracard*c2.lnracard; i++)
		c2.tab[i]=-1;

	cutdata *p  = stepalgorithm (t, *(t.left), g);
	cutdata *q = stepalgorithm (t, *(t.right), g);

	if (p!=NULL){
		cutdata c11=p[0];
		cutdata c12=p[1];

		for (int i=0; i< c11.lracard; i++){
			for (int j= 0;j < c12.lracard; j++){
				for (int k=0; k< c1.lnracard; k++){
					pointset ra = c11.lra[i];
					pointset rb = c12.lra[j];
					pointset rwc = c1.lnra[k];

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
					rac.members= (int*)malloc(c11.lnracard*sizeof(int));

					for (int l=0; l<ua.size; l++){
						int z=0;
						int y=0;
						int x=0;
						for (int m=0; m<c11.nacomp;m++){
							if (c11.pointtorepincomp[2*m]==ua.members[l]){
								z=c11.pointtorepincomp[2*m+1];
								break;
							}
						} 
						for (int m=0; m<c11.nrepincomp; m++){
							if (c11.complementtc[m]==z){
								x=m;
								break;
							}
						}
						for (int m=0;m<c11.lnracard;m++){
							int common=0;
							for (int n=0; n<rac.size; n++){
								for (int o=0; o<c11.lnra[m].size; o++){
									if (c11.lnra[m].members[o]==rac.members[n]){
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

						rac = c11.mcomp[y*c11.nrepincomp+x];
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
					rbc.members= (int*)malloc(c12.lnracard*sizeof(int));

					for (int l=0; l<ub.size; l++){
						int z=0;
						int y=0;
						int x=0;
						for (int m=0; m<c12.nacomp;m++){
							if (c12.pointtorepincomp[2*m]==ub.members[l]){
								z=c12.pointtorepincomp[2*m+1];
								break;
							}
						} 
						for (int m=0;m<c12.nrepincomp;m++){
							if (c12.complementtc[m]==z){
								x=m;
								break;
							}
						}
						for (int m=0;m<c12.lnracard;m++){
							int common=0;
							for (int n=0; n<rbc.size; n++){
								for (int o=0; o<c12.lnra[m].size; o++){
									if (c12.lnra[m].members[o]==rbc.members[n]){
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
						rbc = c12.mcomp[y*c12.nrepincomp+x];
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
					rw.members= (int*)malloc(c1.lracard*sizeof(int));

					for (int l=0; l<uw.size; l++){
						int z=0;
						int y=0;
						int x=0;
						for (int m=0; m<c1.na;m++){
							if (c1.pointtorep[2*m]==uw.members[l]){
								z=c1.pointtorep[2*m+1];
								break;
							}
						} 
						for (int m=0;m<c1.nrep;m++){
							if (c1.tc[m]==z){
								x=m;
								break;
							}
						}
						for (int m=0;m<c1.lracard;m++){
							int common=0;
							for (int n=0; n<rw.size; n++){
								for (int o=0; o<c1.lra[m].size; o++){
									if (c1.lra[m].members[o]==rw.members[n]){
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
						rw = c1.m[y*c1.nrep+x];
					}

					int ain=0;
					int bin=0;
					int win=0;
					int acin=0;
					int bcin=0;
					int wcin=0;
						
					for (int l=0;l<c11.lracard;l++){
						int common=0;
						for (int m=0;m<c11.lra[l].size;m++){
							for (int n=0;n<ra.size; n++){
								if (c11.lra[l].members[m]==ra.members[n]){
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

					for (int l=0;l<c12.lracard;l++){
						int common=0;
						for (int m=0;m<c12.lra[l].size;m++){
							for (int n=0;n<rb.size; n++){
								if (c12.lra[l].members[m]==rb.members[n]){
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


					for (int l=0;l<c1.lracard;l++){
						int common=0;
						for (int m=0;m<c1.lra[l].size;m++){
							for (int n=0;n<rw.size; n++){
								if (c1.lra[l].members[m]==rw.members[n]){
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


					for (int l=0;l<c11.lnracard;l++){
						int common=0;
						for (int m=0;m<c11.lnra[l].size;m++){
							for (int n=0;n<rac.size; n++){
								if (c11.lnra[l].members[m]==rac.members[n]){
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


					for (int l=0;l<c12.lnracard;l++){
						int common=0;
						for (int m=0;m<c12.lnra[l].size;m++){
							for (int n=0;n<rbc.size; n++){
								if (c12.lnra[l].members[m]==rbc.members[n]){
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


					for (int l=0;l<c1.lnracard;l++){
						int common=0;
						for (int m=0;m<c1.lnra[l].size;m++){
							for (int n=0;n<rwc.size; n++){
								if (c1.lnra[l].members[m]==rwc.members[n]){
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
					
					if ((c11.tab[ain*c11.lnracard+acin]!=-1)&&(c12.tab[bin*c12.lnracard+bcin]!=-1)){
						if ((c1.tab[win*c1.lnracard+wcin]==-1)||(c1.tab[win*c1.lnracard+wcin]>c11.tab[ain*c11.lnracard+acin]+c12.tab[bin*c12.lnracard+bcin]))						
							c1.tab[win*c1.lnracard+wcin]=c11.tab[ain*c11.lnracard+acin]+c12.tab[bin*c12.lnracard+bcin];
					}					

	
				}
		
			}
		
		}
	}
	else {
		c1.tab[0]=-1;
		c1.tab[1]=0;
		c1.tab[2]=1;
		c1.tab[3]=1;
	}

if (q!=NULL){
			cutdata c21=q[0];
			cutdata c22=q[1];


			for (int i=0; i< c21.lracard; i++){
				for (int j= 0;j < c22.lracard; j++){
					for (int k=0; k< c2.lnracard; k++){
						pointset ra = c21.lra[i];
						pointset rb = c22.lra[j];
						pointset rwc = c2.lnra[k];

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
						rac.members= (int*)malloc(c21.lnracard*sizeof(int));


						for (int l=0; l<ua.size; l++){
							int z=0;
							int y=0;
							int x =0;
							for (int m=0; m<c21.nacomp;m++){
								if (c21.pointtorepincomp[2*m]==ua.members[l]){
									z=c21.pointtorepincomp[2*m+1];
									break;
								}
							} 
							for (int m=0; m<c21.nrepincomp; m++){
								if (c21.complementtc[m]==z){
									x=m;
									break;
								}
							}
							for (int m=0;m<c21.lnracard;m++){
								int common=0;
								for (int n=0; n<rac.size; n++){
									for (int o=0; o<c21.lnra[m].size; o++){
										if (c21.lnra[m].members[o]==rac.members[n]){
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
							rac = c21.mcomp[y*c21.nrepincomp+x];
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
						rbc.members= (int*)malloc(c22.lnracard*sizeof(int));


						for (int l=0; l<ub.size; l++){

							int z=0;
							int y=0;
							int x=0;
							for (int m=0; m<c22.nacomp;m++){
								if (c22.pointtorepincomp[2*m]==ub.members[l]){
									z=c22.pointtorepincomp[2*m+1];
									break;
								}
							} 
							for (int m=0; m<c22.nrepincomp; m++){
								if (c22.complementtc[m]==z){
									x=m;
									break;
								}
							}

							for (int m=0;m<c22.lnracard;m++){
								int common=0;
								for (int n=0; n<rbc.size; n++){

									for (int o=0; o<c22.lnra[m].size; o++){

										if (c22.lnra[m].members[o]==rbc.members[n]){
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

							rbc = c22.mcomp[y*c22.nrepincomp+x];
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
						rw.members= (int*)malloc(c2.lracard*sizeof(int));


						for (int l=0; l<uw.size; l++){
							int z=0;
							int y=0;
							int x=0;
							for (int m=0; m<c2.na;m++){
								if (c2.pointtorep[2*m]==uw.members[l]){
									z=c2.pointtorep[2*m+1];
									break;
								}
							} 
							for (int m=0; m<c2.nrep; m++){
								if (c2.tc[m]==z){
									x=m;
									break;
								}
							}
							for (int m=0;m<c2.lracard;m++){
								int common=0;
								for (int n=0; n<rw.size; n++){
									for (int o=0; o<c2.lra[m].size; o++){
										if (c2.lra[m].members[o]==rw.members[n]){
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
							rw = c2.m[y*c2.nrep+x];
						}
						
						int ain=0;
						int bin=0;
						int win=0;
						int acin=0;
						int bcin=0;
						int wcin=0;
						
						for (int l=0;l<c21.lracard;l++){
							int common=0;
							for (int m=0;m<c21.lra[l].size;m++){
								for (int n=0;n<ra.size; n++){
									if (c21.lra[l].members[m]==ra.members[n]){
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

						for (int l=0;l<c22.lracard;l++){
							int common=0;
							for (int m=0;m<c22.lra[l].size;m++){
								for (int n=0;n<rb.size; n++){
									if (c22.lra[l].members[m]==rb.members[n]){
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


						for (int l=0;l<c2.lracard;l++){
							int common=0;
							for (int m=0;m<c2.lra[l].size;m++){
								for (int n=0;n<rw.size; n++){
									if (c2.lra[l].members[m]==rw.members[n]){
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


						for (int l=0;l<c21.lnracard;l++){
							int common=0;
							for (int m=0;m<c21.lnra[l].size;m++){
								for (int n=0;n<rac.size; n++){
									if (c21.lnra[l].members[m]==rac.members[n]){
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


						for (int l=0;l<c22.lnracard;l++){
							int common=0;
							for (int m=0;m<c22.lnra[l].size;m++){
								for (int n=0;n<rbc.size; n++){
									if (c22.lnra[l].members[m]==rbc.members[n]){
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


						for (int l=0;l<c2.lnracard;l++){
							int common=0;
							for (int m=0;m<c2.lnra[l].size;m++){
								for (int n=0;n<rwc.size; n++){
									if (c2.lnra[l].members[m]==rwc.members[n]){
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

						if ((c21.tab[ain*c21.lnracard+acin]!=-1)&&(c22.tab[bin*c22.lnracard+bcin]!=-1)){
							if ((c2.tab[win*c2.lnracard+wcin]==-1)||(c2.tab[win*c2.lnracard+wcin]>c21.tab[ain*c21.lnracard+acin]+c22.tab[bin*c22.lnracard+bcin]))						
								c2.tab[win*c2.lnracard+wcin]=c21.tab[ain*c21.lnracard+acin]+c22.tab[bin*c22.lnracard+bcin];
						}					
					}
				}
			
			}
		}
		else {
			c2.tab[0]=-1;
			c2.tab[1]=0;
			c2.tab[2]=1;
			c2.tab[3]=1;

		}

	int size=-1;



	for (int i=0;i<c1.lracard;i++){
		for (int j=0;j<c1.lnracard;j++){
			if ((c1.tab[i*c1.lnracard+j]!=-1)&&(c2.tab[j*c1.lracard+i]!=-1)){
				if ((size==-1)||(size > c1.tab[i*c1.lnracard+j]+c2.tab[j*c1.lracard+i]))
					size=c1.tab[i*c1.lnracard+j]+c2.tab[j*c1.lracard+i];
			}
			
		}
	}
	
	return size;
}

cutdata *stepalgorithm (dectree t, dectree tparticular, graph g){
	cutdata *c;

	if ((tparticular.right==NULL)||(tparticular.left==NULL))
		c=NULL;

	else {
		cutdata c1 = cutThatTree (g, tparticular, 0);

		c1 = firstpreprocess (g,c1);


		c1 = secondpreprocess (c1, g);

		c1 = thirdpreprocess (c1, g);

		c1.tab = (int*)malloc(c1.lracard*c1.lnracard*sizeof(int));
		for (int i= 0; i<c1.lracard*c1.lnracard; i++)
			c1.tab[i]=-1;


		cutdata c2 = cutThatTree (g, tparticular, 1);
		c2 = firstpreprocess (g,c2);
		c2 = secondpreprocess (c2, g);

		c2 = thirdpreprocess (c2, g);

		c2.tab = (int*)malloc(c2.lracard*c2.lnracard*sizeof(int));
		for (int i= 0; i<c2.lracard*c2.lnracard; i++)
			c2.tab[i]=-1;

		cutdata *p  = stepalgorithm (t, *(tparticular.left), g);
		cutdata *q = stepalgorithm (t, *(tparticular.right), g);

		if (p!=NULL){
			cutdata c11=p[0];
			cutdata c12=p[1];

			for (int i=0; i< c11.lracard; i++){
				for (int j= 0;j < c12.lracard; j++){
					for (int k=0; k< c1.lnracard; k++){
						pointset ra = c11.lra[i];
						pointset rb = c12.lra[j];
						pointset rwc = c1.lnra[k];
					
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
						rac.members= (int*)malloc(c11.lnracard*sizeof(int));


						for (int l=0; l<ua.size; l++){
							int z=0;
							int y=0;
							int x=0;
							for (int m=0; m<c11.nacomp;m++){
								if (c11.pointtorepincomp[2*m]==ua.members[l]){
									z=c11.pointtorepincomp[2*m+1];

									break;
								}
							} 
							for (int m=0; m<c11.nrepincomp; m++){
								if (c11.complementtc[m]==z){
									x=m;

									break;
								}
							}
							for (int m=0;m<c11.lnracard;m++){
								int common=0;
								for (int n=0; n<rac.size; n++){
									for (int o=0; o<c11.lnra[m].size; o++){
										if (c11.lnra[m].members[o]==rac.members[n]){
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

							rac = c11.mcomp[y*c11.nrepincomp+x];
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
						rbc.members= (int*)malloc(c12.lnracard*sizeof(int));


						for (int l=0; l<ub.size; l++){
							int z=0;
							int y=0;
							int x=0;
							for (int m=0; m<c12.nacomp;m++){
								if (c12.pointtorepincomp[2*m]==ub.members[l]){
									z=c12.pointtorepincomp[2*m+1];
									break;
								}
							} 
							for (int m=0;m<c12.nrepincomp;m++){
								if (c12.complementtc[m]==z){
									x=m;
									break;
								}
							}
							for (int m=0;m<c12.lnracard;m++){
								int common=0;
								for (int n=0; n<rbc.size; n++){
									for (int o=0; o<c12.lnra[m].size; o++){
										if (c12.lnra[m].members[o]==rbc.members[n]){
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
							rbc = c12.mcomp[y*c12.nrepincomp+x];
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
						rw.members= (int*)malloc(c1.lracard*sizeof(int));


						for (int l=0; l<uw.size; l++){
							int z=0;
							int y=0;
							int x=0;
							for (int m=0; m<c1.na;m++){
								if (c1.pointtorep[2*m]==uw.members[l]){
									z=c1.pointtorep[2*m+1];
									break;
								}
							} 
							for (int m=0;m<c1.nrep;m++){
								if (c1.tc[m]==z){
									x=m;
									break;
								}
							}
							for (int m=0;m<c1.lracard;m++){
								int common=0;
								for (int n=0; n<rw.size; n++){
									for (int o=0; o<c1.lra[m].size; o++){
										if (c1.lra[m].members[o]==rw.members[n]){
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

							rw = c1.m[y*c1.nrep+x];

						}
						
						int ain=0;
						int bin=0;
						int win=0;
						int acin=0;
						int bcin=0;
						int wcin=0;
						
						for (int l=0;l<c11.lracard;l++){
							int common=0;
							for (int m=0;m<c11.lra[l].size;m++){
								for (int n=0;n<ra.size; n++){
									if (c11.lra[l].members[m]==ra.members[n]){
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

						for (int l=0;l<c12.lracard;l++){
							int common=0;
							for (int m=0;m<c12.lra[l].size;m++){
								for (int n=0;n<rb.size; n++){
									if (c12.lra[l].members[m]==rb.members[n]){
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


						for (int l=0;l<c1.lracard;l++){
							int common=0;
							for (int m=0;m<c1.lra[l].size;m++){
								for (int n=0;n<rw.size; n++){
									if (c1.lra[l].members[m]==rw.members[n]){
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


						for (int l=0;l<c11.lnracard;l++){
							int common=0;
							for (int m=0;m<c11.lnra[l].size;m++){
								for (int n=0;n<rac.size; n++){
									if (c11.lnra[l].members[m]==rac.members[n]){
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


						for (int l=0;l<c12.lnracard;l++){
							int common=0;
							for (int m=0;m<c12.lnra[l].size;m++){
								for (int n=0;n<rbc.size; n++){
									if (c12.lnra[l].members[m]==rbc.members[n]){
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


						for (int l=0;l<c1.lnracard;l++){
							int common=0;
							for (int m=0;m<c1.lnra[l].size;m++){
								for (int n=0;n<rwc.size; n++){
									if (c1.lnra[l].members[m]==rwc.members[n]){
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

						if ((c11.tab[ain*c11.lnracard+acin]!=-1)&&(c12.tab[bin*c12.lnracard+bcin]!=-1)){
							if ((c1.tab[win*c1.lnracard+wcin]==-1)||(c1.tab[win*c1.lnracard+wcin]>c11.tab[ain*c11.lnracard+acin]+c12.tab[bin*c12.lnracard+bcin]))						
								c1.tab[win*c1.lnracard+wcin]=c11.tab[ain*c11.lnracard+acin]+c12.tab[bin*c12.lnracard+bcin];
						}					
					}

				}
			
			}
		}
		else {
			c1.tab[0]=-1;
			c1.tab[1]=0;
			c1.tab[2]=1;
			c1.tab[3]=1;
		}

		
		if (q!=NULL){
			cutdata c21=q[0];
			cutdata c22=q[1];


			for (int i=0; i< c21.lracard; i++){
				for (int j= 0;j < c22.lracard; j++){
					for (int k=0; k< c2.lnracard; k++){
						pointset ra = c21.lra[i];
						pointset rb = c22.lra[j];
						pointset rwc = c2.lnra[k];

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
						rac.members= (int*)malloc(c21.lnracard*sizeof(int));


						for (int l=0; l<ua.size; l++){
							int z=0;
							int y=0;
							int x =0;
							for (int m=0; m<c21.nacomp;m++){
								if (c21.pointtorepincomp[2*m]==ua.members[l]){
									z=c21.pointtorepincomp[2*m+1];
									break;
								}
							} 
							for (int m=0; m<c21.nrepincomp; m++){
								if (c21.complementtc[m]==z){
									x=m;
									break;
								}
							}
							for (int m=0;m<c21.lnracard;m++){
								int common=0;
								for (int n=0; n<rac.size; n++){
									for (int o=0; o<c21.lnra[m].size; o++){
										if (c21.lnra[m].members[o]==rac.members[n]){
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
							rac = c21.mcomp[y*c21.nrepincomp+x];
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
						rbc.members= (int*)malloc(c22.lnracard*sizeof(int));


						for (int l=0; l<ub.size; l++){

							int z=0;
							int y=0;
							int x=0;
							for (int m=0; m<c22.nacomp;m++){
								if (c22.pointtorepincomp[2*m]==ub.members[l]){
									z=c22.pointtorepincomp[2*m+1];
									break;
								}
							} 
							for (int m=0; m<c21.nrepincomp; m++){
								if (c21.complementtc[m]==z){
									x=m;
									break;
								}
							}

							for (int m=0;m<c22.lnracard;m++){
								int common=0;
								for (int n=0; n<rbc.size; n++){

									for (int o=0; o<c22.lnra[m].size; o++){

										if (c22.lnra[m].members[o]==rbc.members[n]){
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

							rbc = c22.mcomp[y*c22.nrepincomp+x];
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
						rw.members= (int*)malloc(c2.lracard*sizeof(int));


						for (int l=0; l<uw.size; l++){
							int z=0;
							int y=0;
							int x=0;
							for (int m=0; m<c2.na;m++){
								if (c2.pointtorep[2*m]==uw.members[l]){
									z=c2.pointtorep[2*m+1];
									break;
								}
							} 
							for (int m=0; m<c2.nrep; m++){
								if (c2.tc[m]==z){
									x=m;
									break;
								}
							}
							for (int m=0;m<c2.lracard;m++){
								int common=0;
								for (int n=0; n<rw.size; n++){
									for (int o=0; o<c2.lra[m].size; o++){
										if (c2.lra[m].members[o]==rw.members[n]){
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
							rw = c2.m[y*c2.nrep+x];
						}
						
						int ain=0;
						int bin=0;
						int win=0;
						int acin=0;
						int bcin=0;
						int wcin=0;
						
						for (int l=0;l<c21.lracard;l++){
							int common=0;
							for (int m=0;m<c21.lra[l].size;m++){
								for (int n=0;n<ra.size; n++){
									if (c21.lra[l].members[m]==ra.members[n]){
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

						for (int l=0;l<c22.lracard;l++){
							int common=0;
							for (int m=0;m<c22.lra[l].size;m++){
								for (int n=0;n<rb.size; n++){
									if (c22.lra[l].members[m]==rb.members[n]){
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


						for (int l=0;l<c2.lracard;l++){
							int common=0;
							for (int m=0;m<c2.lra[l].size;m++){
								for (int n=0;n<rw.size; n++){
									if (c2.lra[l].members[m]==rw.members[n]){
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


						for (int l=0;l<c21.lnracard;l++){
							int common=0;
							for (int m=0;m<c21.lnra[l].size;m++){
								for (int n=0;n<rac.size; n++){
									if (c21.lnra[l].members[m]==rac.members[n]){
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


						for (int l=0;l<c22.lnracard;l++){
							int common=0;
							for (int m=0;m<c22.lnra[l].size;m++){
								for (int n=0;n<rbc.size; n++){
									if (c22.lnra[l].members[m]==rbc.members[n]){
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


						for (int l=0;l<c2.lnracard;l++){
							int common=0;
							for (int m=0;m<c2.lnra[l].size;m++){
								for (int n=0;n<rwc.size; n++){
									if (c2.lnra[l].members[m]==rwc.members[n]){
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

						if ((c21.tab[ain*c21.lnracard+acin]!=-1)&&(c22.tab[bin*c22.lnracard+bcin]!=-1)){
							if ((c2.tab[win*c2.lnracard+wcin]==-1)||(c2.tab[win*c2.lnracard+wcin]>c21.tab[ain*c21.lnracard+acin]+c22.tab[bin*c22.lnracard+bcin]))						
								c2.tab[win*c2.lnracard+wcin]=c21.tab[ain*c21.lnracard+acin]+c22.tab[bin*c22.lnracard+bcin];
						}					
					}
				}
			
			}
		}
		else {
			c2.tab[0]=-1;
			c2.tab[1]=0;
			c2.tab[2]=1;
			c2.tab[3]=1;
		}


		c=(cutdata*)malloc(2*sizeof(cutdata));
		c[0]=c1;
		c[1]=c2;
	}

	return c;
}
