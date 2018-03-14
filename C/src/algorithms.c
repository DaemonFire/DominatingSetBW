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
// TODO : HOW ARE WE GONNA KEEP THE POINTS OF FATHER TREE IN??? :'( :'( :'(
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
			c.matrixrevisited[i*c.na+j]=g.matrix[c.a[i]*g.size+c.acomp[j]];
		}

	}

	return c;
}


/*
this first pre-process takes a graph, its decomposition tree and a cut of this decomposition tree as input and computes its equivalency 
classes, the output being the lists of representants of the equivalency classes of both primary and secondary sub-graphs, their numbers
and two lists which associate each point of a subgraph to the representant of its equivalency class 
*/
int firstpreprocess(graph g, dectree t, cutdata c){

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
	realloc(c.tc,cursor*sizeof(int));						// we reallocate our pointer in order not to waste memory


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
	realloc(t.complementtc,cursor*sizeof(int));


	/*
	 now that we have reduced the graph to equivalency classes we'd like to reduce the matrix to the representants of those equivalency
	classes in order to speed up future computations ever further
	*/
	int *tmp = (int*) malloc(c.nrep*c.nrepincomp*sizeof(int));
	for (int i=0;i<c.nrep;i++){
		for (int j=0;j<c.nrepincomp;j++){
			tmp[i*c.nrepincomp+j]=c.matrixrevisited[c.tc[i]*c.nacomp+c.complementtc[j]];
		}
	}

	c.matrixrevisited=tmp;

	return EXIT_SUCCESS;

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



/* 
This is the second pre-processing function. We give it a tree and a cutdata on which the first preprocess has been carried on 
and it uses informations on representants to compute the list of representant sets, their respectives neighboorhoods and associations
between those two sets of data, to be stored in the cutdata
*/
int secondpreprocess (dectree t, cutdata c){

	c.lra=NULL;
	c.lnra=NULL;
	pointset *nextLevel=NULL;
	pointset s;
	s.size=0;						
	pointset *lastLevel=(int*)malloc(sizeof(pointset));
	lastLevel[0]=s;
	int sizeoflast=1;
	int sizeofnext=0;
	
	while (lastLevel!=NULL) {								// we initialize lastLevel as a set containing only an empty set of points
		for (int i=0; i<sizeoflast;i++){					// and we loop until lastLevel is an empty set itself
			pointset r = lastLevel[i];						// for each set of points in lastLevel
			for (int j=0; j<c.nrep; j++){					// we try adding each representants of equivalency classes of primary subgraph	
				pointset rprime;
				rprime.size=0;
				int alreadyin = 0;

				for (int k=0;k<r.size;k++){
					rprime.size++;
					realloc (rprime.members, rprime.size*sizeof(int));
					rprime.members[k]=r[k];
					if (rprime.members[k]==c.tc[j])
						alreadyin = 1;
				}
				
				if (already==0) {							// if the vertex was already in the pointset, we have nothing to add, so we pass
					rprime.members[r.size+1]=c.tc[j];		// on to the next representant. Else, we add it to a pointset and compute the
					pointset n;								// neighboorhood of this new set of points
					n.size=0;
					for (int k=0;k<rprime.size;k++){

						for (int l=0; l<c.nrepincomp; l++){	// for each point of this new set we get neighboorhoods

							if(c.matrixrevisited[rprime.members[l]*c.nrepincomp + l]==1){
															// for each representant of equivalence classes in the secondary subgraph, we
															// check if they are neighboors of the point computed

								int alreadyin=0;			// we check that this point has not yet been added to the neighboorhood
								for (int m=0; m<n.size;m++){

									if (n.members[m]==c.complementtc[l]){
										alreadyin=1;
										break;
									}

								}

								if (alreadyin==0){			// else, we add it
									realloc	(n.members, (n.size+1)*sizeof(int));
									n.members[n.size]=c.complementtc[l];
									n.size++;
								}
							}
						}
					}

					// we then check if this neighboorhood is already in lnra
					int alreadyin = 0;
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

					if (already == 0) {									// if n isn't in lnra, we going to add it and the rprime associated
						c.lracard++;
						realloc (c.lra, c.lracard*sizeof(pointset));
						c.lra[c.lracard-1]=rprime;

						sizeofnext++;
						realloc(nextLevel, sizeofnext*sizeof(pointset));
						nextLevel[sizeofnext-1]=rprime;

						c.lnracard++;
						realloc(c.lnra, c.lnracard*sizeof(pointset));
						c.lnra[c.lnracard-1]=n;

						realloc(c.assoc, c.lracard*2*sizeof(pointset));
						c.assoc[c.lracard*2-2]=rprime;
						c.assoc[c.lracard*2-1]=n;
					}
					

				}

				
			}
		}
		lastLevel=nextLevel;											// we loop back until no representating set and no neighboorhoods 
		nextLevel=NULL;													// are interesting enough to be added
	}

	return EXIT_SUCCESS;
}
