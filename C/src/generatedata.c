#include "../include/generatedata.h"	

#include <time.h>
#include <stdlib.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

// Generating a geometric graph with nbpoint points ranging from 0 to xrange in x and from 0 to yrange in y. Points are connected by an edge if their distance is lower than the threshold

graph generategraph (int nbpoint, int xrange, int yrange, int threshold){

	srand (time(NULL));

	graph g;
	g.size = nbpoint;
	g.matrix = (int*)malloc(nbpoint*nbpoint*sizeof(int));
	g.pos = (int*)malloc(nbpoint*3*sizeof(int));

	for (int i = 0; i < nbpoint; i++){
		int x = rand() % xrange;				// we start by giving each node random x and y and storing this data in the position matrix
		int y = rand() % yrange;
		g.pos[3*i] = i;
		g.pos[3*i+1] = x;
		g.pos[3*i+2] = y;		 
	}
	generateEdges(g,threshold);

	return g;
}



// Once we've generated a graph g, this function allows you to save this graph in a file, the path being name, for future usage
int storegraph (graph g, char* name){
	int fd = open (name, O_CREAT|O_TRUNC|O_RDWR);						// TODO : add a control structure to prevent crashes due to file access error
	chmod (name, 777);
	char* str = (char*) malloc (100);

	for (int i=0; i<g.size; i++){
		sprintf(str, "%d %d\n", g.pos[3*i+1], g.pos[3*i+2]);			// we use a format in which each line is related to a point of our graph, stating first its x value and then its y value
		write (fd, str, strlen(str));
	}

	return fd;
}



// We're not going to construct the graph in static code, or we'll never be able to run tests on large data sets. This function allows to load a file at path and build the graphs containing all points given in the file and using the threshold to compute edges
graph loadgraph (char* path, int threshold) {
	FILE *f;
	char buffer[21];
	char* subtoken;

	graph g;
	g.size = computeSizeFromFile(path);								// we have to call this function to get size in order to allocate memory for our matrices, even if that forces us to parse the whole file twice
	if ((f=fopen(path,"r")) == NULL) {
		fprintf(stderr,"Error fopen");
		perror("fopen");
	}
	g.pos = (int*)malloc(3*g.size*sizeof(int));
	g.matrix = (int*)malloc(g.size*g.size*sizeof(int));



	for (int i=0; i<g.size; i++){
		fgets(buffer,80,f);
		int a = atoi(subtoken=strtok(buffer," "));
		int b = atoi(subtoken=strtok(NULL," "));
		g.pos[3*i]=i;												// our position matrix containing a field to state the label of the points for which the two succeessing fields give coordinates, we use the variable of incrementation to label the point from 0 to nbpoints
		g.pos[3*i+1]=a;
		g.pos[3*i+2]=b;
	}

	fclose(f);
	generateEdges(g,threshold);										// once the file has been parsed and points are in good order, we call upon generateEdges function to compute all edges

	return g;
}



// Whether you're generating a graph or loading one, this function will allow you to compute edges for points of g, with respect to the given threshold
void generateEdges (graph g, int threshold){

	for (int i = 0; i<g.size; i++) {
		for (int j= 0; j<i; j++) {
			if ((g.pos[3*i+1]-g.pos[3*j+1])*(g.pos[3*i+1]-g.pos[3*j+1])+(g.pos[3*i+2]-g.pos[3*j+2])*(g.pos[3*i+2]-g.pos[3*j+2])<=(threshold*threshold)){				// then we add an edge between points that are near enough (i.e their distance squared is lower
				g.matrix[g.size*i+j]=1;																																// than the threshold squared. Adding an edge means putting 1 in the adjacency matrix
				g.matrix[g.size*j+i]=1;
			}
			else {
				g.matrix[g.size*i+j]=0;																																// if distance between the two points squared is greater than the threshold squared, there are no
				g.matrix[g.size*i+j]=0;																																// edges between those two points and we put a 0 in the adjacency matrix
			}														// as we deal with geometric non-oriented non-reflexive graphs, the adjacency matrix is symetric and we can fill two boxes at the same time, using this symetry
		}
		g.matrix[g.size*i+i]=0;										// in non-reflexive graphs, no node is his own neighboor so we put 0 on the entire diagonal 
	}
}


// Once your graph is loaded, it's wise to load up the decomposition tree for this graph, which will be used by the algorithm we implement here.
dectree loadtree (char* path) {
	FILE *f;
	char buffer[21];
	char* subtoken;

	dectree t;

	if ((f=fopen(path,"r")) == NULL) {
		fprintf(stderr,"Error fopen");
		perror("fopen");
	}



	(fgets(buffer,80,f));
	char* a = strtok(buffer," ");								// our format for the tree file give a node of the tree at each line, giving its label and its corresponding point in the graph. We start by getting the label

	int b = atoi(subtoken=strtok(NULL," "));					// then the corresponding point
	if (b==-1){													// if this point happens to be -1, that means that this tree node does not represent any point of the graph but rather a set of points. In this case, this node will have two children which would be partitions
		char* c = strtok(NULL," ");								// of this larger partition. The file will state each of those children and that's what we're gathering now
		char* e = strtok(NULL," ");
		char* d = (char*)malloc(strlen(e)-1);
		strncpy(d,e,strlen(e)-1);								// as the last point is at the end of the line, we've gotten the line terminator character also. So we escape it to keep only the label we're interested in

		dectree q = lookTree(f, c);								// we then call upon lookTree for each son, which is basically the same function as this one but with a specific node for a target, and without the initialization we've done in this one
		dectree r = lookTree(f,d);
		t.label=-1;
		t.left=(dectree*)malloc(sizeof(dectree));
		t.right=(dectree*)malloc(sizeof(dectree));
		*(t.left) = q;
		*(t.right) = r;											// both dectree return by lookTree are registered by pointers as sons of the current node

	}

	else {														// if the label is not -1, then the current node is a leaf, which represents a point of the graph. It has no sons
		t.label=b;
		t.right=NULL;
		t.left=NULL;
	}

		fclose(f);

	return t;
}


// Basically the same as the previous one without the initializations and specific to one node
dectree lookTree (FILE *f, char *node){
	int found=0;
	dectree t;

	char buffer[81];
	char* subtoken;
	fseek(f,0,SEEK_SET);										// we start reading the file again at the beggining to find the line describing the node we're computing here

	while ((found==0)&&(fgets(buffer,80,f)!=NULL)){				// we read lines until we find the node we're looking for or reach the end of the file. TODO: Had an error if we reach the end of the file without finding the node we're looking for
		
		char* a = strtok(buffer," ");
		if (strcmp(a,node)==0){									// if the first word of the line is the label we're looking for, we're at
																// the right line and we then process the data exactly as in the previous
																// function

			int b = atoi(subtoken=strtok(NULL," "));
			if (b==-1){
				char* c = strtok(NULL," ");
				char* e = strtok(NULL," ");
				char* d = (char*)malloc(strlen(e)-1);
				strncpy(d,e,strlen(e)-1);

				dectree q = lookTree(f, c);
				dectree r = lookTree(f,d);
				t.label=-1;
				t.left=(dectree*)malloc(sizeof(dectree));
				t.right=(dectree*)malloc(sizeof(dectree));
				*(t.left) = q;
				*(t.right) = r;

			}

			else {
				t.label=b;
				t.right=NULL;
				t.left=NULL;
			}
			found=1;

		}
	}

	return t;

}


// This function computes the number of points of the graph by counting lines in the data file until it reaches the end of the document
int computeSizeFromFile(char* path){
	FILE *f;
	char buffer[80];

	if ((f=fopen(path,"r"))==NULL) {
		fprintf(stderr,"Error fopen");
		perror("fopen");
	}

	int i = 0;

	while (fgets(buffer,80,f)!=NULL){
		i++;
	}

	fclose(f);
	return i;
}

int storetree (dectree t, FILE* f, char* node){

	char* str = (char*) malloc (100);
	if (t.label!=-1){
		sprintf(str, "%s %d\n", node, t.label);		
		fwrite ( str, strlen(str), 1, f);	
	}
	else {


		char* tmp1 = (char*)malloc(strlen(node)+1);
		char* tmp2 = (char*)malloc(strlen(node)+1);
		for (int i=0;i<strlen(node);i++){
			tmp1[i]=node[i];
			tmp2[i]=node[i];
		}
		tmp1[strlen(node)]='1';
		tmp2[strlen(node)]='2';
		sprintf(str, "%s -1 %s %s\n", node, tmp1, tmp2);

		fwrite (str, strlen(str), 1, f);	
	//printf("YP %s %s\n",tmp1,tmp2);
		storetree(*(t.left), f, tmp1);
		storetree(*(t.right), f, tmp2);
	}

	return EXIT_SUCCESS;
}
