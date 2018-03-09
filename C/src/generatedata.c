#include "../include/generatedata.h"	

#include <time.h>
#include <stdlib.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

graph generategraph (int nbpoint, int xrange, int yrange, int threshold){

	srand (time(NULL));

	graph g;
	g.size = nbpoint;
	g.matrix = (int*)malloc(nbpoint*nbpoint*sizeof(int));
	g.pos = (int*)malloc(nbpoint*3*sizeof(int));

	for (int i = 0; i < nbpoint; i++){
		int x = rand() % xrange;				// we start by giving each node random x and y et storing this date in the position matrix
		int y = rand() % yrange;
		g.pos[3*i] = i;
		g.pos[3*i+1] = x;
		g.pos[3*i+2] = y;		 
	}
	generateEdges(g,threshold);

	return g;
}

// maybe we'd like to save our graphs for future reviews
int storegraph (graph g, char* name){
	int fd = open (name, O_CREAT|O_TRUNC|O_RDWR);
	chmod (name, 777);
	char* str = (char*) malloc (100);
/*	
	sprintf(str, "%d\n\n", g.size);
	write (fd, str, strlen(str));
*/
	for (int i=0; i<g.size; i++){
		sprintf(str, "%d %d\n", g.pos[3*i+1], g.pos[3*i+2]);
		write (fd, str, strlen(str));
	}

/*
	for (int i=0; i< g.size;i++){
		for (int j=0; j<g.size; j++){
			sprintf (str, "%d|", g.matrix[g.size*i+j]);
			write (fd, str, strlen(str));
		}
		sprintf(str,"\n");
		write (fd, str, strlen(str));
	}
*/
	return fd;
}

//what if we could load-up some graph from file?
graph loadgraph (char* path, int threshold) {
	FILE *f;
	char buffer[21];
	char* subtoken;

	graph g;
	g.size = computeSizeFromFile(path);
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
		g.pos[3*i]=i;
		g.pos[3*i+1]=a;
		g.pos[3*i+2]=b;
	}
/*
	fgets(buffer,80,f);
	
	for (int i = 0; i<g.size; i ++) {
		fgets(buffer,2000,f);
			printf("%s\n",buffer);
		for (int j = 0; j<g.size;j++) {

			int a;
			if (j==0)
				a = atoi(subtoken=strtok(buffer,"|"));
			
			else 
				a = atoi(subtoken=strtok(NULL,"|"));
			printf("a = %d, i = %d, j = %d\n",a, i, j);
			g.matrix[g.size*i+j]=a;
		}
	}
*/
	fclose(f);
	generateEdges(g,threshold);

	return g;
}


void generateEdges (graph g, int threshold){
	for (int i = 0; i<g.size; i++) {
		for (int j= 0; j<i; j++) {
			if ((g.pos[3*i+1]-g.pos[3*j+1])*(g.pos[3*i+1]-g.pos[3*j+1])+(g.pos[3*i+2]-g.pos[3*j+2])*(g.pos[3*i+2]-g.pos[3*j+2])<(threshold*threshold)){				// then we add an edge between points that are near enough (i.e their distance is lower
				g.matrix[g.size*i+j]=1;				// than the threshold squared. Adding an edge means putting 1 in the adjency matrix
				g.matrix[g.size*j+i]=1;
			}
			else {
				g.matrix[g.size*i+j]=0;
				g.matrix[g.size*i+j]=0;
			}
		}
		g.matrix[g.size*i+i]=0;						// no node is his own neighboor
	}
}

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
	char* a = strtok(buffer," ");
	int b = atoi(subtoken=strtok(NULL," "));
	if (b==-1){
		char* c = strtok(NULL," ");
		char * d = strtok(NULL," ");

		t.left = lookTree(f, c);
		t.right = lookTree(f,d);
	}

	else {

		t.label=b;
	}

		fclose(f);

	return t;
}

dectree* lookTree (FILE *f, char *node){
	int found=0;
	dectree t;
	printf("Looking for %s\n",node);
	char buffer[21];
	char* subtoken;
	fseek(f,0,SEEK_SET);

	while ((found==0)&&(fgets(buffer,80,f)!=NULL)){	
		
		char* a = strtok(buffer," ");
		//printf("buffer is %d",strcmp(a,node));
		if (strcmp(a,node)==0){
			printf("Found %s\n",node);
			int b = atoi(subtoken=strtok(NULL," "));
			if (b==-1){
				char* c = strtok(NULL," ");
				char * d = strtok(NULL," ");

				t.left = lookTree(f, c);
				printf("gonna try the right\n");
				t.right = lookTree(f,d);
			}
			else 
				t.label=b;
			
			found=1;

		}
	}
	return &t;

}

int computeSizeFromFile(char* path){
	FILE *f;
	char buffer[21];
	
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
