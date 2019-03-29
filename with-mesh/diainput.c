#include <math.h>
#include "common.h"

void diaIn(char *infile, double *parDia, double *parPos, int *np){
	char filename[20];
	strcpy(filename, infile);
	strcat(filename ,".in"); 
	FILE *pDiaFile = fopen(filename, "rt");

	if (pDiaFile == NULL){
		fprintf(stderr, "Can't open the parameter file %s! \n", filename);
		char c = getchar();
		exit(1);
	}
	int num = 0;

    //printf("No of par %d\n",np);
	findRec(pDiaFile, "PARTICLE");
	for(int i=0; i<*np; i++){
		fscanf(pDiaFile, "%lf", &parDia[i]);
		fscanf(pDiaFile, "%lf", &parPos[i*dim+0]);
		fscanf(pDiaFile, "%lf", &parPos[i*dim+1]);
		fscanf(pDiaFile, "%lf", &parPos[i*dim+2]);
		//for(int i=0; i<dim; i++)fscanf(pDiaFile, "%lf", &parPos[i*dim+i]);
		
		printf("D[%d]: %lf, %lf\n", i, parDia[i], parPos[i*dim+0]);
	
	}

	// printf("num of diameters: %d. \n", num);
	// int k = 0;
	// for (int i=0; i<Tnp/num+1; i++){
	// 	for (int j=0; j<num; j++){
	// 		if (k < Tnp){
	// 			// radius
	// 			hRad[k] = diaTemp[j] / 2.0f; 
	// 			k = k + 1;
	// 		}
	// 	}
	// }

	fclose(pDiaFile);
}



