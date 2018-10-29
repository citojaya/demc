#include "common.h"

/*--- Find the reading position in the file--*/
void findRec(FILE *inFile, char* strDest){
	int nbytes = 256;
	char* strSrc;
	strSrc = (char *)malloc(nbytes+1);

	rewind(inFile);
	int n=strlen(strDest);
	while(!feof(inFile)){
		fgets(strSrc, 256, inFile);
		strSrc[n]='\0';
		if (strcmp(strDest, strSrc) == 0){
			break;
		}
	}

	if(strcmp(strDest, strSrc) != 0){
		//free(strSrc);
		//printf("Unable to find relevant info of: %s \n", strDest);
		exit(1);
	}
	free(strSrc);
}

/*---Read input data from a file ----*/
void readData(char *infile, int *np){
	// input file reading
	char filename[20];
	strcpy(filename, infile);
	strcat(filename ,".in"); 
	FILE *InFile = fopen(filename, "rt");

	if (InFile == NULL){
		fprintf(stderr, "Can't open the parameter file %s! \n", filename);
		char c = getchar();
		exit(1);
	}
	// Log file
	char LogjName[20];
	strcpy(LogjName, infile);
	strcat(LogjName ,"_Log.dat"); 
	FILE *LogFile = fopen(LogjName, "a");

	// expand size for computing area
	double exComp = 0.0;

	// loading phase
	double DieDepth = 1.0;
	double UnDepth  = 0.0;
	double BdDepth  = 0.0;

	double parDia = 0.0;
	
	findRec(InFile, "PAR_NUMBER");
	fscanf(InFile, "%d",  np);
	//printf("np: %d\n",np);
	//findRec(InFile, "packing_TIME");
	//fscanf(InFile, "%lf", &h_load->tPacking);	// packing time

    fclose(InFile);
	fclose(LogFile);
}

