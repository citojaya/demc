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

void readParticleData(char *infile){
	char com[21];
	int ret;
	//int np = 0;
	FILE *injFile = fopen(infile, "r");
	ret = fscanf(injFile,"%d",&np);
	double *x;

	//while(ret!=EOF)
	for(int i=0; i<np; i++)
	{
		ret = fscanf(injFile,"%lf",&demPart[i].pos[0]);
		ret = fscanf(injFile,"%lf",&demPart[i].pos[1]);
		ret = fscanf(injFile,"%lf",&demPart[i].pos[2]);
		ret = fscanf(injFile,"%lf",&demPart[i].dia);
	}
	fclose(injFile);
	printf("No of particles in the simulation %d \n",np);
}

/*---Read input data from a file ----*/
void readInput(char *infile, int *np, double *dens, double *ymod, 
			double *pois, double *sfc, double *rec, double *dmpn, double *rf, 
			double *cyldia, double *dt, int *nW, int *updateDPM, double *haConst){
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

	findRec(InFile, "MATERIAL");
	fscanf(InFile, "%lf", dens);
	fscanf(InFile, "%lf", ymod);
	fscanf(InFile, "%lf", pois);
	fscanf(InFile, "%lf", sfc);
	fscanf(InFile, "%lf", dmpn);
	fscanf(InFile, "%lf", rf);
	fscanf(InFile, "%lf", haConst);

	findRec(InFile, "CylinderBC");
	fscanf(InFile, "%lf", cyldia);
	fprintf(LogFile,"Domain size\n");
	fprintf(LogFile,"xDiv,yDiv,zDiv: %d,%d,%d\n :",xDiv,yDiv,zDiv);

	findRec(InFile, "SIMULATION");
	fscanf(InFile, "%lf", dt);	

	findRec(InFile, "WALLS");
	fscanf(InFile, "%d", nW);	

	findRec(InFile, "DPM");
	fscanf(InFile, "%d", updateDPM);	
	//fprintf(LogFile,"xmin,xmax %lf,%lf\n :",xmin,xmax);
	//fprintf(LogFile,"ymin,ymax %lf,%lf\n :",ymin,ymax);
	//fprintf(LogFile,"zmin,zmax %lf,%lf\n :",zmin,zmax);

    fclose(InFile);
	fclose(LogFile);
}

void writeLogNum(char *infile, char *line, double num){
	FILE *LogFile = fopen(infile, "a");
	//fprintf(LogFile,line);
	fprintf(LogFile,"%lf\n",num);
	fclose(LogFile);
}

void writeLogLine(char *infile, char *line){
	// FILE *LogFile = fopen(infile, "a");
	// fprintf(LogFile,line);
	// fclose(LogFile);
}

// void writeInjectionFile(char *infile){
// 	char filename[20];
// 	strcpy(filename, infile);
// 	FILE *injFile = fopen(filename, "w");
// 	Injection *I;
//     Injection *Ilist = Get_dpm_injections();

//     loop(I,Ilist)
//     {
//     	Particle *p;
//         loop(p,I->p)
//         { 
// 			char line[100];
// 			strcpy(line, "((");
// 			//strcat(filename ,".in"); 
// 			//fprintf(injFile,"((");
// 			char pX[8];
// 			sprintf(pX, "%f", demPart[p->part_id].pos[0]/lengthFactor);
// 			strcat(line ,pX); 
// 			strcat(line ," "); 
// 			char pY[8];
// 			sprintf(pY, "%f", demPart[p->part_id].pos[1]/lengthFactor);
// 			strcat(line ,pY); 
// 			strcat(line ," ");
// 			char pZ[8];
// 			sprintf(pZ, "%f", demPart[p->part_id].pos[2]/lengthFactor);
// 			strcat(line ,pZ);
// 			strcat(line ," "); 
// 			char pVelX[8];
// 			sprintf(pVelX, "%f", demPart[p->part_id].vel[0]/velocityFactor);
// 			strcat(line ,pVelX);
// 			strcat(line ," ");
// 			char pVelY[8];
// 			sprintf(pVelY, "%f", demPart[p->part_id].vel[1]/velocityFactor);
// 			strcat(line ,pVelY);
// 			strcat(line ," ");
// 			char pVelZ[8];
// 			sprintf(pVelZ, "%f", demPart[p->part_id].vel[2]/velocityFactor);
// 			strcat(line ,pVelZ);
// 			strcat(line ," ");
// 			char pDia[8];
// 			sprintf(pDia, "%f", demPart[p->part_id].dia/lengthFactor);
// 			strcat(line ,pDia);
// 			strcat(line ," 0.0 1.0))\n");
// 			fprintf(injFile, line);
// 		}
// 	}

// 	fclose(injFile);
// }

void diaInput(char *infile, struct demParticle *par, int *np){
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

	double pDia, pX, pY, pZ;
    //printf("No of par %d\n",np);
	findRec(pDiaFile, "PARTICLE");
	for(int i=0; i<*np; i++){
		fscanf(pDiaFile, "%lf", &pDia);
		fscanf(pDiaFile, "%lf", &pX);
		fscanf(pDiaFile, "%lf", &pY);
		fscanf(pDiaFile, "%lf", &pZ);
		demPart[i].dia = pDia;
		demPart[i].pos[0] = pX;
		demPart[i].pos[1] = pY;
		demPart[i].pos[2] = pZ;
		//printf("D[%d]: %lf, %lf, %lf, %lf\n", i, demPart[i].dia, demPart[i].posX, demPart[i].posY, demPart[i].posZ);
	}
	fclose(pDiaFile);
}

void readWalls(char *infile, int *walls){
	char filename[20];
	strcpy(filename, infile);
	strcat(filename ,".in"); 
	FILE *wFile = fopen(filename, "rt");

	if (wFile == NULL){
		fprintf(stderr, "Can't open the parameter file %s! \n", filename);
		char c = getchar();
		exit(1);
	}

	findRec(wFile, "WALL_NO");
	for(int i=0; i<noOfWalls; i++){
		int wallNo;
		fscanf(wFile, "%d", &wallNo);
		walls[i] = wallNo;
	}
	fclose(wFile);
}

void demSave(){
	//FILE *outfile; 
	char filename[20];
	sprintf(filename, "particle.dat");
	FILE *outfile = fopen(filename, "a");
	fprintf(outfile, "TIME = %lf\n",demTime);

  	// Injection *I;
  	// Injection *Ilist = Get_dpm_injections();
  
  	// Update FLUENT particle postion and velocity 
  	int ip = 0; 
	 for(int i=0; i<np; i++)
  	 {
     	//Particle *p;
     	// loop(p,I->p)
     	// {
	 		// fprintf(outfile, "%11.5lf   %11.5lf   %11.5lf   %11.5f  %11.5lf   %11.5lf  %11.5lf\n",
	 		// P_POS(p)[0]/conversion, P_POS(p)[1]/conversion,P_POS(p)[2]/conversion,
	 		// P_VEL(p)[0],P_VEL(p)[1],P_VEL(p)[2], P_DIAM(p)/conversion);

		 	fprintf(outfile, "%11.5lf   %11.5lf   %11.5lf   %11.5f  %11.5lf   %11.5lf  %11.5lf\n",
			demPart[i].pos[0]/(lengthFactor*conversion), demPart[i].pos[1]/(lengthFactor*conversion),
			demPart[i].pos[2]/(lengthFactor*conversion),
			demPart[i].vel[0]/(velocityFactor),demPart[i].vel[1]/(velocityFactor),demPart[i].vel[2]/(velocityFactor),
			demPart[i].dia/(lengthFactor*conversion));
	 	//}
			
     }
 
	fclose(outfile);
	
}

// void writeTec(double *pPosX, double *parPosY, double *parPosZ){
// 	FILE *outfile; 
// 	char filename[20];
// 	sprintf(filename, "particle_info.dat");

// 	if (h_dem->Outs == 0)
// 	{
// 		outfile = fopen(filename, "wt");

// 		fprintf(outfile, "TITLE = \" PARTICLE INFORMATION \" \n");
// 		fprintf(outfile, "VARIABLES = X   Y   Z   R   VX   VY   VZ   W   F\n");
// 		fclose(outfile);
// 	}

// 	outfile = fopen(filename, "a");
// 	fprintf(outfile, "ZONE T= \" %12.6lf s \" \n", h_dem->ctime * Rdu->rtunit);
// 	for (int ip = 0; ip<TotalParticle; ip++)
// 	{
// 		fprintf(outfile, "%11.5lf   %11.5lf   %11.5lf   %11.5f  %11.5lf   %11.5lf   %11.5lf   %11.5lf   %11.5lf\n", 
// 		                  hPos[ip].x,   hPos[ip].y,   hPos[ip].z,  hRad[ip],
// 		                  hVel[ip].x,   hVel[ip].y,   hVel[ip].z, 
// 		                  length(hAngVel[ip]), length(hForce[ip]));
// 	}

// 	fclose(outfile);
// }

