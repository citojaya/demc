
#ifndef COMMON_H
#define COMMON_H

//#include "udf.h"
//#include "dpm.h"
//#include "dem.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

/*--------------Common definitions---------------*/
#define NUM_MAT 2    // two types material
#define PI 3.1415926f
#define gravity 9.8f //ms^-2
//#define CUTGAP 0.5f //mm
#define conversion 1.0e-3f //convert length values to meters
#define overlapLimit 0.005f //ovelap cuttoff value

/*------------ Particle information -------------*/
#define largestParDia 0.130f //mm
#define largestParDensity 883.0f //kgm^-3
#define multif3 2.2 //multification factor used in bounding box divisions

#define DIM 3 // 3D problem
//#define NO_OF_PARTICLES 100 //Particle array size
#define NBSIZE 40 //size of neighbourlist
#define NO_OF_FLUID_CELLS 1 //number of fluid cells in a bounding box
#define NO_OF_FACES 20 //maximum number of faces in a wall mesh
#define NO_OF_EDGES 5 //maximum number of edges in a wal mesh 
#define NO_OF_PARTICLES_IN_BDCELL 30
//Testing parameters
#define FVOLF 0.5f //fluid volume fraction
#define BETA 0.4f //interphase momentum exchange coefficient
#define yMinBottom  0.0f//Y coordinate of bottom surface

/*----------- Fluid properties -------------------*/
#define VISCOSITY 1.7893e-5f //kg/m-s
#define POROSITY 0.5f
//typedef unsigned int uint;

// err detection utility
#define FPRINTF(a) fprintf a


//Domain *fd; //fluid domain
//Thread *t; //fluid cell thread
/*-----------------------------------------------*/
/*--------------Global declaration---------------*/ 
/*-----------------------------------------------*/
//extern char* genfile;
//static int np;
//double *sortedList; //sorted array for recording particle start and end positions
//double *nrmDisp;
FILE *LogFile;
int prevCPUTime;
int time_count;
int particle_counter;
int noOfWallSurfaces, noOfVertices;
//int noOfWalls;
//int *walls;
int np, parArraySize; //number of particles, particle array size 
int cycleCount;
unsigned int initialized;
//int *sortedParIndex; //sorted array for recording particle index 
//int *cellSE; //keeps a record of start and end positions for all particle (start=1, end=2)
//int *parNb; //partilce neighbourlist
//int *parCIndex; //particle start and end positons of sortedList
//int *parNoOfNb; //no of neighbours in a particle
//double *parPosX,// *parPosY, *parPosZ, *parDia, *parInert, *parMass, cutGap; //particle parameters
double cutGap; //Neighbour region
//double largestParDia; //Largest particle diameter
double cellRadius; //Radius defined by bounding box cell given by 0.5*sqrt(dx^2+dy*2)
double refLength,refDensity,lengthFactor,volumeFactor,massFactor,timeFactor,
	densityFactor, forceFactor, pressureFactor, StressFactor, energyFactor, momentFactor,
	powerFactor, velocityFactor, accFactor, angVelFactor, angAccFactor, freqFactor, inertiaFactor;
int xDiv, yDiv, zDiv; //Number of divisions in orthogonal domain cells



// Wall mesh 
struct BdBox *wallMeshBox;
int wmxDiv, wmyDiv, wmzDiv;
double wallMxmin, wallMxmax, wallMymin, wallMymax, wallMzmin, wallMzmax; //wall mesh boundary
double wallMdx, wallMdy, wallMdz;

//DEM mesh
struct BdBox *bdBox;
double domainDx, domainDy, domainDz; //Domain cell size
double xmax, ymax, zmax; //Max values of domain boundary 
double xmin, ymin, zmin; //Min values of domain boundary

/*
rIn = largestParDia
rOut = 1.55*rIn
allowedDisp = (rOut-rIn)/2
*/
double rIn, rOut, allowedDisp; //if particle displacement > allowedDisp -> update neighbourlist

double *uVec, *ipRVec, *jpRVec, *ijVec, *rotVel;
double *ipCntPntVel, *jpCntPntVel, *cntPntVel;
double dens, ymod, pois, sfc, rf, rec, dmpn, elasticMod, haConst; //particle material property
double dsmaxCff, dti, dd, dsmax, fdt; //used in contact force calculation
double cyldia; //cylinder diameter
double timeStep, demTime, maxTime;


struct wallFace *face;
struct line *edge;
struct demParticle *demPart;
struct position *xPos;
struct position *yPos;
struct position *zPos;

//Tracked_Particle *dpmList[2];
int updateDPM;
int saveDEM; //counter used for saving DEM particle for TECPLOT

//Position details of xstart, xend (used in neighbourlist update)
struct position{
	int partId;
	double value;
	int type;
};


//Walls face
struct wallFace{
	double node1[DIM];
	double node2[DIM];
	double node3[DIM];
	double centroid[DIM];
};

struct line{
	double node1[DIM];
	double node2[DIM];
	double center[DIM];
};

//Bounding box which holds CFD cells and DEM particles
struct BdBox{
	double fluidVelX, fluidVelY, fluidVelZ;
	double pGradX, pGradY, pGradZ;
	double fluidVolF;
	double dragFX, dragFY, dragFZ;
	int totalFaces; //number of contact faces
	int noOfParticles; //number of DEM particles
	int parts[NO_OF_PARTICLES_IN_BDCELL];
	//struct wallFace *face;
	int noOfFaces, noOfEdges;
	int wallFaceId[NO_OF_FACES];
	int edgeId[NO_OF_EDGES];
	int isBoundary, isEdgeBoundary;
	//double faceCentroid[NO_OF_FACES*DIM];
};


//Particle
struct demParticle{
	double dt; //particle time step
	double currentTime; //current time
	double displacement; //if displacement > rMax update neighbourlist
	double dia, inert, mass, nrmDisp;
	double *pos, *angVel, *vel, *hisDisp, *force, *momentum;
	double *faceNode1, *faceNode2, *faceNode3, *surfNorm;
	int neigh[NBSIZE]; // array for neighbour list
	int contDone[NBSIZE]; //keeps a reocrd of neignbour particles which has already force calculation done
	int startIndex[DIM]; //start index of position array 
	int endIndex[DIM]; //end index of position array
	//int noOfCntF;

	//int neighCntCalculated[NBSIZE];
	int noOfNeigh;
	int contacts;//Number of neighbour contacts completed
	int prevCellIndex;
	//int insertable;
	double ha; //Hamaker constant
	//int pCntFcalculated;

};

//parIndex		- Array index of particles
//parIndexNL 	- Array for storing particle indices for corresponding neighbour list
//parCord 		- cordinates of particles
//parDia		- particle diameter

//static int num_of_mats = 5;
//static int np; //Total number of particles in the system
double solidFraction(int ip);
void writeLogNum(char *infile, char *line, double num);
//void test(Tracked_Particle *p, Thread *t);
//void allocateMat(struct MatType *mt);
int *allocateIntArray(int size);
double *allocateDoubleArray(int size);

struct BdBox *allocateBdBoxArray(int size);
struct demParticle *allocatePar(int np);
//void insertCellToBBox(cell_t c, double x[]);
double partVol(int p);
void insertToBdBox(int p, int cI);
int isInside(double *n1, double *n2, double *n3, double *pt);
double area(double *vec1, double *vec2);
//void addToBdBox();
void addFacesAndEdges();
void checkFaceContact(int p);
int faceExist(int cI, int fId);
int edgeExist(int cI, int eId);
void readInput(char *infile, int *np, double *dens, double *ymod, 
			double *pois, double *sfc, double *rec, double *dmpn, double *rf, double *cyldia, 
			double *dt, int *updateDPM, double *haConst);
void findRec(FILE *inFile, char* strDest);
void diaInput(char *diaFile, struct demParticle *par, int *np);
//void readWalls(char *infile, int *walls);
void readParticleData(char *infile);
void readSurface(char *infile);
void readEdges(char *infile);
void writeTec();
//void demInit(int *xD, int *yD, int *zD, double xmin,double xmax,
//					double ymin,double ymax,double zmin,double zmax);
void demInit();
void buildDEMCFDCellMap();
void copyDEMInfo();
void demLoop();
void demSave();
void allocate();

//void updateForce(Tracked_Particle *p);
void updateForce(int p);
void partContactForce(int ip, int jp, double nrmDsp);
void boundaryContactForce(int pI, double *n1, double *n2, double *n3, double *uVec);
void ppVWForce(int ip, int jp, double vGap);
void pWallVWForce(int p, double vGap, double *uVec);
//void dragForce(int pI);
void assignGravity();
double getOverlap(double *parPos, double dia, double *n1, double *n2, double *n3, double *uVec);

int contactNotDone(int ip, int jp);
void neighbourContactForce(int pI);
void surfaceContactForce(int p, double nrmDisp, double *uVec);
void vecAdd(double *v1, double *v2, double *vec);
void crossProd(double *v1, double *v2, double *vec);
void vecSub(double *v1, double *v2, double *vec);
void unitVec(double *v, double *vec);
void sclMult(double scl, double *vec);
void sclVecMult(double scl, double *inVec, double *outVec);
double vecMag(double *vec);
void projVec(double *v1, double *n, double *vec, int type);
double dotProduct(double *v1, double *v2);
void getUnitVector(double *v1, double *v2, double *v3, double *uVec);

void initialize(double *sortedList, int *sortedParIndex, int *cellSE, int np,
    double *pos, double *parDia);
//void insertionSort(double *nbArray, int size, int *parIArray, int *cellSEArray, int firstTime);
//void assignNeighbours(double *sortedList, int *sortedParIndex, int *cellSE, int size);
int insertable(int ip, int jp);
void addNeighbour(int  ip, int jp);
void updateNeighbourList(int p);
int insertable (int ip, int jp);
void deleteNeighbour(int ip, int jp);
void update(double *pX, int np);

double getCenterDist(int ip, int jp);
void setReduceUnits();
void updateBBFluidVel();
void insertToBdBox(int p, int cI);
void deleteParticle(int p, int cI);
void forceCalculation(int p);
void updatePosition(int p);
void checkXContact(int p, double xMin, double xMax);
void checkYContact(int p, double yMin, double yMax);
void checkZContact(int p, double zMin, double zMax);
//void insertionSort(double arr[], int n);
void insertionSort(struct position *xp, int n, int type);
void printArray(struct position *xp, int n);
//void printArray(double arr[], int n);
void run();
void deleteAll();

#endif 
