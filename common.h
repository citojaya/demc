#ifndef COMMON_H
#define COMMON_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*--------------Common definitions---------------*/
#define NUM_MAT 2    // two types material
#define PI  3.1415926f
#define gravity 9.8f //ms^-2
#define CUTGAP 0.1f //mm
#define conversion 1.0e-3f //convert length values to meters


/*------------ Particle information -------------*/
#define largestParDia 5.0f //mm
#define largestParDensity 2800.0f //kgm^-3
#define arrSize 5
#define dim 3 // 3D problem
#define nbSize 10 //size of neighbourlist

typedef unsigned int uint;

// err detection utility
#define FPRINTF(a) fprintf a

/*-----------------------------------------------*/
/*--------------Global declaration---------------*/ 
/*-----------------------------------------------*/
//extern char* genfile;
//static int np;
double *sortedList; //sorted array for recording particle start and end positions
int *sortedParIndex; //sorted array for recording particle index 
int *cellSE; //keeps a record of start and end positions for all particle (start=1, end=2)
int *parNb; //partilce neighbourlist
int *parNoOfNb; //no of neighbours in a particle
double *parPosX, *parPosY, *parPosZ, *parDia, cutGap;
double refLength,refDensity,lengthFactor,volumeFactor,massFactor,timeFactor,
	densityFactor, forceFactor, pressureFactor, StressFactor, energyFactor, momentFactor,
	powerFactor, velocityFactor, accFactor, angVelFactor, angAccFactor, freqFactor, inertiaFactor;

//int *np;
//*--- Boundary Condition ---*// 
struct CylinderBC { 
	double cir;     // the position of axial and radius (X, Y)
	double R;        // radius
	double Tw;       // top position
	double Bw;       // bottom position
	double topv, btmv;            
};

// material properties
struct MatType {
	double density;
	double emod, ymod, pois;
	double dmpn;
	double sfrc, rfrc;
	double yldp;        
};

// // Particle properties
// struct Particle {
// 	double *center;
// 	center = malloc(sizeof(double)*3);
// 	double dia; // particle diameter
// }



//parIndex		- Array index of particles
//parIndexNL 	- Array for storing particle indices for corresponding neighbour list
//parCord 		- cordinates of particles
//parDia		- particle diameter

//static int num_of_mats = 5;
//static int np; //Total number of particles in the system

void allocateMat(struct MatType *mt);
int *allocateIntArray(int size);
double *allocateDoubleArray(int size);
void readData(char *infile, int *np);
void findRec(FILE *inFile, char* strDest);
void diaInput(char *diaFile, double *parDia, double *parPosX, 
					double *parPosY, double *parPosZ, int *np);

void initialize(double *sortedList, int *sortedParIndex, int *cellSE, int np,
    double *pos, double *parDia);
void insertionSort(double *nbArray, int size, int *parIArray, int *cellSEArray, int firstTime);
void assignNeighbours(double *sortedList, int *sortedParIndex, int *cellSE, int size);
void addNeighbour(int ip, int jp);
void deleteNeighbour(int ip, int jp);


void updateParPosition();

double getCenterDist(int ip, int jp);
void setScaleFactors();


void test();
void run();


// //*--- particle information ---*//
// extern double3 *hPos, *hVel, *hAngVel, *hForce;
// extern double *hRad;
// extern double3 *dPos, *dVel, *dAngVel, *dForce, *dMom;   //*dAngPos, 
// extern double  *dRMom, *dRad;

#endif 
