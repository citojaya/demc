#ifndef COMMON_H
#define COMMON_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/*--------------Common definitions---------------*/
#define NUM_MAT 2    // two types material
#define PI  3.1415926
#define gravity 9.8f

#define arrSize 5
#define dim 3 // 3D problem 

typedef unsigned int uint;

// err detection utility
#define FPRINTF(a) fprintf a

/*-----------------------------------------------*/
/*--------------Global declaration---------------*/ 
/*-----------------------------------------------*/
//extern char* genfile;
//static int np;
double *oldNL, *newNL; //sorted array for recording particle start and end positions
int *parIndexNL; //sorted array for recording particle index 
int *cellStart; //start cell index for particle
int *cellEnd; //end cell index for particle
double *parPos, *parDia;

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
// 	int dim;; //array dimension (Ex. for 2D dim=2, 3D dim=3)
// 	double *center;
// 	double dia; // particle diameter
// 	int *cStart, *cEnd; //start and end cell indices of neighbour list
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
void diaInput(char *diaFile, double *parDia, double *parPos, int *np);

void insertionSort(double *array, int size);
void updateParPosition();


void test();
void run();


// //*--- particle information ---*//
// extern double3 *hPos, *hVel, *hAngVel, *hForce;
// extern double *hRad;
// extern double3 *dPos, *dVel, *dAngVel, *dForce, *dMom;   //*dAngPos, 
// extern double  *dRMom, *dRad;

#endif 
