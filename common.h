#ifndef COMMON_H
#define COMMON_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "double_math.h"

/*--------------Common definitions---------------*/
#define NUM_MAT 2    // two types material
#define PI  3.1415926
#define gravity 9.8f

#define arrSize 5
#define dim 3
#define noOfPar 5
// #define threadsPerBlock  128
// #define IY  index+d_Params.NP
// #define IZ  index+2*d_Params.NP
typedef unsigned int uint;

// err detection utility
#define FPRINTF(a) fprintf a

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


/*-----------------------------------------------*/
/*--------------Global declaration---------------*/ 
/*-----------------------------------------------*/
//extern char* genfile;

static int num_of_mats = 5;

void allocateMat(struct MatType *mt);
int *allocateIntArray(int size);
double *allocateDoubleArray(int size);
void test();
void run();
// //*--- particle information ---*//
// extern double3 *hPos, *hVel, *hAngVel, *hForce;
// extern double *hRad;
// extern double3 *dPos, *dVel, *dAngVel, *dForce, *dMom;   //*dAngPos, 
// extern double  *dRMom, *dRad;

#endif 
