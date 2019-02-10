#include "common.h"

/*
Initially neighbourlist array is empty. Fill neighbourlist. 
*/
void initialize(double *nbList, int *parIndex, int *cellSE, int np,
    double *pos, double *parDia){
    int j=0;
    for (int i=0; i<np; i++){
        sortedList[j] = pos[i] - 0.5*parDia[i];
        sortedList[j+1] = pos[i] + 0.5*parDia[i];
        parIndex[j] = i; //particle index
        parIndex[j+1] = i; //particle index
        cellSE[j] = 1; //cell start
        cellSE[j+1] = 2; //cell end
        parCIndex[j] = j; //cell start index for each particle
        parCIndex[j+1] = j+1; //cell end index for each particle
        j += 2;
    }
}

/*
Assign scale factors
*/
void setScaleFactors(){
    //Scale factors for reduced units
	refLength = largestParDia*conversion;
	refDensity = largestParDensity;
	lengthFactor = 1.0/refLength;
	volumeFactor = pow(lengthFactor,3);
	massFactor = 6.0/(PI*pow(refLength,3)*refDensity);
	timeFactor = sqrt(gravity/refLength);
	densityFactor = 6.0/(PI*refDensity);
	forceFactor = 6.0/(gravity*PI*pow(refLength,3)*refDensity);
	pressureFactor = 6.0/(gravity*PI*refLength*refDensity);
	StressFactor = pressureFactor;
	energyFactor = 6.0/(gravity*PI*pow(refLength,4)*refDensity);
	momentFactor = energyFactor;
	powerFactor = 6.0/(pow(gravity,1.5)*PI*pow((double)refLength,3.5)*refDensity);
	velocityFactor = 1.0/sqrt(refLength*gravity);
	accFactor = 1.0/gravity;
	angVelFactor = sqrt(refLength/gravity);
	angAccFactor = refLength/gravity;
	freqFactor = sqrt(refLength/gravity);
	inertiaFactor = 6.0/(PI*pow(refLength,5)*refDensity);

    cutGap = CUTGAP*conversion*lengthFactor;
    //printf("refLength  %lf\n ",refLength);
}

/*
void allocateMat(struct MatType *m){
    printf("allocation.c\n");
    //struct MatType *mt = (struct MatType *)malloc(num_of_mats*sizeof(struct MatType));
    //m = mt;    
}*/

/*
Allocate an integer type array
return: int* 
*/
int *allocateIntArray(int size){
    int *val = (int*)malloc(size*sizeof(int));
    memset(val,0,size*sizeof(int));
    return val;
}

/*
Allocate a double* type array
return: double* 
*/
double *allocateDoubleArray(int size){
    double *val = (double*)malloc(size*sizeof(double));
    memset(val,0.0,size*sizeof(double));
    return val;
}

/*
Allocate a char* type array
return: char* 
*/
char *allocateCharArray(int size){
    char *val = (char*)malloc(size*sizeof(char));
    return val;
}

