#include "common.h"

/*
Initialize neighbourlist array for the first time
Initially neighbourlist array is empty
*/
void initialize(double *nbList, int *parIndex, int *cellSE, int np,
    double *pos, double *parDia){
    int j=0;
    for (int i=0; i<np; i++){
        nbList[j] = pos[i] - 0.5*parDia[i];
        nbList[j+1] = pos[i] + 0.5*parDia[i];
        j += 2;
    }
}


void allocateMat(struct MatType *m){
    printf("allocation.c\n");
    //struct MatType *mt = (struct MatType *)malloc(num_of_mats*sizeof(struct MatType));
    //m = mt;    
}

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

