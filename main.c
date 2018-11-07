#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"

//#include "double_math.h"


int main(void){
    run();
}

void test(){
    // struct MatType *mat;
    // allocateMat(mat);

    // mat[2].density = 1000.0;
    // printf("Density %f\n",mat[1].density);

    // int *num;
    // num  = allocateNumArray(num);
    // for(int i=0; i<num_of_mats; i++) printf("%i\n", num[i]);
    // int *nebListIndex = allocateIntArray(num_of_mats); 
    // nebListIndex[1] = 1;
    // //for(int i=0; i<num_of_mats; i++) printf("%f\n", nebList[i]);  

    // int *boundX;
    // boundX = &nebListIndex[1];
    // printf("Before %i\n",*boundX);
    // nebListIndex[1] = 50;
    // printf("After %i\n",*boundX);
}
void run(){
    // Read input file 
    int np;
    readData("infile", &np);
    printf("PAR  %d\n",np);
    nbList = allocateDoubleArray(np*2);
    parIndex = allocateIntArray(np*2);
    parPosX = allocateDoubleArray(np);
    parPosY = allocateDoubleArray(np);
    parPosZ = allocateDoubleArray(np);
    cellSE = allocateIntArray(np*2);
    parDia = allocateDoubleArray(np);
    cellSE = allocateIntArray(np*3);

    // Read particle information
    diaInput("pardia", parDia, parPosX, parPosY, parPosZ, &np);
    
    // Initialize neighbourlist array for the first time
    initialize(nbList, parIndex, cellSE, np, parPosX, parDia);

    //updateParPosition();
    //parIndex = allocateIntArray(np);
 

    // Delete dynamic memeory allocation
    free(parPosX);
    free(parPosY);
    free(parPosZ);
    free(nbList);
    free(parIndex);
    free(cellSE);
    //free(parIndex);
    free(parDia);
    printf("All good!\n");
    // free(nebListIndex);
}
