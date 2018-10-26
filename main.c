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
    readData("infile");

    oldNL = allocateDoubleArray(np*2);
    newNL = allocateDoubleArray(np*2);
    parIndexNL = allocateIntArray(np*2);
    parPos = allocateDoubleArray(np*3);
    parDia = allocateDoubleArray(np);

    // Read particle information
    diaInput("pardia", parDia, parPos);
    //parIndex = allocateIntArray(np);
 

    // Delete dynamic memeory allocation
    free(oldNL);
    free(newNL);
    free(parIndexNL);
    free(parPos);
    //free(parIndex);
    free(parDia);
    printf("All good!\n");
    // free(nebListIndex);
}
