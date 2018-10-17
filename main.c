#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"

//#include "double_math.h"

#define const int arrSize = 5;

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
}
void run(){
    double *nebList = allocateDoubleArray(num_of_mats);
    int *nebListIndex = allocateIntArray(num_of_mats); 
    for(int i=0; i<num_of_mats; i++) printf("%f\n", nebList[i]);  

    int *boundX;
    boundX = nebListIndex[1];
  
    // Delete dynamic memeory allocation
    //free (num);
    //free (mat);
    free(nebList);
    free(nebListIndex);
}
