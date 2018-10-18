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
}
void run(){
    double *nebList = allocateDoubleArray(num_of_mats);
    int *nebListIndex = allocateIntArray(num_of_mats); 
    nebListIndex[1] = 1;
    //for(int i=0; i<num_of_mats; i++) printf("%f\n", nebList[i]);  

    int *boundX;
    double *parCenter = (double *)malloc(noOfPar*dim*sizeof(double));
    memset(parCenter,0.0,noOfPar*dim*sizeof(double));
    int count = 0;
    for (int i=0; i<noOfPar; i++){
        for (int j=0; j<dim; j++){
            *(parCenter+i*dim+j) = ++count;
        }
    }


     printf("partilce center %f",*(parCenter+5+1));

    boundX = &nebListIndex[1];
    printf("Before %i\n",*boundX);
    nebListIndex[1] = 50;
    printf("After %i\n",*boundX);
  
    // Delete dynamic memeory allocation
    //free (num);
    //free (mat);
    free(nebList);
    free(nebListIndex);
    free (parCenter);
}
