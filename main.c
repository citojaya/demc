#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"

//#include "double_math.h"

#define const int arrSize = 5;
int main(void){
    num_of_mats = 3;


    struct MatType *mat;
    allocateMat(mat);

    mat[2].density = 1000.0;
    printf("Density %f\n",mat[1].density);

    int *num;
    allocateNumArray(num);
    //int x = num[1];
    for(int i=0; i<5; i++) printf("%i\n", num[i]);
    
    // Delete dynamic memeory allocation
    free (num);
    free (mat);

}
