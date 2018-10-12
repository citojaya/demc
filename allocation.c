#include "common.h"

void allocateMat(struct MatType *&mt){
    printf("allocation.c\n");
    mt = (struct MatType *)malloc(num_of_mats*sizeof(struct MatType));
}

void allocateNumArray(int *&num){
    num = (int*)malloc(5*sizeof(int));
    memset(num,0,5*sizeof(int));
}
