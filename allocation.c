#include "common.h"

void allocateAll(struct MatType *mt){
    printf("allocation.c\n");
    mt = (struct MatType *)calloc(5, sizeof(struct MatType));
    //free(mt);
}

void allocateNumArray(int *num){
    num = (int *)malloc(5*sizeof(int));
    // memset(num,0,5*sizeof(int));
    printf("%i\n",num[2]);
}