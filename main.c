#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"

//#include "double_math.h"

#define const int arrSize = 5;
int main(void){
    // int *array = malloc(5*sizeof(int));
    // memset(array,1,5*sizeof(int));
    // printf("%i\n",array[2]);

    
    //allocateAll(h_Mat);


    // struct MatType	*h_Mat = (struct MatType *)calloc(5, sizeof(struct MatType));

    // printf("%f",h_Mat->density);
    // free(h_Mat);
    // printf("completed");

    int *num;
    allocateNumArray(num);
    int x = *num[1];
    //printf("%i\n", (int)*num[1]);

}