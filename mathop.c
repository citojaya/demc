#include "common.h"

void initialSort(double *nbArray, int size, int *parIArray, int *cellSEArray){
   // loop through all numbers 
   for(int i = 1; i < size; i++) { 
   }   
}


void insertionSort(double *nbArray, int size, int *parIArray, int *cellSEArray) {

   double nbValToInsert;
   int pIValToInsert, cSEValToInsert; 
   int holePosition;
   int i;
  
   // loop through all numbers 
   for(i = 1; i < size; i++) { 
        // select a value to be inserted. 
        nbValToInsert = nbArray[i];
        pIValToInsert = parIArray[i];
        cSEValToInsert = cellSEArray[i];
        // select the hole position where number is to be inserted 
        holePosition = i;
		
      // check if previous no. is larger than value to be inserted 
    while (holePosition > 0 && nbArray[holePosition-1] > nbValToInsert) {
        nbArray[holePosition] = nbArray[holePosition-1];
        parIArray[holePosition] = parIArray[holePosition-1];
        cellSEArray[holePosition] = cellSEArray[holePosition-1];
        holePosition--;

        /*bmbn => bnbm neighbourlist unchanged
          bmen => enbm remove from neighbourlist
          emen => enem neighbourlist unchanged
          embn => bnem add to neighbourlist */
        //check for begining position of particle
        if(cellSEArray[holePosition] == 1){
            if(cellSEArray[holePosition-1] == 2){
                printf("removed\n");
                //remove from neighbourlist
            }
        }
        //check for end position of particle
        else if(cellSEArray[holePosition] == 2){
            if(cellSEArray[holePosition-1] == 1){
                printf("added\n");
                //add to neighbourlist
            }
        }
         //printf(" item moved : %lf\n" , nbArray[holePosition]);
    }

    if(holePosition != i) {
         //printf(" item inserted : %lf, at position : %d\n" , nbValToInsert,holePosition);
         // insert the number at hole position 
         nbArray[holePosition] = nbValToInsert;
         parIArray[holePosition] = pIValToInsert;
         cellSEArray[holePosition] = cSEValToInsert;
    }
      //printf("Iteration %d#:",i);		
    }  
}

