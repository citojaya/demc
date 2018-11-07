#include "common.h"

void insertionSort(double *nbArray, int size, int *parIArray, int *cellSEArray) {

   double nbValueToInsert;
   int pIValue, cSEValue; 
   int holePosition;
   int i;
  
   // loop through all numbers 
   for(i = 1; i < size; i++) { 
	
      // select a value to be inserted. 
      nbValueToInsert = nbArray[i];
      pIValue = parIArray[i];
      // select the hole position where number is to be inserted 
      holePosition = i;
		
      // check if previous no. is larger than value to be inserted 
      while (holePosition > 0 && nbArray[holePosition-1] > nbValueToInsert) {
         nbArray[holePosition] = nbArray[holePosition-1];
         parIArray[holePosition] = parIArray[holePosition-1];
         cellSEArray[holePosition] = cellSEArray[holePosition-1];
         holePosition--;
         printf(" item moved : %lf\n" , nbArray[holePosition]);
      }

      if(holePosition != i) {
         printf(" item inserted : %lf, at position : %d\n" , nbValueToInsert,holePosition);
         // insert the number at hole position 
         nbArray[holePosition] = nbValueToInsert;
         parIArray[holePosition] = pIValue;
         cellSEArray[holePosition] = cSEValue;
      }

      printf("Iteration %d#:",i);		
   }  
}

