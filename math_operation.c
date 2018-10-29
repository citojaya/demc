#include "common.h"

void insertionSort(double *array, int size) {

   int valueToInsert;
   int holePosition;
   int i;
  
   // loop through all numbers 
   for(i = 1; i < size; i++) { 
	
      // select a value to be inserted. 
      valueToInsert = array[i];
		
      // select the hole position where number is to be inserted 
      holePosition = i;
		
      // check if previous no. is larger than value to be inserted 
      while (holePosition > 0 && array[holePosition-1] > valueToInsert) {
         array[holePosition] = array[holePosition-1];
         holePosition--;
         printf(" item moved : %lf\n" , array[holePosition]);
      }

      if(holePosition != i) {
         printf(" item inserted : %d, at position : %d\n" , valueToInsert,holePosition);
         // insert the number at hole position 
         array[holePosition] = valueToInsert;
      }

      printf("Iteration %d#:",i);		
   }  
}

