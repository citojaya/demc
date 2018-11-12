#include "common.h"

/*
Sort neighbourlist for each simulation step
*/
void insertionSort(double *nbArray, int size, int *parIArray, int *cellSEArray, int firstTime) {
   double nbValToInsert;
   int pIValToInsert, cSEValToInsert; 
   int holePosition;
   int i, ip, jp;
  
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
            
            //Execute only after first neighboulist sort iteration
            if(!firstTime){
                ip = sortedParIndex[holePosition];
                jp = sortedParIndex[holePosition-1];
                if(cellSEArray[holePosition] == 1 && cellSEArray[holePosition-1] == 2){
                    deleteNeighbour(ip,jp);
                    deleteNeighbour(jp,ip);
                    //printf("removed\n");
                    //remove from neighbourlist
                }
                //check for end position of particle
                else if(cellSEArray[holePosition] == 2 && cellSEArray[holePosition-1] == 1){
                    addNeighbour(ip,jp);
                    addNeighbour(jp,ip);
                    //printf("added\n");
                    //add to neighbourlist
                }
            }
        }
         //printf(" item moved : %lf\n" , nbArray[holePosition]);
        if(holePosition != i) {
            //printf(" item inserted : %lf, at position : %d\n" , nbValToInsert,holePosition);
            // insert the number at hole position 
            nbArray[holePosition] = nbValToInsert;
            parIArray[holePosition] = pIValToInsert;
            cellSEArray[holePosition] = cSEValToInsert;
        }
   }
}

/*
Add a new neighbour to the neighbour list
*/
void addNeighbour(int ip, int jp){
    parNb[ip*nbSize+parNoOfNb[ip]] = jp;
    parNoOfNb[ip]++;
}

/*
Delete neighbour
*/
void deleteNeighbour(int ip, int jp){
    for (int i=0; i<parNoOfNb[ip]; i++){
        if(parNb[ip*nbSize+i] == jp){
            parNb[ip*nbSize+i] = parNb[ip*nbSize+parNoOfNb[ip]-1];
            parNoOfNb[ip]--;
        }
    }

}

/*
At the begining assign neighbourlist for all particles 
*/	
void assignNeighbours(double *sortedList, int *sortedParIndex, int *cellSE, int size){
    int i, holePosition, ip, jp;
    for (i=1; i<size; i++){
        holePosition = i;
        while(holePosition < size){
            if(cellSE[i] == 1 && cellSE[holePosition+1] == 2){
                //printf("CHECK \n");
                ip = sortedParIndex[i];
                jp = sortedParIndex[holePosition];
                if(getCenterDist(ip,jp)-0.5*(parDia[ip]-parDia[jp]) < cutGap){
                    parNb[ip*nbSize+parNoOfNb[ip]] = jp;
                    parNoOfNb[ip]++;
                    parNb[jp*nbSize+parNoOfNb[jp]] = ip;
                    parNoOfNb[jp]++;
                    //printf("CUTGAP  %lf\n ",cutGap);
                };
                    //check for neighbour region 
            }
            holePosition++;
        }
    }
}

