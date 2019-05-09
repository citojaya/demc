#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"

//#include "double_math.h"

/* Program begining*/
int main(void){
    printf("START\n");
    run();
}

/* Run the program*/
void run(){
    demInit();
    
    readParticleData("initial.inj");

    //Read particle information
    
    // //Set mass
    // for(int i=0; i<np; i++){
    //     demPart[i].mass = (4.0/3.0)*PI*pow((0.5*demPart[i].dia),3.0)*dens;
    // }

    //Setup DEM scaling 
    setReduceUnits(); 

    printf("Timestep %lf\n",timeStep);

    //Assign to bDBox cells

    //addToBdBox();

    for(int i=0; i<np; i++){
        // writeLog3Num("logfile4.log", "xmin,y,ymin,zmin",xmin/lengthFactor,ymin/lengthFactor,zmin/lengthFactor);
        // writeLog3Num("logfile4.log", "Particle position ", 1e3*demPart[i].pos[0]/lengthFactor,
        //1e3*demPart[i].pos[1]/lengthFactor,1e3*demPart[i].pos[2]/lengthFactor);
      int iIndex = ceil((demPart[i].pos[0]-xmin)/domainDx); //ceil gives upper value
      int jIndex = ceil((demPart[i].pos[1]-ymin)/domainDy); //but we need lower value
      int kIndex = ceil((demPart[i].pos[2]-zmin)/domainDz); // therefore -1
    //   writeLog3Num("logfile4.log", "i,j,k ",iIndex,jIndex,kIndex);
      int cI = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv; 
      //printf("%d\n",cI);
      insertToBdBox(i, cI);
    }

    for(int i=0; i<np; i++){
      updateNeighbourList(i);
    }

    // int count = 0;
    for(int i=0; i<50000; i++){
        //demLoop();
        for(int i=0; i<np; i++){
            forceCalculation(i);
        }
        for(int i=0; i<np; i++){
            updatePosition(i);
        } 
        for(int i=0; i<np; i++){
            updateNeighbourList(i);
        }   
        
        if(cycleCount%500 == 0){
            clock_t CPU_time_1 = clock();
            printf("%lf ",demTime/timeFactor);
            printf(" CPU time : %lu \n", CPU_time_1-prevCPUTime);
            prevCPUTime = CPU_time_1;
            demSave();
            //cycleCount = 0;
        }
        cycleCount++;
        demTime += timeStep;
 
    }

  // Delete dynamic memeory
  free(uVec);
  free(ipRVec);
  free(jpRVec);
  free(bdBox);
  free(ijVec);
  free(rotVel);
  free(ipCntPntVel);
  free(jpCntPntVel);
  free(cntPntVel);

  free(demPart);
    
  printf("All good!\n");
    // free(nebListIndex);
}
