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
    //***** Test code ***********
    //largestParDia = 1.0; //(mm)
    xmin = 0.0; //(m)
    xmax = 0.250;
    ymin = -0.005;
    ymax = 0.005;
    zmin = -0.050;
    zmax = 0.050;

    xDiv = floor((xmax-xmin)/(largestParDia*conversion*multif3));
    yDiv = floor((ymax-ymin)/(largestParDia*conversion*multif3));
    zDiv = floor((zmax-zmin)/(largestParDia*conversion*multif3));

    domainDx = (xmax-xmin)/xDiv;
    domainDy = (ymax-ymin)/yDiv;
    domainDz = (zmax-zmin)/zDiv;

    xmin -= 2.*domainDx;
    ymin -= 2.*domainDy;
    zmin -= 2.*domainDz;

    xmax += 2.*domainDx;
    ymax += 2.*domainDy;
    zmax += 2.*domainDz;

    xDiv = floor((xmax-xmin)/(largestParDia*conversion*multif3));
    yDiv = floor((ymax-ymin)/(largestParDia*conversion*multif3));
    zDiv = floor((zmax-zmin)/(largestParDia*conversion*multif3));

    domainDx = (xmax-xmin)/xDiv;
    domainDy = (ymax-ymin)/yDiv;
    domainDz = (zmax-zmin)/zDiv;

    //writeLogLine("DOMAIN domainDx, domainDy, domainDz\n");
    printf("Domain dx,dy,dz %lf, %lf, %lf\n",domainDx, domainDy, domainDz);
    printf("xDiv, yDiv, zDiv %d, %d, %d \n",xDiv, yDiv, zDiv);
    //writeLogNum("logfile.log", "domainDx ", domainDx);
    //printf("%d,%d,%d,%lf,%lf,%lf\n",xDiv,yDiv,zDiv,domainDx,domainDy,domainDz,xmin*conversion,ymin*conversion,zmin*conversion);
    demInit();


    //Testing pointer operation
    //double val1 = 1.0;
    //*val1 = 1.0;
    //*val1 = 1.0;
    //double val2 = 2.0;

    // double *xPos[2];

    // xPos[0] = &val1;
    // xPos[1] = xPos[0];

    // printf("VALUE of 0 position %lf", *xPos[0]);

   
    readParticleData("initial.inj");

    //Read particle information
    
    //Set mass
    for(int i=0; i<np; i++){
        demPart[i].mass = (4.0/3.0)*PI*pow((0.5*demPart[i].dia),3.0)*dens;
        //printf("POSITION %lf,%lf,%lf\n",demPart[i].pos[0],demPart[i].pos[1],demPart[i].pos[2]);
    }

    //Setup DEM scaling 
    setReduceUnits(); 
    //printf("PAR POS %lf,%lf, %lf\n", demPart[0].pos[0]/lengthFactor, demPart[0].pos[1]/lengthFactor,demPart[0].pos[2]/lengthFactor);

    //Set contact surface
 
    printf("Timestep %lf\n",timeStep);
    //Assign to bDBox cells
    
    for(int i=0; i<np; i++){
        int iIndex = ceil((demPart[i].pos[0]-xmin)/domainDx); //ceil gives upper value
        int jIndex = ceil((demPart[i].pos[1]-ymin)/domainDy); //but we need lower value
        int kIndex = ceil((demPart[i].pos[2]-zmin)/domainDz); // therefore -1
        int cI = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;  
        insertToBdBox(i, cI);
    }

    for(int i=0; i<np; i++){
        updateNeighbourList(i);
    }

    //Set DEM gravitational force
    //assignGravity();

/*****************************************************************************************
 * NEW CODE FOR NEIGHBOURLIST
 * **************************************************************************************/
    //particle information

    for(int i=0; i<np; i++){
        printf("Particle %d array position %d, %d\n", i,demPart[i].startIndex[0],demPart[i].endIndex[0]);
        printf("No of neighs %d", demPart[i].noOfNeigh);
        for(int j=0; j<demPart[i].noOfNeigh; j++){
            printf(" neigh %d ", demPart[i].neigh[j]);
        }
        printf("\n");
        printf("\n");
    }
    //printArray(xPos, np*2);

    printf("----------\n");
    // double arr[] = {12.5, 11.1, 13.7, 5.6, 6.6}; 
    // int n = sizeof(arr)/sizeof(arr[0]); 
  
    

    // /***** UPDATE PARTICLE POSITION ************/
    // double tempX = demPart[0].pos[0];
    // demPart[0].pos[0] = demPart[2].pos[0];
    // demPart[2].pos[0] = tempX;

    // xPos[demPart[0].startIndex[0]].value = demPart[0].pos[0]-demPart[0].dia*0.5;
    // xPos[demPart[0].endIndex[0]].value = demPart[0].pos[0]+demPart[0].dia*0.5;

    // xPos[demPart[2].startIndex[0]].value = demPart[2].pos[0]-demPart[2].dia*0.5;
    // xPos[demPart[2].endIndex[0]].value = demPart[2].pos[0]+demPart[2].dia*0.5;
    // /*******************************************/

    //insertionSort(xPos, np*2, 0);
    //insertionSort(yPos, np*2, 1);
    //insertionSort(zPos, np*2, 2);


    // demPart[0].pos[0] = 0.128*lengthFactor;
    // xPos[demPart[0].startIndex[0]].value = demPart[0].pos[0]-demPart[0].dia*0.5;
    // xPos[demPart[0].endIndex[0]].value = demPart[0].pos[0]+demPart[0].dia*0.5;
    // insertionSort(xPos, np*2, 0);
 
    // insertionSort(yPos, np*2, 1);
    // insertionSort(zPos, np*2, 2);

    //printArray(xPos, np*2);

    for(int i=0; i<np; i++){
        printf("Particle %d array position %d, %d\n", i,demPart[i].startIndex[0],demPart[i].endIndex[0]);
    
        printf("No of neighs %d", demPart[i].noOfNeigh);
        for(int j=0; j<demPart[i].noOfNeigh; j++){
            printf(" neigh %d ", demPart[i].neigh[j]);
        }
        printf("\n");
    }
    // //printArray(arr, n); 

    // /****** END OF NEW NEIGHBOUR CODE *********************************************/

    //exit(0);



    printf("CUTGAP %lf\n", cutGap/lengthFactor);
    //Update neighbour list
    for (int i=0; i<np; i++){
        updateNeighbourList(i);
    }
    //exit(0);

    // int count = 0;
    for(int i=0; i<600000; i++){
        demLoop();
   
        cycleCount++;
        if(cycleCount > 500){
            clock_t CPU_time_1 = clock();
            printf("%lf ",demTime/timeFactor);
            int t = CPU_time_1-prevCPUTime;
            printf(" CPU time : %d \n", t);
            prevCPUTime = CPU_time_1;
            demSave();
            cycleCount = 0;
        }
 
    }

    

    /******* Project vector testing code ********/
   /* 
    double *u = allocateDoubleArray(3);
    double *p = allocateDoubleArray(3);

    double *a = allocateDoubleArray(3);
    double *b = allocateDoubleArray(3);
    double *c = allocateDoubleArray(3);
    double *res = allocateDoubleArray(3);
    double *n1 = allocateDoubleArray(3);
    double *n2 = allocateDoubleArray(3);
    double *n3 = allocateDoubleArray(3);

    a[0] = 5.0;
    a[1] = 0.0;
    a[2] = 0.0;

    b[0] = 0.0;
    b[1] = 5.0;
    b[2] = 0.0;

    p[0] = 5.0;
    p[1] = 5.0;
    p[2] = 0.4;

    n1[0] = 0.0;
    n1[1] = 0.0;
    n1[2] = 0.0;

    n2[0] = 1.0;
    n2[1] = 0.0;
    n2[2] = 0.0;

    n3[0] = 0.0;
    n3[1] = 2.0;
    n3[2] = 0.0;

    

    crossProd(a,b,u);
    printf("cross prod %lf,%lf,%lf \n",u[0],u[1],u[2]);
    unitVec(u,u);
    //printf("unit vector %lf,%lf,%lf \n",u[0],u[1],u[2]);

    printf("OVERLAP %lf\n",getOverlap(p,1.0,n3,n2,n1,u));

    //projVec(p,u,res,0);
    //printf("proj vector %lf,%lf,%lf \n",res[0],res[1],res[2]);

    free(res);
    free(u);
    free(p);
    free(a);
    free(b);
    free(c);
    free(n1);
    free(n2);
    free(n3);
    

    for(int i=0; i<5; i++){
        demLoop();
    }
    demSave();
    */


 /******* Testing code for neighbourlist, currently not in use *******/
 
/*********************************************************************************/
    
    // Delete dynamic memeory

  // Delete dynamic memeory
  free(uVec);
  free(ipRVec);
  free(jpRVec);

  free(bdBox);
    
  free(ijVec);
    //free(tempVec);
  free(rotVel);
  free(ipCntPntVel);
  free(jpCntPntVel);
  free(cntPntVel);

  for(int i=0; i<np; i++){
    free(demPart[i].pos);
    free(demPart[i].angVel);
    free(demPart[i].vel);
    free(demPart[i].hisDisp);
    free(demPart[i].force);
    free(demPart[i].momentum);
    free(demPart[i].faceNode1);
    free(demPart[i].faceNode2);
    free(demPart[i].faceNode3);
    free(demPart[i].surfNorm);
  }
  free(demPart);
    
    printf("All good!\n");
    // free(nebListIndex);
}
