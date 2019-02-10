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
    largestParDia = 0.1; //(mm)
    xmin = 0.042; //(m)
    xmax = 0.057;
    ymin = -0.001;
    ymax = 0.001;
    zmin = -0.007;
    zmax = 0.014;
    // xmin = xmin - 4.0*largestParDia*conversion;
    // xmax = xmax + 4.0*largestParDia*conversion;
    // ymin = ymin - 4.0*largestParDia*conversion;
    // ymax = ymax + 4.0*largestParDia*conversion;
    // zmin = zmin - 4.0*largestParDia*conversion;
    // zmax = zmax + 4.0*largestParDia*conversion;
    xDiv = floor((xmax-xmin)/(largestParDia*conversion*multif3));
    yDiv = floor((ymax-ymin)/(largestParDia*conversion*multif3));
    zDiv = floor((zmax-zmin)/(largestParDia*conversion*multif3));

    domainDx = (xmax-xmin)/xDiv;
    domainDy = (ymax-ymin)/yDiv;
    domainDz = (zmax-zmin)/zDiv;

    xmin -= domainDx;
    ymin -= domainDy;
    zmin -= domainDz;

    xmax += domainDx;
    ymax += domainDy;
    zmax += domainDz;

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
    readParticleData("initial.inj");

    //Read particle information
    
    //Set mass
    for(int i=0; i<np; i++){
        demPart[i].mass = (4.0/3.0)*PI*pow((0.5*demPart[i].dia),3.0)*dens;
        //printf("POSITION %lf,%lf,%lf\n",demPart[i].pos[0],demPart[i].pos[1],demPart[i].pos[2]);
    }

    //Setup DEM scaling 
    setReduceUnits(); 

    //Set contact surface
 
    printf("Timestep %lf\n",timeStep);
    //Assign to bDBox cells
    addToBdBox();

    //Set DEM gravitational force
    //assignGravity();
    
    printf("CUTGAP %lf\n", cutGap/lengthFactor);
    //Update neighbour list
    for (int i=0; i<np; i++){
        updateNeighbourList(i);
    }

    // int count = 0;
    for(int i=0; i<1000000; i++){
        demLoop();
   
        cycleCount++;
        if(cycleCount > 5000){
            printf("%lf\n",demTime/timeFactor);
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

    n2[0] = 8.0;
    n2[1] = 0.0;
    n2[2] = 0.0;

    n3[0] = 0.0;
    n3[1] = 6.0;
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
    */
/*
    for(int i=0; i<5; i++){
        demLoop();
    }
    // demSave();
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
