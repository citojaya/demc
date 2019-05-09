#include "common.h"


/*Calculate overlap with contact surface
param:
parPos - particle center
dia - particle diameter
node1, node2, node3 - contact surface nodes
uVec - unit vector normal to surface
*/
double getOverlap(double *parPos, double dia, double *n1, double *n2, double *n3, double *uVec){ 
    double *v1 = allocateDoubleArray(DIM); //vector running from node1 to particles center
    double *v2 = allocateDoubleArray(DIM); //vector running from node1 to particles center
    double *ppDash = allocateDoubleArray(DIM); //vector running from particle center to projection
    vecSub(parPos,n1,v1);  
    projVec(v1, uVec, v2, 0);
    vecSub(v1,v2,ppDash);
    uVec[0] = ppDash[0];
    uVec[1] = ppDash[1];
    uVec[2] = ppDash[2];

    unitVec(uVec,uVec);
    double overlap = vecMag(ppDash)-0.5*dia;

    free(v1);
    free(v2);
    free(ppDash);
    return overlap;
}


/*Find contact forces with all neighbour particles*/
void neighbourContactForce(int pI){
    //Loop through all neighbour particles
    int count = 0;
    
    for(int i=0; i<demPart[pI].noOfNeigh; i++){
        //printf("Neighbours\n");
        int jp = demPart[pI].neigh[i];
        double gap = getCenterDist(pI,jp)-(demPart[pI].dia+demPart[jp].dia)*0.5;
        //if(pI < jp){
            // if(gap < demPart[pI].dia){
            //      ppElectForce(pI, jp, getCenterDist(pI,jp));
            // }
            if(gap < 0.0){
                count += 1;
                partContactForce(pI,jp, -gap);
                //partContactForce(jp,pI, -gap);
            }
            if(gap < 100.e-9*lengthFactor)//activate vanderwaal force when gap<100nm
            {
                //ppVWForce(pI, jp, gap);
                //ppVWForce(jp, pI, gap);
            }
        //}
    }
    //demPart[pI].cordNo = count;
}


/*Find particle-wall contact from solid boundary*/
void findContactFromBoundary(int p){
    
//working code
    checkYContact(p, ductymin*lengthFactor, ductymax*lengthFactor);
    //check top bound contact
    if(demPart[p].pos[2] > (ductzmax*lengthFactor-0.5*demPart[p].dia)){
        checkZContact(p, 0, ductzmax*lengthFactor);
    }
    else if(demPart[p].pos[0] <= ductxedge1*lengthFactor || 
        demPart[p].pos[0] >= ductxedge2*lengthFactor){
        checkZContact(p, 0,ductzmax*lengthFactor);
    }
    else if(demPart[p].pos[2] <= 0){
         checkXContact(p, ductxedge1*lengthFactor, ductxedge2*lengthFactor);
         checkZContact(p, ductzmin*lengthFactor, ductzmax*lengthFactor);
    }
 
    else{    
        //check for x=100 edge
        if(demPart[p].pos[0] > ductxedge1*lengthFactor && 
            demPart[p].pos[0] < ductxedge1*lengthFactor+0.5*demPart[p].dia){
            if(demPart[p].pos[2] > 0 && demPart[p].pos[2] < 0.5*demPart[p].dia){
                double *iVec = allocateDoubleArray(DIM);
                double *jVec = allocateDoubleArray(DIM);
                
                iVec[0] = demPart[p].pos[0];
                iVec[1] = 0.0;
                iVec[2] = demPart[p].pos[2];
                jVec[0] = ductxedge1*lengthFactor;
                jVec[1] = 0.0;//demPart[p].pos[1];
                jVec[2] = 0.0;
            
                double *ijVec = allocateDoubleArray(DIM);
                vecSub(iVec, jVec, ijVec);
                double gap = vecMag(ijVec);
                unitVec(ijVec, uVec);

                if(gap < 0.5*demPart[p].dia){
                    surfaceContactForce(p, (0.5*demPart[p].dia-gap), uVec);
                }
                free(ijVec);
                free(iVec);
                free(jVec);                 
            }
        }
     
        //check for x=150 edge
        else if(demPart[p].pos[0] < ductxedge2*lengthFactor && 
            demPart[p].pos[0] > ductxedge2*lengthFactor-0.5*demPart[p].dia){
            
            if(demPart[p].pos[2] > 0 && demPart[p].pos[2] < 0.5*demPart[p].dia){
                double *iVec = allocateDoubleArray(DIM);
                double *jVec = allocateDoubleArray(DIM);
                
                iVec[0] = demPart[p].pos[0];
                iVec[1] = 0.0;
                iVec[2] = demPart[p].pos[2];
                jVec[0] = ductxedge2*lengthFactor;
                jVec[1] = 0.0;//demPart[p].pos[1];
                jVec[2] = 0.0;
       
                double *ijVec = allocateDoubleArray(DIM);
                vecSub(iVec, jVec, ijVec);
                double gap = vecMag(ijVec);
                unitVec(ijVec, uVec);

                if(gap < 0.5*demPart[p].dia){
                    surfaceContactForce(p, (0.5*demPart[p].dia-gap), uVec);
                 }
                free(ijVec);
                free(iVec);
                free(jVec); 
            }
        }
    }
}

/*Find particle-wall contact and calculate force using trianglusr boundary Mesh*/
/*void findContactFromMesh(int p){
    int cI = demPart[p].prevCellIndex;
    int nearestFaceIndex = -1;
    double minDist = 100000;
    double ipX = demPart[p].pos[0];
    double ipY = demPart[p].pos[1];
    double ipZ = demPart[p].pos[2];
    writeLogNum("logfile2.log", "particle cell index ", cI);
    for(int i=0; i<bdBox[cI].totalFaces; i++){
        
        double jpX = bdBox[cI].face[i].centroid[0]*lengthFactor;
        double jpY = bdBox[cI].face[i].centroid[1]*lengthFactor;
        double jpZ = bdBox[cI].face[i].centroid[2]*lengthFactor;

        double dist = sqrt((ipX-jpX)*(ipX-jpX) + (ipY-jpY)*(ipY-jpY) + (ipZ-jpZ)*(ipZ-jpZ));
        if(dist < minDist){
            nearestFaceIndex = i;
            minDist = dist;
        }
    }
    if(nearestFaceIndex != -1){
        bdBox[cI].face[nearestFaceIndex].node1[0] = bdBox[cI].face[nearestFaceIndex].node1[0]*lengthFactor;
        bdBox[cI].face[nearestFaceIndex].node1[1] = bdBox[cI].face[nearestFaceIndex].node1[1]*lengthFactor;
        bdBox[cI].face[nearestFaceIndex].node1[2] = bdBox[cI].face[nearestFaceIndex].node1[2]*lengthFactor;

        double *uVec = allocateDoubleArray(DIM);
        getUnitVector(bdBox[cI].face[nearestFaceIndex].node1,bdBox[cI].face[nearestFaceIndex].node2,
                            bdBox[cI].face[nearestFaceIndex].node3,uVec);

        double gap = getOverlap(demPart[p].pos,demPart[p].dia, bdBox[cI].face[nearestFaceIndex].node1,  
                                bdBox[cI].face[nearestFaceIndex].node2, bdBox[cI].face[nearestFaceIndex].node3, uVec);
   
        if(gap < 0) //If contact exists calculate contact force
        {
             surfaceContactForce(p, -gap, uVec);
        } 
        free(uVec);            
    }
}*/

void checkXContact(int p, double xMin, double xMax){
    //Contact with xMin
    double gap = demPart[p].pos[0] - xMin - demPart[p].dia*0.5;
    uVec[0] = 1.0;
    uVec[1] = 0.0;
    uVec[2] = 0.0;  
    
    if(gap < 0) //If contact exists calculate contact force
    {
       surfaceContactForce(p, -gap, uVec);
    }  
    if(gap < 100.e-9*lengthFactor){pWallVWForce(p, gap, uVec);}
    //Contact with xMax
    gap = xMax - demPart[p].pos[0] - demPart[p].dia*0.5;
    uVec[0] = -1.0;
    uVec[1] = 0.0;
    uVec[2] = 0.0;  
    if(gap < 0) //If contact exists calculate contact force
    {
       surfaceContactForce(p, -gap, uVec);
    }
    if(gap < 100.e-9*lengthFactor){pWallVWForce(p, gap, uVec);}
    //if(gap < demPart[p].dia*0.05){pWallVWForce(p, gap, uVec);}  
}

void checkZContact(int p, double zMin, double zMax){
    //Contact with zMin
    double gap = -(zMin - demPart[p].pos[2]) - demPart[p].dia*0.5;
    uVec[0] = 0.0;
    uVec[1] = 0.0;
    uVec[2] = 1.0;  
    if(gap < 0) //If contact exists calculate contact force
    {
       surfaceContactForce(p, -gap, uVec);
    }  
    if(gap < 100.e-9*lengthFactor){pWallVWForce(p, gap, uVec);}
    //Contact with z=88mm wall
    gap = zMax - demPart[p].pos[2] - demPart[p].dia*0.5;
    uVec[0] = 0.0;
    uVec[1] = 0.0;
    uVec[2] = -1.0;  
    if(gap < 0) //If contact exists calculate contact force
    {
       surfaceContactForce(p, -gap, uVec);
    } 
    if(gap < 100.e-9*lengthFactor){pWallVWForce(p, gap, uVec);} 
}

void checkYContact(int p, double yMin, double yMax){
    // Contact with yMin
    double gap = -(yMin - demPart[p].pos[1]) - demPart[p].dia*0.5;
    uVec[0] = 0.0;
    uVec[1] = 1.0;
    uVec[2] = 0.0;  

    if(gap < 0) //If contact exists calculate contact force
    {
        //demPart[p].force[1] += -demPart[p].mass;
        surfaceContactForce(p, -gap, uVec);
    }
    if(gap < 100.e-9*lengthFactor){pWallVWForce(p, gap, uVec);}  

    // Contact with yMax
    gap = yMax - (demPart[p].pos[1] + demPart[p].dia*0.5);
    uVec[0] = 0.0;
    uVec[1] = -1.0;
    uVec[2] = 0.0;  

    if(gap < 0) //If contact exists calculate contact force
    {
        surfaceContactForce(p, -gap, uVec);
    }
    if(gap < 100.e-9*lengthFactor){pWallVWForce(p, gap, uVec);}
}




