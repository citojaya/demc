#include "common.h"

/* 
DEM simulation
At every DPM time step this method is called by DEFINE_DPM_SOURCE
*/
void demLoop(){
    demTime += timeStep;
    for(int i=0; i<np; i++){
    //Find contact forces
        forceCalculation(i);
    }
    
    for(int i=0; i<np; i++){
        //Update position
        updatePosition(i);
    }

    for(int i=0; i<np; i++){
        //Update position
        updateNeighbourList(i);  
        demPart[i].currentTime += demPart[i].dt;     
    }
    
}

/* Assign graviational force*/
void assignGravity(){
    for(int i=0; i<np; i++){
        double gForce = demPart[i].mass; //gravitational acceleration
        demPart[i].force[1] = -gForce; // vertical force component
    }
}

/* Find contact force */
// void updateForce(Tracked_Particle *p){
//     demPart[p].force[0] = 0.0;
//     demPart[p].force[1] = -P_MASS(p)*massFactor; //Gravitional force
//     demPart[p].force[2] = 0.0;
// }


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
    for(int i=0; i<demPart[pI].noOfNeigh; i++){
        //printf("Neighbours\n");
        int jp = demPart[pI].neigh[i];
        //if(contactNotDone(pI, jp)){
            double gap = getCenterDist(pI,jp)-(demPart[pI].dia+demPart[jp].dia)*0.5;
            if(gap < 0.0){
                partContactForce(pI,jp, -gap);
                demPart[pI].contDone[demPart[pI].contacts] = jp; //record contact on jp
                demPart[pI].contacts++;
                demPart[jp].contDone[demPart[jp].contacts] = pI; //record contact on jp
                demPart[jp].contacts++;
            }
            double vGap = gap; //use to calcualte vanderwal force
        
            if(gap < (demPart[pI].dia+demPart[jp].dia)*0.05){ppVWForce(pI, jp, gap);}
        //}
    }
}

/*Particle-particle vanderwal force*/
void ppVWForce(int ip, int jp, double vGap){
    
    double ipDia = demPart[ip].dia;
    double jpDia = demPart[jp].dia;
    //double vGapMn = 0.5*(ipDia+jpDia)*0.01;
    double vGapMn = 1.0e-9*lengthFactor;
    
    double ijHa = sqrt(demPart[ip].ha*demPart[jp].ha);
    if(vGap < vGapMn){vGap = vGapMn;}
    double fv = 0;
    //if(vGap < (ipDia+jpDia)*0.01){
    fv = -ijHa*pow((ipDia*jpDia),3)*(vGap + 0.5*(ipDia + jpDia))
        /pow(((pow(vGap,2) + vGap*ipDia + vGap*jpDia)*(pow(vGap,2) + ipDia*vGap + jpDia*vGap + ipDia*jpDia)),2);       
    //}
    double *uVec = allocateDoubleArray(DIM);
    vecSub(demPart[ip].pos, demPart[jp].pos, uVec);
    unitVec(uVec, uVec);

    //writeLogNum("logfile2.log","VWForce ",fv/(lengthFactor*forceFactor));
    demPart[ip].force[0] += uVec[0]*fv;
    demPart[ip].force[1] += uVec[1]*fv;
    demPart[ip].force[2] += uVec[2]*fv;
    
    //demPart[jp].force[0] += -uVec[0]*fv;
    //demPart[jp].force[1] += -uVec[1]*fv;
    //demPart[jp].force[2] += -uVec[2]*fv;

   
    free(uVec);
   /* if(fabs(fv) > 1e-6){
        writeLogNum("logfile2.log", "VW force ", fv/(forceFactor*lengthFactor));
    }*/
}

/*Particle-wall vanderwal force*/
void pWallVWForce(int p, double vGap, double *uVec){
    double vGapMn = 1.0e-9*lengthFactor;
    if(vGap < vGapMn){vGap = vGapMn;}
    double fv = -demPart[p].ha*pow(demPart[p].dia,3)*0.5/pow((vGap*(vGap+demPart[p].dia)),2);

    demPart[p].force[0] += uVec[0]*fv;
    demPart[p].force[1] += uVec[1]*fv;
    demPart[p].force[2] += uVec[2]*fv;   
}

void checkXContact(int p, double xMin, double xMax){
    //Contact with xMin
    double gap = demPart[p].pos[0] - xMin - demPart[p].dia*0.5;
    uVec[0] = 1.0;
    uVec[1] = 0.0;
    uVec[2] = 0.0;  
    
    if(gap < 0) //If contact exists calculate contact force
    {
       surfaceContactForce(p, -gap, uVec);
       if(gap < demPart[p].dia*0.05){pWallVWForce(p, gap, uVec);}
    }  

    //Contact with xMax
    gap = xMax - demPart[p].pos[0] - demPart[p].dia*0.5;
    uVec[0] = -1.0;
    uVec[1] = 0.0;
    uVec[2] = 0.0;  
    if(gap < 0) //If contact exists calculate contact force
    {
       surfaceContactForce(p, -gap, uVec);
       if(gap < demPart[p].dia*0.05){pWallVWForce(p, gap, uVec);}
    }  
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
       if(gap < demPart[p].dia*0.05){pWallVWForce(p, gap, uVec);}
    }  

    //Contact with z=88mm wall
    gap = zMax - demPart[p].pos[2] - demPart[p].dia*0.5;
    uVec[0] = 0.0;
    uVec[1] = 0.0;
    uVec[2] = -1.0;  
    if(gap < 0) //If contact exists calculate contact force
    {
       surfaceContactForce(p, -gap, uVec);
       if(gap < demPart[p].dia*0.05){pWallVWForce(p, gap, uVec);}
    }  
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
        if(gap < demPart[p].dia*0.05){pWallVWForce(p, gap, uVec);}
    }  

    // Contact with yMax
    gap = yMax - (demPart[p].pos[1] + demPart[p].dia*0.5);
    uVec[0] = 0.0;
    uVec[1] = -1.0;
    uVec[2] = 0.0;  

    if(gap < 0) //If contact exists calculate contact force
    {
        surfaceContactForce(p, -gap, uVec);
        if(gap < demPart[p].dia*0.05){pWallVWForce(p, gap, uVec);}
    }
}

/*
Calculate particle-wall contact forces
param:
ip - ith particle
nrmDisp - overlap
uVec - unit vector normal to contact surface
*/
void surfaceContactForce(int p, double nrmDisp, double *uVec){
    double rStar = 0.5*demPart[p].dia;
    sclVecMult(-0.5*demPart[p].dia,uVec,ipRVec);
    
    crossProd(demPart[p].angVel,ipRVec,rotVel);
    vecAdd(demPart[p].vel,rotVel,ipCntPntVel);

    double *relVel = allocateDoubleArray(DIM);

    sclVecMult(1.0,demPart[p].vel,relVel);
    
    double nrmVel = dotProduct(relVel,uVec);
    free(relVel);
        
    sclVecMult(1.0,ipCntPntVel,cntPntVel);

    double *totalForce = allocateDoubleArray(DIM);
    double *momentum = allocateDoubleArray(DIM);

    double nrmCntForce = elasticMod*sqrt(rStar*nrmDisp)*nrmDisp;
    double nrmDampForce = -dmpn*elasticMod*sqrt(rStar*nrmDisp)*nrmVel;

    double *cntPntDisp = allocateDoubleArray(DIM);
    sclVecMult(timeStep,cntPntVel,cntPntDisp);
    double *tipCntPntDisp = allocateDoubleArray(DIM);
    projVec(cntPntDisp, uVec, tipCntPntDisp, 0);
    double slidingDisp = vecMag(tipCntPntDisp);

    double *disp = allocateDoubleArray(DIM);
    projVec(demPart[p].hisDisp, uVec, disp, 1);
    vecAdd(disp, tipCntPntDisp, disp);
    dsmax = dsmaxCff*nrmDisp;
    dd = vecMag(disp);
    fdt = sfc*nrmCntForce;
    dti = vecMag(tipCntPntDisp);

    double *fdtVec = allocateDoubleArray(DIM);
    //double *tngUVec = allocateDoubleArray(DIM);

    // unitVec(disp, tngUVec);
    // tngUVec[0] = -tngUVec[0];
    // tngUVec[1] = -tngUVec[1];
    // tngUVec[2] = -tngUVec[2];
    //sclVecMult(fdt,tngUVec,fdtVec);
    
    if(dd < 1e-6){
        dd = 0.0;
    }
    if(dti < 1e-6){
        dti = 0.0;
    }
    if(dd >= dsmax){
        sclMult(dsmax/dd,disp);
        dti = vecMag(tipCntPntDisp);
        if(dti != 0){
            sclVecMult(fdt/dti,tipCntPntDisp,fdtVec);
        }
        else{
            sclMult(0.0,fdtVec);
        }
    }
    else{
        if(dd != 0.0){
            fdt = fdt*(1.0 - pow((1.0 - dd/dsmax),1.5));
            sclVecMult(-fdt/dd,disp,fdtVec);
        }
        else{
            sclMult(0.0,fdtVec);
        }
    }

    sclVecMult(1.0,disp, demPart[p].hisDisp);

    //sum of forces
    double nrmForce = (nrmCntForce + nrmDampForce);
    sclVecMult(nrmForce, uVec, totalForce);

    vecAdd(totalForce, fdtVec, totalForce);
    crossProd(ipRVec, totalForce, momentum);
    double *rotMom = allocateDoubleArray(DIM);
    sclVecMult(0.5*rf*demPart[p].dia*nrmCntForce, demPart[p].angVel, rotMom);
    vecAdd(momentum, rotMom, momentum);

    //Update force and momentum on particle
    vecAdd(demPart[p].force, totalForce, demPart[p].force);
    vecAdd(demPart[p].momentum, momentum, demPart[p].momentum);

    free(rotMom);   
    free(fdtVec);
    free(tipCntPntDisp);
    free(cntPntDisp);
    free(totalForce);
    free(momentum);
    free(disp);
}

/*
Calculate interparticle forces
param:
ip - ith particle
jp - neighbour particle
nrmDisp - overlap
*/
void partContactForce(int ip, int jp, double nrmDisp){
    double rStar = 0.5*demPart[ip].dia*demPart[jp].dia/(demPart[ip].dia+demPart[jp].dia);
    vecSub(demPart[ip].pos, demPart[jp].pos, ijVec);
    unitVec(ijVec,uVec);

    //writeLog("logfile2.log","GAP ",nrmDisp);
    //double *relVel = allocateDoubleArray(DIM);
    //double *pVel = allocateDoubleArray(DIM);

    sclVecMult(-0.5*demPart[ip].dia,uVec,ipRVec);
    sclVecMult(0.5*demPart[jp].dia,uVec,jpRVec);

    crossProd(demPart[ip].angVel,ipRVec,rotVel);
    vecAdd(demPart[ip].vel,rotVel,ipCntPntVel);

    crossProd(demPart[jp].angVel,jpRVec,rotVel);
    vecAdd(demPart[jp].vel,rotVel,jpCntPntVel);

    double *relVel = allocateDoubleArray(DIM);
    vecSub(demPart[ip].vel,demPart[jp].vel,relVel);

    double nrmVel = dotProduct(relVel,uVec);
    vecSub(ipCntPntVel,jpCntPntVel,cntPntVel);
    free(relVel);
 
    double *totalForce = allocateDoubleArray(DIM);
    double *momentum = allocateDoubleArray(DIM);

    double nrmCntForce = elasticMod*sqrt(rStar*nrmDisp)*nrmDisp;
    double nrmDampForce = -dmpn*elasticMod*sqrt(rStar*nrmDisp)*nrmVel;

    double *cntPntDisp = allocateDoubleArray(DIM);
    sclVecMult(timeStep,cntPntVel,cntPntDisp);
    double *tipCntPntDisp = allocateDoubleArray(DIM);
    projVec(cntPntDisp, uVec, tipCntPntDisp, 0);
    double slidingDisp = vecMag(tipCntPntDisp);

    double *disp = allocateDoubleArray(DIM);

    sclVecMult(1.0,tipCntPntDisp,disp);

    //projVec(demPart[ip].hisDisp, uVec, disp, 1);
    //vecAdd(disp, tipCntPntDisp, disp);
    dsmax = dsmaxCff*nrmDisp;
    dd = vecMag(disp);
    dti = vecMag(tipCntPntDisp);
    fdt = sfc*nrmCntForce;

    double *fdtVec = allocateDoubleArray(DIM);
 
    if(dd < 1e-6){
        dd = 0.0;
    }
    if(dti < 1e-6){
        dti = 0.0;
    }
    if(dd >= dsmax){
        sclMult(dsmax/dd,disp);
        if(dti != 0){
            sclVecMult(-fdt/dti,tipCntPntDisp,fdtVec);
        }
        else{
            sclMult(0.0,fdtVec);
        }
    }
    else{
        //sclVecMult(1.0,disp, demPart[ip].hisDisp);
        if(dd != 0.0){
            fdt = fdt*(1.0 - pow((1.0 - dd/dsmax),1.5));
            sclVecMult(-fdt/dd,disp,fdtVec);
        }
        else{
            sclMult(0.0,fdtVec);
        }
    }

     //sum of forces
    double nrmForce = (nrmCntForce + nrmDampForce);
    //writeLog("logfile2.log","nrmForce ",nrmForce);
    sclVecMult(nrmForce, uVec, totalForce);
    vecAdd(totalForce, fdtVec, totalForce);
    crossProd(ipRVec, totalForce, momentum);
    double *rotMom = allocateDoubleArray(DIM);
    sclVecMult(0.5*rf*demPart[ip].dia*nrmCntForce, demPart[ip].angVel, rotMom);
    vecAdd(momentum, rotMom, momentum);

    //Update force and momentum on particle
    vecAdd(demPart[ip].force, totalForce, demPart[ip].force);
    vecAdd(demPart[ip].momentum, momentum, demPart[ip].momentum);


    //Add forces to jp particle
    sclVecMult(-1.0, totalForce, totalForce);
    sclVecMult(-1.0, momentum, momentum);
    //vecAdd(demPart[jp].force, totalForce, demPart[jp].force);
    //vecAdd(demPart[jp].momentum, momentum, demPart[jp].momentum);   


    free(rotMom);
    free(fdtVec);
    free(tipCntPntDisp);
    free(cntPntDisp);
    free(totalForce);
    free(momentum);
    free(disp);
    
}


/* Convert DEM Lagrangian values to Euler values & then makes a copy to Fluent 
secondary phase
*/
void copyDEMInfo(){
    printf("copyDEMInfo %d, %d, %d\n",xDiv,yDiv,zDiv);
}


/* Update particle position*/
void forceCalculation(int p)
{
    double dt = timeStep;
    //writeLog("logfile2.log","MOVE DEM PART time ",P_TIME(p));

    demPart[p].force[0] = 0.0;
    demPart[p].force[1] = 0.0;
    demPart[p].force[2] = -demPart[p].mass; //gravitational force

    //Add Bouynacy

    int nearestFaceIndex = -1;
    int nearestFaceCellIndex = -1;
    double minDist = 100000.0;
    double ipX = demPart[p].pos[0];
    double ipY = demPart[p].pos[1];
    double ipZ = demPart[p].pos[2];

    int iIndex = ceil((demPart[p].pos[0]-xmin)/domainDx);
    int jIndex = ceil((demPart[p].pos[1]-ymin)/domainDy);
    int kIndex = ceil((demPart[p].pos[2]-zmin)/domainDz);
 /* 
    for(int rr=kIndex-1; rr<kIndex+2; rr++){
        for(int qq=jIndex-1; qq<jIndex+2; qq++){
            for(int pp=iIndex-1; pp<iIndex+2; pp++){
                int neighCellIndex = pp + qq*xDiv + rr*xDiv*yDiv;
                for(int i=0; i<bdBox[neighCellIndex].totalFaces; i++){
                    double jpX = bdBox[neighCellIndex].face[i].centroid[0]*lengthFactor;
                    double jpY = bdBox[neighCellIndex].face[i].centroid[1]*lengthFactor;
                    double jpZ = bdBox[neighCellIndex].face[i].centroid[2]*lengthFactor;

                    double dist = sqrt((ipX-jpX)*(ipX-jpX) + (ipY-jpY)*(ipY-jpY) + (ipZ-jpZ)*(ipZ-jpZ));
                    if(dist < minDist){
                        nearestFaceIndex = i;
                        nearestFaceCellIndex = neighCellIndex;
                        minDist = dist;
                    }
                }
            }
        }
    }   


    if(nearestFaceIndex != -1){
        demPart[p].faceNode1[0] = bdBox[nearestFaceCellIndex].face[nearestFaceIndex].node1[0]*lengthFactor;
        demPart[p].faceNode1[1] = bdBox[nearestFaceCellIndex].face[nearestFaceIndex].node1[1]*lengthFactor;
        demPart[p].faceNode1[2] = bdBox[nearestFaceCellIndex].face[nearestFaceIndex].node1[2]*lengthFactor;

        demPart[p].faceNode2[0] = bdBox[nearestFaceCellIndex].face[nearestFaceIndex].node2[0]*lengthFactor;
        demPart[p].faceNode2[1] = bdBox[nearestFaceCellIndex].face[nearestFaceIndex].node2[1]*lengthFactor;
        demPart[p].faceNode2[2] = bdBox[nearestFaceCellIndex].face[nearestFaceIndex].node2[2]*lengthFactor;

        demPart[p].faceNode3[0] = bdBox[nearestFaceCellIndex].face[nearestFaceIndex].node3[0]*lengthFactor;
        demPart[p].faceNode3[1] = bdBox[nearestFaceCellIndex].face[nearestFaceIndex].node3[1]*lengthFactor;
        demPart[p].faceNode3[2] = bdBox[nearestFaceCellIndex].face[nearestFaceIndex].node3[2]*lengthFactor;

        //writeLogNum("logfile2.log", "nearestFaceIndex ", cI);
            //double gap = demPart[p].pos[1] - demPart[p].dia*0.5;
        double *uVec = allocateDoubleArray(DIM);
        //double *tempN1 = allocateDoubleArray(DIM);
        getUnitVector(demPart[p].faceNode1,demPart[p].faceNode2,
                            demPart[p].faceNode3,uVec);
 
        double gap = getOverlap(demPart[p].pos,demPart[p].dia, demPart[p].faceNode1,  
                                demPart[p].faceNode2, demPart[p].faceNode2, uVec);
   
        // writeLogNum("logfile2.log","uVec1 ",uVec[0]);
        // writeLogNum("logfile2.log","uVec2 ",uVec[1]);
        // writeLogNum("logfile2.log","uVec3 ",uVec[2]);
        //writeLogNum("logfile2.log","gap ",gap/lengthFactor);
        if(gap < 0) //If contact exists calculate contact force
        {
            surfaceContactForce(p, -gap, uVec);
        } 
        free(uVec);
    */
    
//working code
    //check contact with Y face
    checkYContact(p, -0.005*lengthFactor, 0.005*lengthFactor);

    //check top bound contact
    if(demPart[p].pos[2] > (0.05*lengthFactor-0.5*demPart[p].dia)){
        checkZContact(p, 0, 0.05*lengthFactor);
    }
    else if(demPart[p].pos[0] <= 0.1*lengthFactor || 
        demPart[p].pos[0] >= 0.15*lengthFactor){
        checkZContact(p, 0,0.05*lengthFactor);
    }
    else if(demPart[p].pos[2] <= 0){
        checkXContact(p, 0.1*lengthFactor, 0.15*lengthFactor);
        checkZContact(p, -0.05*lengthFactor, 0.05*lengthFactor);
    }
 
    else{    
        //check for x=100 edge
        if(demPart[p].pos[0] > 0.1*lengthFactor && 
            demPart[p].pos[0] < 0.1*lengthFactor+0.5*demPart[p].dia){
            if(demPart[p].pos[2] > 0 && demPart[p].pos[2] < 0.5*demPart[p].dia){
                double *iVec = allocateDoubleArray(DIM);
                double *jVec = allocateDoubleArray(DIM);
                
                iVec[0] = demPart[p].pos[0];
                iVec[1] = 0.0;
                iVec[2] = demPart[p].pos[2];
                jVec[0] = 0.1*lengthFactor;
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
        else if(demPart[p].pos[0] < 0.15*lengthFactor && 
            demPart[p].pos[0] > 0.15*lengthFactor-0.5*demPart[p].dia){
            
            if(demPart[p].pos[2] > 0 && demPart[p].pos[2] < 0.5*demPart[p].dia){
                double *iVec = allocateDoubleArray(DIM);
                double *jVec = allocateDoubleArray(DIM);
                
                iVec[0] = demPart[p].pos[0];
                iVec[1] = 0.0;
                iVec[2] = demPart[p].pos[2];
                jVec[0] = 0.15*lengthFactor;
                jVec[1] = 0.0;//demPart[p].pos[1];
                jVec[2] = 0.0;
       
                double *ijVec = allocateDoubleArray(DIM);
                vecSub(iVec, jVec, ijVec);
                double gap = vecMag(ijVec);
                unitVec(ijVec, uVec);

                if(gap < 0.5*demPart[p].dia){
                    //writeLogNum("logfile2.log",  "EDGE 2 OVERLAP ",gap);
                    surfaceContactForce(p, (0.5*demPart[p].dia-gap), uVec);
               }
                free(ijVec);
                free(iVec);
                free(jVec); 
                 
            }
             
        }
    }

    //Find particle-particle contact force and particle-particle Vanderwal force
    neighbourContactForce(p);
}

void updatePosition(int p){
    
    //demPart[p].noOfCntF = 0; //reset
    //demPart[p].insertable = 1; //reset
    double dxDot = demPart[p].force[0]*timeStep/demPart[p].mass;
    double dyDot = demPart[p].force[1]*timeStep/demPart[p].mass;
    double dzDot = demPart[p].force[2]*timeStep/demPart[p].mass;
    
    demPart[p].vel[0] += dxDot;
    demPart[p].vel[1] += dyDot;
    demPart[p].vel[2] += dzDot;

    double dx = demPart[p].vel[0]*timeStep;
    double dy = demPart[p].vel[1]*timeStep;
    double dz = demPart[p].vel[2]*timeStep;

    demPart[p].pos[0] += dx;
    demPart[p].pos[1] += dy;
    demPart[p].pos[2] += dz;
 
    demPart[p].currentTime += timeStep;

    demPart[p].displacement += sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
    demPart[p].contacts = 0; //reset neighbour contact calculation to zero

}





