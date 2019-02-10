#include "common.h"

/* 
DEM simulation
At every DPM time step this method is called by DEFINE_DPM_SOURCE
*/
void demLoop(){
    demTime += timeStep;
    for(int i=0; i<np; i++){
        //deleteParticle(i);
        forceCalculation(i);
        updatePosition(i);
        //addToBdBox(i);
        demPart[i].currentTime += demPart[i].dt;
    }
    addToBdBox();
    
    for(int i=0; i<np; i++){
        if(demPart[i].displacement > allowedDisp){
            updateNeighbourList(i);
        }
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
//     demPart[p->part_id].force[0] = 0.0;
//     demPart[p->part_id].force[1] = -P_MASS(p)*massFactor; //Gravitional force
//     demPart[p->part_id].force[2] = 0.0;
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
        double gap = getCenterDist(pI,jp)-(demPart[pI].dia+demPart[jp].dia)*0.5;
        if(gap < 0.0){
            printf("OVERLAP %lf\n",gap*1e3/lengthFactor);
            partContactForce(pI,jp, -gap);
        }
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
    //double rStar = 0.5*P_DIAM(p)*lengthFactor;
    double rStar = 0.5*demPart[p].dia;
    sclVecMult(-0.5*demPart[p].dia,uVec,ipRVec);
    
    crossProd(demPart[p].angVel,ipRVec,rotVel);
    vecAdd(demPart[p].vel,rotVel,ipCntPntVel);

    double *relVel = allocateDoubleArray(DIM);
    // double *pVel = allocateDoubleArray(DIM);
    // pVel[0] = P_VEL(p)[0]*velocityFactor;
    // pVel[1] = P_VEL(p)[1]*velocityFactor;
    // pVel[2] = P_VEL(p)[2]*velocityFactor;

    
    //sclVecMult(1.0,pVel,relVel);
    sclVecMult(1.0,demPart[p].vel,relVel);
    
    double nrmVel = dotProduct(relVel,uVec);
    free(relVel);
    //free(pVel);
    
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

    // writeLogNum("logfile2.log"," disp X ",disp[0]);
    // writeLogNum("logfile2.log"," disp Y ",disp[1]);
    // writeLogNum("logfile2.log"," disp Z ",disp[2]);

    double *fdtVec = allocateDoubleArray(DIM);
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

    // writeLogNum("logfile2.log"," FDT X ",fdtVec[0]);
    // writeLogNum("logfile2.log"," FDT Y ",fdtVec[1]);
    // writeLogNum("logfile2.log"," FDT Z ",fdtVec[2]);
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
    projVec(demPart[ip].hisDisp, uVec, disp, 1);
    vecAdd(disp, tipCntPntDisp, disp);
    dsmax = dsmaxCff*nrmDisp;
    dd = vecMag(disp);
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

    sclVecMult(1.0,disp, demPart[ip].hisDisp);
   
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
    //Thread *tc = P_CELL_THREAD(p);
    //cell_t c = P_CELL(p);
    //demPart[p->part_id].force[1] += partVol(p->part_id)*C_R(c,tc)*densityFactor;

    //Find contact force with wall

    //Domain *d = Get_Domain(1);
    //Thread *t = Lookup_Thread(d,0);

    //Thread **tph = THREAD_SUB_THREADS(t);
    //writeLog("logfile2.log","CELL ",C_VOF(P_CELL(p),t));
    //int i;
    //face_t f;

    
    //writeLogNum("logfile2.log","NO OF WALLS ",noOfWalls);
    //writeLogNum("logfile2.log","WALLS ",walls[0]);
    /*
    for(int k=0; k<noOfWalls; k++)
    {
        //if(t_id == walls[k])
 
            c_face_loop(c,tc,i)//Loop over faces
            {
                int t_id = THREAD_ID(C_FACE_THREAD(c,tc,i));

                if(t_id == walls[k])
                {
        // if(p->part_id == 1 && t_id == 5)
        // {
        //     writeLogNum("logfile2.log","FACE ID ",P_POS(p)[1]);
        // }
        // writeLogNum("logfile2.log","Y ",P_POS(p)[1]);
        // writeLogNum("logfile2.log","face thread ID ",t_id);
        

            //writeLogNum("logfile2.log","t_id ",t_id);
            // writeLogNum("logfile2.log","face thread ID ",t_id);
            //writeLogNum("logfile2.log","Y ",P_POS(p)[1]);
            // writeLogLine("logfile2.log","---------------\n");
            //writeLogNum("logfile2.log","face thread ID ",t_id);
            f = C_FACE(c,tc,i);
            Thread *tf = C_FACE_THREAD(c,tc,i);

            //int face_t_id = THREAD_ID(C_FACE_THREAD(cc,t,i));
            int n;
            int j=0;
            f_node_loop(f,tf,n)//Loop over face nodes
            {
            Node *node = F_NODE(f,tf,n);
            if(j==0)
            {
                demPart[p->part_id].faceNode1[0] = NODE_X(node)*lengthFactor;
                demPart[p->part_id].faceNode1[1] = NODE_Y(node)*lengthFactor;
                demPart[p->part_id].faceNode1[2] = NODE_Z(node)*lengthFactor;
            }
            if(j==1)
            {
                demPart[p->part_id].faceNode2[0] = NODE_X(node)*lengthFactor;
                demPart[p->part_id].faceNode2[1] = NODE_Y(node)*lengthFactor;
                demPart[p->part_id].faceNode2[2] = NODE_Z(node)*lengthFactor;
            }
            if(j==2)
            {
                demPart[p->part_id].faceNode3[0] = NODE_X(node)*lengthFactor; 
                demPart[p->part_id].faceNode3[1] = NODE_Y(node)*lengthFactor;
                demPart[p->part_id].faceNode3[2] = NODE_Z(node)*lengthFactor;
            }
            j++;
            }//end of loop over face nodes
            //double gap = demPart[p->part_id].pos[1] - demPart[p->part_id].dia*0.5;
            double *uVec = allocateDoubleArray(DIM);
            getUnitVector(demPart[p->part_id].faceNode1,demPart[p->part_id].faceNode2,
                            demPart[p->part_id].faceNode3,uVec);
            double *parPos = allocateDoubleArray(DIM);
            parPos[0] = P_POS(p)[0]*lengthFactor;
            parPos[1] = P_POS(p)[1]*lengthFactor;
            parPos[2] = P_POS(p)[2]*lengthFactor;
//            double gap = getOverlap(demPart[p->part_id].pos,demPart[p->part_id].dia, demPart[p->part_id].faceNode1,  
//                                demPart[p->part_id].faceNode2, demPart[p->part_id].faceNode3, uVec);
            double gap = getOverlap(parPos,demPart[p->part_id].dia, demPart[p->part_id].faceNode1,  
                                demPart[p->part_id].faceNode2, demPart[p->part_id].faceNode3, uVec);        
            //writeLog("logfile2.log","PAR TIME ",P_TIME(p));
            if(gap < 0) //If contact exists calculate contact force
            {
                //writeLogNum("logfile2.log", "PAR ", p->part_id);
                //
                //if(gap < -demPart[p->part_id].dia*0.05){gap = -demPart[p->part_id].dia*0.002;}
                //writeLogNum("logfile2.log", "gap ", 1e3*gap/lengthFactor);
                surfaceContactForce(p->part_id, -gap, uVec);
            } 
            free(uVec);
            free(parPos);
            }//if(t_id == walls[i])
        }//end of loop over faces
    }//end of nOfWalls
    */

    //Fluid drag force calculation
    // double velFX = C_U(c,tc)*velocityFactor;  
    // double velFY = C_V(c,tc)*velocityFactor;  
    // double velFZ = C_W(c,tc)*velocityFactor;    
    // double pGX = C_P_G(c,tc)[0]*pressureFactor/lengthFactor;
    // double pGY = C_P_G(c,tc)[1]*pressureFactor/lengthFactor;
    // double pGZ = C_P_G(c,tc)[2]*pressureFactor/lengthFactor;
    // double velPX = demPart[p->part_id].vel[0];
    // double velPY = demPart[p->part_id].vel[1];
    // double velPZ = demPart[p->part_id].vel[2];
    //double density = C_R(c,tc)*densityFactor;
    //double visc = VISCOSITY*massFactor/(lengthFactor*timeFactor);

    //double instPor = 1.0-solidFraction(p->part_id);//Solid fraction;
    //writeLogNum("logfile2.log", "PORO  ",instPor);
    //--- Test code ---------------------------------------------------------
    // double Re = instPor*demPart[p->part_id].dia*sqrt((velFX-velPX)*(velFX-velPX)  + (velFY-velPY)*(velFY-velPY) 
    //                 + (velFZ-velPZ)*(velFZ-velPZ))*density/visc;
    
    // double coeff, dCoeff,modPor;

    //-- Standard equations -------------------------------------------------
	// if (Re < 1e-6)//if Reynolds number is too small following parameters cannot be defined
	// {
	// 	coeff = 0.0;
	// 	dCoeff = 0.0;
	// 	modPor = 0.0;
	// }
	// else 
	// {
	// 	coeff = 3.7-0.65*exp(-(1.5-log10(Re))*(1.5-log10(Re))*0.5);
	// 	dCoeff = (0.63+4.8/sqrt(Re))*(0.63+4.8/sqrt(Re));
	// 	modPor = pow(instPor,(-coeff));
	// }

    // double pfFX = modPor*dCoeff*density*instPor*(velFX-velPX)*fabs(velFX-velPX)*PI*(pow(demPart[p->part_id].dia*0.5,2));        
    // double pfFY = modPor*dCoeff*density*instPor*(velFY-velPY)*fabs(velFY-velPY)*PI*(pow(demPart[p->part_id].dia*0.5,2));        
    // double pfFZ = modPor*dCoeff*density*instPor*(velFZ-velPZ)*fabs(velFZ-velPZ)*PI*(pow(demPart[p->part_id].dia*0.5,2));        

	//------------------------------------------------------------------------

    // double relVelMag = sqrt((velFX-velPX)*(velFX-velPX)+(velFY-velPY)*(velFY-velPY)+(velFZ-velPZ)*(velFZ-velPZ));
    // double Re = instPor*demPart[p->part_id].dia*relVelMag*density/visc;
    
    // double dragCoeff, beta;
    // if(Re <= 1000){
    //     dragCoeff = 0.5*(24.0/Re)*(1.0+0.15*pow(Re,0.687));
    // }
    // else if(Re > 1000){
    //     dragCoeff = 0.43;
    // }
    
    //writeLogNum("logfile2.log", "DRAG ", dragCoeff);
    // if(instPor <= 0.8){
    //     double b1 = (1.0-instPor)/(demPart[p->part_id].dia*pow(instPor, 2));
    //     double b2 = 150.0*(1.0-instPor)*visc/demPart[p->part_id].dia;
    //     double b3 = 1.75*density*instPor*relVelMag;
    //     beta = b1*(b2 + b3);
    // }
    // else if(instPor > 0.8){
    //     beta = (3.0/4.0)*dragCoeff*density*relVelMag*(1.0-instPor)*pow(instPor,-2.7)/demPart[p->part_id].dia;
    // }

    // double pfFX = (velFX-velPX)*partVol(p->part_id)*beta/(1.0-instPor);
    // double pfFY = (velFY-velPY)*partVol(p->part_id)*beta/(1.0-instPor);
    // double pfFZ = (velFZ-velPZ)*partVol(p->part_id)*beta/(1.0-instPor);

    // double pGFX = -pGX*PI*pow(demPart[p->part_id].dia,3)/6.0;
    // double pGFY = -pGY*PI*pow(demPart[p->part_id].dia,3)/6.0;
    // double pGFZ = -pGZ*PI*pow(demPart[p->part_id].dia,3)/6.0;

    // writeLogNum("logfile2.log","pGFX",pGFX);
    // writeLogNum("logfile2.log","pGFY",pGFY);
    // writeLogNum("logfile2.log","pGFZ",pGFZ);

    // writeLogLine("logfile2.log","--------------------");


    // writeLogNum("logfile2.log","pfFX",pfFX);
    // writeLogNum("logfile2.log","pfFY",pfFY);
    // writeLogNum("logfile2.log","pfFZ",pfFZ);

    // writeLogLine("logfile2.log","--------------------");

    // End of fluid drag force calculation
    


    double xMax = 0.055*lengthFactor;
    double zMax = 0.010*lengthFactor;
    double yMax = 0.0005*lengthFactor;
    double zMin = -0.005*lengthFactor;
    double xMin = 0.045*lengthFactor;
    double yMin =  -0.0005*lengthFactor;


    //Contact with xMin
    double gap = demPart[p].pos[0] - xMin - demPart[p].dia*0.5;
    uVec[0] = 1.0;
    uVec[1] = 0.0;
    uVec[2] = 0.0;  
    
    if(gap < 0) //If contact exists calculate contact force
    {
       surfaceContactForce(p, -gap, uVec);
    }  

    //Contact with xMax
    gap = xMax - demPart[p].pos[0] - demPart[p].dia*0.5;
    uVec[0] = -1.0;
    uVec[1] = 0.0;
    uVec[2] = 0.0;  
    if(gap < 0) //If contact exists calculate contact force
    {
       surfaceContactForce(p, -gap, uVec);
    }  

    //Contact with zZMin
    gap = -(zMin - demPart[p].pos[2]) - demPart[p].dia*0.5;
    uVec[0] = 0.0;
    uVec[1] = 0.0;
    uVec[2] = 1.0;  
    if(gap < 0) //If contact exists calculate contact force
    {
       surfaceContactForce(p, -gap, uVec);
    }  
/*
    //Contact with z=88mm wall
    gap = zMax - demPart[p->part_id].pos[2] - demPart[p->part_id].dia*0.5;
    uVec[0] = 0.0;
    uVec[1] = 0.0;
    uVec[2] = -1.0;  
    if(gap < 0) //If contact exists calculate contact force
    {
       //writeLogNum("logfile2.log","GAP ",gap/lengthFactor);
       surfaceContactForce(p->part_id, -gap, uVec);
    }  
*/
    // Contact with yMin
    gap = -(yMin - demPart[p].pos[1]) - demPart[p].dia*0.5;
    uVec[0] = 0.0;
    uVec[1] = 1.0;
    uVec[2] = 0.0;  

    if(gap < 0) //If contact exists calculate contact force
    {
        //demPart[p->part_id].force[1] += -demPart[p->part_id].mass;
        surfaceContactForce(p, -gap, uVec);
    }  

    // Contact with yYMax
    gap = yMax - (demPart[p].pos[1] + demPart[p].dia*0.5);
    uVec[0] = 0.0;
    uVec[1] = -1.0;
    uVec[2] = 0.0;  

    if(gap < 0) //If contact exists calculate contact force
    {
        surfaceContactForce(p, -gap, uVec);
    }

    //Find particle-particle contact force
    neighbourContactForce(p);

    // writeLogNum("logfile2.log","pfFY ",pfFY);
    // writeLogNum("logfile2.log","pGFY ",pGFY);
    // demPart[p].force[0]  += pfFX + pGFX;
    // demPart[p].force[1]  += pfFY + pGFY + partVol(p)*density;
    // demPart[p].force[2]  += pfFZ + pGFZ;


}

void updatePosition(int p){
    
    double dxDot = demPart[p].force[0]*timeStep/demPart[p].mass;
    double dyDot = demPart[p].force[1]*timeStep/demPart[p].mass;
    double dzDot = demPart[p].force[2]*timeStep/demPart[p].mass;
    
    demPart[p].vel[0] += dxDot;
    demPart[p].vel[1] += dyDot;
    demPart[p].vel[2] += dzDot;

    demPart[p].pos[0] += demPart[p].vel[0]*timeStep;
    demPart[p].pos[1] += demPart[p].vel[1]*timeStep;
    demPart[p].pos[2] += demPart[p].vel[2]*timeStep;
   
    demPart[p].currentTime += timeStep;

    demPart[p].displacement += sqrt(pow(dxDot,2)+pow(dyDot,2)+pow(dzDot,2));

}





