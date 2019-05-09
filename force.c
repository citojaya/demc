#include "common.h"

/* Calculate particle-particle, particle-wall contact forces and drag force due to fluid*/
void forceCalculation(int p)
{
    double dt = timeStep;
    demPart[p].force[0] = 0.0; //gravitational force
    demPart[p].force[1] = 0.0;//gravitational force
    demPart[p].force[2] = -demPart[p].mass; 
    demPart[p].momentum[0] = 0.0;
    demPart[p].momentum[1] = 0.0;
    demPart[p].momentum[2] = 0.0;

    //findContactFromDPM(p);
    //findContactFromMesh(p);
    findContactFromBoundary(p);
 
    //Find particle-particle contact force
    neighbourContactForce(p);

    //Find drag force on particles from fluid 
    calculateDragForce(p);
}

/*Particle-particle vanderwal force*/
void ppVWForce(int ip, int jp, double vGap){
    
    double ipDia = demPart[ip].dia;
    double jpDia = demPart[jp].dia;
    double vGapMn = 1.0e-9*lengthFactor;
    
    double ijHa = sqrt(demPart[ip].haPp*demPart[jp].haPp);
    if(vGap < vGapMn){vGap = vGapMn;}
    
    double rStar = ipDia*jpDia/(ipDia+jpDia);
    double fv = -ijHa*rStar/(6.*vGap*vGap);

    if(ip == 1){
        writeLogNum("logfile3.log","MG ",demPart[ip].mass*9.81);
        writeLog3Num("logfile3.log","COHE-FORCE ",vGap*1e6/lengthFactor, 1e6*fv/forceFactor, 1e6*demPart[ip].mass/forceFactor);
    }
    double *uVec = allocateDoubleArray(DIM);
    vecSub(demPart[ip].pos, demPart[jp].pos, uVec);
    unitVec(uVec, uVec);

    demPart[ip].force[0] += uVec[0]*fv;
    demPart[ip].force[1] += uVec[1]*fv;
    demPart[ip].force[2] += uVec[2]*fv;
    free(uVec);
}

/*Particle-wall vanderwal force*/
void pWallVWForce(int p, double vGap, double *uVec){
    double vGapMn = 1.0e-9*lengthFactor;
    if(vGap < vGapMn){vGap = vGapMn;}

    double fv = -demPart[p].haPp*demPart[p].dia*0.5/(6.*vGap*vGap);

    demPart[p].force[0] += uVec[0]*fv;
    demPart[p].force[1] += uVec[1]*fv;
    demPart[p].force[2] += uVec[2]*fv;   
}

/*Particle-particle electrostatic force*/
void ppElectForce(int ip, int jp, double gap){
    double gap1 = gap/lengthFactor;
    double *uVec = allocateDoubleArray(DIM);
    vecSub(demPart[ip].pos, demPart[jp].pos, uVec);
    unitVec(uVec, uVec);

    double fe = -demPart[ip].eCharge*demPart[jp].eCharge*forceFactor/(4.*PI*permitivity*gap1*gap1);
    // if(ip == 1)
    //    writeLogNum("logfile3.log","E-FORCE ",1e6*fe/forceFactor);
    demPart[ip].force[0] += uVec[0]*fe;
    demPart[ip].force[1] += uVec[1]*fe;
    demPart[ip].force[2] += uVec[2]*fe;
    free(uVec);        
}

/*Partcile-particle capillary force*/
/*void ppCapillaryForce(int ip, int jp, double gap){
    
    double s_rupture = (1+0.5*cont_ang)*pow(liq_vol,1/3);
    if(gap/lengthFactor < s_rupture){
        double *uVec = allocateDoubleArray(DIM);
        vecSub(demPart[ip].pos, demPart[jp].pos, uVec);
        unitVec(uVec, uVec);
        double sepMin = 5.0e-6; //((Hornbaker, Albert et al. 1997, Nase, Vargas et al. 2001)
        double separation = fmax(sepMin, gap/lengthFactor);
        double rStar = (demPart[ip].dia*demPart[jp].dia/(demPart[ip].dia+demPart[jp].dia))/lengthFactor;
        double capf = -2.*PI*rStar*surf_tens*cos(cont_ang)/(1+1/sqrt(1+2.*liq_vol/(PI*rStar*pow(separation,2))-separation));
        // if(ip == 1)
        //    writeLogNum("logfile3.log","E-FORCE ",1e6*fe/forceFactor);
        demPart[ip].force[0] += uVec[0]*capf*forceFactor;
        demPart[ip].force[1] += uVec[1]*capf*forceFactor;
        demPart[ip].force[2] += uVec[2]*capf*forceFactor;
        free(uVec);  
    }         
}*/

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

    //projVec(demPart[p].hisDisp, uVec, disp, 1);
    //vecAdd(disp, tipCntPntDisp, disp);
    sclVecMult(1.0,tipCntPntDisp,disp);
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
        if(dd != 0.0){
            fdt = fdt*(1.0 - pow((1.0 - dd/dsmax),1.5));
            sclVecMult(-fdt/dd,disp,fdtVec);
        }
        else{
            sclMult(0.0,fdtVec);
        }
    }

    //sclVecMult(1.0,disp, demPart[p].hisDisp);

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

    free(rotMom);
    free(fdtVec);
    free(tipCntPntDisp);
    free(cntPntDisp);
    free(totalForce);
    free(momentum);
    free(disp);
}
/*
void calculateDragForce(int p){
    double x[ND_ND];
    Thread *tc = P_CELL_THREAD(p);
    cell_t c = P_CELL(p);
    //Fluid drag force calculation
    double velFX = C_U(c,tc)*velocityFactor;  
    double velFY = C_V(c,tc)*velocityFactor;  
    double velFZ = C_W(c,tc)*velocityFactor;    
    double pGX = C_P_G(c,tc)[0]*pressureFactor/lengthFactor;
    double pGY = C_P_G(c,tc)[1]*pressureFactor/lengthFactor;
    double pGZ = C_P_G(c,tc)[2]*pressureFactor/lengthFactor;
    double velPX = demPart[p].vel[0];
    double velPY = demPart[p].vel[1];
    double velPZ = demPart[p].vel[2];
    double density = C_R(c,tc)*densityFactor;
    double visc = C_MU_L(c,tc)*massFactor/(lengthFactor*timeFactor);

    C_CENTROID(x,c,tc);
    int cellIndex = C_UDMI(c,tc,0);
   
    double instPor = cfdcell[cellIndex].porosity;
    if(instPor > 0.9){instPor = 0.9;}
    
    double relVelMag = sqrt((velFX-velPX)*(velFX-velPX)+(velFY-velPY)*(velFY-velPY)+(velFZ-velPZ)*(velFZ-velPZ));
    double Re = instPor*demPart[p].dia*relVelMag*density/visc;
  
    double beta, dCoeff;

    // //-- Standard equations --//
	if (Re < 1000){
		dCoeff = 24.*(1.+0.15*pow(Re, 0.687))/Re;
	}
	else {
		dCoeff = 0.44;
	}

    // -  Second method  Yi He paper-//

    if(instPor <= 0.8){
        beta = 150.*(1.-instPor)*(1.-instPor)*visc/(instPor*pow(demPart[p].dia,2)) 
                + 1.75*(1.-instPor)*density*relVelMag/demPart[p].dia;
    }
    else if(instPor > 0.8){
        beta = (3.0/4.0)*dCoeff*relVelMag*density*(1.0-instPor)*pow(instPor,-2.7)/demPart[p].dia;
    }

    double pfFX = 0.0;
    double pfFY = 0.0;
    double pfFZ = 0.0;
    
    pfFX = (velFX-velPX)*partVol(p)*beta/(1.0-instPor); //always instPor != 0
    pfFY = (velFY-velPY)*partVol(p)*beta/(1.0-instPor);
    pfFZ = (velFZ-velPZ)*partVol(p)*beta/(1.0-instPor);

    double pGFX = -pGX*PI*pow(demPart[p].dia,3)/6.0;
    double pGFY = -pGY*PI*pow(demPart[p].dia,3)/6.0;
    double pGFZ = -pGZ*PI*pow(demPart[p].dia,3)/6.0;

    //Update force on particles
    demPart[p].force[0]  += pfFX + pGFX;
    demPart[p].force[1]  += pfFY + pGFY; 
    demPart[p].force[2]  += pfFZ + pGFZ + partVol(p)*density;

    //Store drag force on particle for later use in calculating source term for fluid
    demPart[p].dragFX = pfFX + pGFX;
    demPart[p].dragFY = pfFY + pGFY;
    demPart[p].dragFZ = pfFZ + pGFZ + partVol(p)*density;

    if(updateSourceTerm == 1){
        cfdcell[cellIndex].solidVol += partVol(p)/volumeFactor;
        cfdcell[cellIndex].noOfParts += 1;
        cfdcell[cellIndex].dragFX +=  -demPart[p].dragFX/forceFactor; 
        cfdcell[cellIndex].dragFY +=  -demPart[p].dragFY/forceFactor; 
        cfdcell[cellIndex].dragFZ +=  -demPart[p].dragFZ/forceFactor;
    }

}*/


