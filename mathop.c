#include "common.h"

/* Return particle volume*/
double partVol(int p){
    return (4.0/3.0)*PI*pow(demPart[p].dia*0.5,3);
}

/*Find solid fraction within cell radius*/
double solidFraction(int ip){
    //Find cell center
    int iIndex = ceil((demPart[ip].pos[0] - xmin)/domainDx);
    int jIndex = ceil((demPart[ip].pos[1] - ymin)/domainDy);
    int kIndex = ceil((demPart[ip].pos[2] - zmin)/domainDz);
    int cellIndex = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;     
    
    double cX = iIndex*domainDx + 0.5*domainDx;
    double cY = jIndex*domainDy + 0.5*domainDy;
    double cZ = kIndex*domainDz + 0.5*domainDz;
    double vol = 0.0;
    double totVol = 0.0;
    for(int i=0; i<bdBox[cellIndex].noOfParticles; i++){
        int jp = bdBox[cellIndex].parts[i];
        double r = 0.5*demPart[jp].dia;
        double jpX = demPart[jp].pos[0];
        double jpY = demPart[jp].pos[1];
        double jpZ = demPart[jp].pos[2];
        
        double rR = sqrt((jpX-cX)*(jpX-cX) + (jpY-cY)*(jpY-cY) + (jpZ-cZ)*(jpZ-cZ)); 
        
        if( rR <= (0.5*domainDx-r)){
            vol = 0.8*partVol(jp);
        }
        else if(rR > (cellRadius-r)){
            vol = partVol(jp);
        }
        else if(rR <= (cellRadius-r) && rR  > (domainDx*0.5-r)){
            vol = 0.5*partVol(jp); 
        }
        totVol += vol;       
    }
    //writeLogNum("logfile3.log","VOL ",totVol/(lengthFactor*lengthFactor*lengthFactor));
    return totVol/((4.0/3.0)*PI*pow(cellRadius,3));
}

/* Return center distance of two particles
param:
p1 - particle 1
p2 - particle 2 */
double getCenterDist(int ip, int jp){
    double ipX = demPart[ip].pos[0];
    double ipY = demPart[ip].pos[1];
    double ipZ = demPart[ip].pos[2];
    double jpX = demPart[jp].pos[0];
    double jpY = demPart[jp].pos[1];
    double jpZ = demPart[jp].pos[2];

    double val = sqrt((ipX-jpX)*(ipX-jpX) + (ipY-jpY)*(ipY-jpY) + (ipZ-jpZ)*(ipZ-jpZ));
    return val;
 }

/* Add two vecotrs
param:
v1 - input vector 1
v2 - input vector 2
vec - resultant vector */
void vecAdd(double *v1, double *v2, double *vec){
   vec[0] = v1[0]+v2[0];
   vec[1] = v1[1]+v2[1];
   vec[2] = v1[2]+v2[2];
}

/* Substract two vecotrs
param:
v1 - input vector 1
v2 - input vector 2
vec - resultant vector */
void vecSub(double *v1, double *v2, double *vec){
   vec[0] = v1[0]-v2[0];
   vec[1] = v1[1]-v2[1];
   vec[2] = v1[2]-v2[2];
}

/* Find vector cross product
param:
v1 - input vector 1
v2 - input vector 2
vec - resultant vector */
void crossProd(double *v1, double *v2, double *vec){
    double temp1 = v1[1]*v2[2];
    double temp2 = v2[1]*v1[2];
    vec[0] = temp1 - temp2;

    temp1 = v2[0]*v1[2];
    temp2 = v1[0]*v2[2];
    vec[1] = temp1 - temp2;

    temp1 = v1[0]*v2[1];
    temp2 = v2[0]*v1[1];
    vec[2] = temp1 - temp2;   
}

/*
Dot product of two vectors
param:
v1 - vector 1
v2 - vector 2
*/
double dotProduct(double *v1, double *v2){
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

/*
Find normal unit vector to the surface defined by vector v1v2 and v2v3
*/
void getUnitVector(double *v1, double *v2, double *v3, double *uVec){
    double *v1v2 = allocateDoubleArray(DIM);
    double *v1v3 = allocateDoubleArray(DIM);

    vecSub(v1, v2, v1v2);
    vecSub(v1, v3, v1v3);
    crossProd(v1v2, v1v3, uVec);
    unitVec(uVec, uVec);

    free(v1v2);
    free(v1v3);
}

/* Find unit vector
param: 
v - input vector
vec - unit vector */
void unitVec(double *v, double *vec){
    if(vecMag(v) > 0){
        double temp1 = v[0]/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        double temp2 = v[1]/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        double temp3 = v[2]/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        vec[0] = temp1;
        vec[1] = temp2;
        vec[2] = temp3;
    }
    else{
        vec[0] = 0.0;
        vec[1] = 0.0;
        vec[2] = 0.0;
    }
}

/* Find magnitude of a vector*/
double vecMag(double *vec){
    return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
}

/* Multiply input vector by a scaler
param: 
scl - scaler
vec - vector to be multiplied*/
void sclMult(double scl, double *vec){
    vec[0] = scl*vec[0];
    vec[1] = scl*vec[1];
    vec[2] = scl*vec[2];
}

/* Multiply input vector by a scaler and assign to a vector
param: 
scl - scaler
inVec - input vector to be multiplied
outVec - reusltant output vector*/
void sclVecMult(double scl, double *inVec, double *outVec){
    outVec[0] = scl*inVec[0];
    outVec[1] = scl*inVec[1];
    outVec[2] = scl*inVec[2];
}


/* Returns the project vector on the plane defined by normal unit vector
param:
v - input vector
n - unit vector
vec - resultant project vector
type - either 0 or 1, 0 gives project vector, 1 gives project vector scaled by input vector */
void projVec(double *v, double *n, double *vec, int type){
    double *tV1 = malloc(DIM*sizeof(double));
    double *tV2 = malloc(DIM*sizeof(double));
    if(type == 0){
        crossProd(n,v, tV1);
        tV2[0] = -n[0];
        tV2[1] = -n[1];
        tV2[2] = -n[2];
        //printf("tV1 %lf,%lf,%lf\n ",tV1[0],tV1[1],tV1[2]);
        //printf("tV2 %lf,%lf,%lf\n ",tV2[0],tV2[1],tV2[2]);
        crossProd(tV2, tV1, vec);
        //printf("Type 0\n");
    }
    else{
        double *tV3 = malloc(DIM*sizeof(double));
        crossProd(n, v, tV1);
        tV2[0] = -n[0];
        tV2[1] = -n[1];
        tV2[2] = -n[2];
        crossProd(tV2, tV1, tV3);
        double temp;
        if((v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) != 0.0){
            temp = sqrt((tV3[0]*tV3[0]+tV3[1]*tV3[1]+tV3[2]*tV3[2])/
                            (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));  
        }
        else{
            temp = 0.0;
        }
        vec[0] = tV3[0]*temp;
        vec[1] = tV3[1]*temp;
        vec[2] = tV3[2]*temp;      
        free(tV3);
    }
    free(tV1);
    free(tV2);
}



