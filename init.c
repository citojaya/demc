#include "common.h"

//allocated = 0;
/*Read particle information and domain information
*/
void demInit(){
    prevCPUTime = 0;
    cycleCount = 0;
    saveDEM = 0;
    //updateDPM= 0;
    time_count = 0;// time counter used for saving output file
    particle_counter = 0; //counter for keeping the track of number of particles
    demTime = 0.0;

    readDomain("infile");

    xDiv = ceil((xmax-xmin)/(largestParDia*multif3));
    yDiv = ceil((ymax-ymin)/(largestParDia*multif3));
    zDiv = ceil((zmax-zmin)/(largestParDia*multif3));

    domainDx = (xmax-xmin)/xDiv;
    domainDy = (ymax-ymin)/yDiv;
    domainDz = (zmax-zmin)/zDiv;

    //writeLogLine("DOMAIN domainDx, domainDy, domainDz\n");
    printf("Domain dx,dy,dz %lf, %lf, %lf\n",domainDx, domainDy, domainDz);
    printf("xDiv, yDiv, zDiv %d, %d, %d \n",xDiv, yDiv, zDiv);

    if(LogFile){
        fclose(LogFile);
    }
    //LogFile = fopen("logfile.log", "a");
    if(sizeof(walls) != 0){
        //printf("SIZE OF WALLS %d\n",sizeof(walls));
        free(walls);
        free(uVec);
        free(ipRVec);
        free(jpRVec);
        free(bdBox);
        free(ijVec);
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
    } 


    //Read material data
    readInput("infile", &parArraySize, &dens, &ymod, &pois, &sfc, &rec, &dmpn, &rf, &cyldia, &timeStep, 
            &noOfWalls, &haConst);

    allocate();
    //Read particle-wall contact surfaces
    //readWalls("infile", walls);
    //printf("TIME STEP %lf\n",timeStep/timeFactor);
    readGeom("infile");

    //writeLogNum("logfile.log","Update DPM ",updateDPM);
    printf("Particle Array Size %d\n",parArraySize);
    writeLogNum("logfile3.log","Particle Array Size ",parArraySize);
    writeLog3Num("logfile3.log", "Initial min ",xmin,ymin,zmin);
    writeLog3Num("logfile3.log", "Initial max ",xmax,ymax,zmax);
    writeLogNum("logfile3.log","ductxmin ",ductxmin);
    writeLogNum("logfile3.log","ductxmax ",ductxmax);
    writeLogNum("logfile3.log","ductxedge1 ",ductxedge1);
    writeLogNum("logfile3.log","ductxedge2 ",ductxedge2);
    writeLogNum("logfile3.log","ductymin ",ductymin);
    writeLogNum("logfile3.log","ductymax ",ductymax);
    writeLogNum("logfile3.log","ductzmin ",ductzmin);
    writeLogNum("logfile3.log","ductzmax ",ductzmax);
    writeLogNum("logfile3.log","ductzedge ",ductzedge);
    writeLogNum("logfile3.log","maxVel ",maxVel);
    writeLogNum("logfile3.log","PP haConst ",haConst*1.0e20);
    writeLogNum("logfile3.log","lamda1 ",lamda1*1e9);
    writeLogNum("logfile3.log","lamda2 ",lamda2*1e9);
    writeLogNum("logfile3.log","rms1 ",rms1*1e9);
    writeLogNum("logfile3.log","rms2 ",rms2*1e9);
    
    writeLogNum("logfile3.log","Particle Array Size ",parArraySize);
    writeLogNum("logfile3.log","density ",dens);
    writeLogNum("logfile3.log","Youngs Modulus ",ymod);
    writeLogNum("logfile3.log","Timestep ",timeStep);
}

/* Allocate arrays */
void allocate(){
    uVec = allocateDoubleArray(DIM);
    ipRVec = allocateDoubleArray(DIM);
    jpRVec = allocateDoubleArray(DIM);
    ijVec = allocateDoubleArray(DIM);
    rotVel = allocateDoubleArray(DIM);
    ipCntPntVel = allocateDoubleArray(DIM);
    jpCntPntVel = allocateDoubleArray(DIM);
    cntPntVel = allocateDoubleArray(DIM);

    //parNb = allocateIntArray(np*nbSize); //neighbourlist of each particle
    //parNoOfNb = allocateIntArray(np); //no of neighbours in each particle
    bdBox = allocateBdBoxArray(xDiv*yDiv*zDiv); //bounding box array
    
    for(int i=0; i<xDiv*yDiv*zDiv; i++){
        //bdBox[i].noOfFaces = 0;
        bdBox[i].noOfParticles = 0;
        //bdBox[i].noOfFluidCells = 0;
        bdBox[i].fluidVelX = 0.0;
        bdBox[i].fluidVelY = 0.0;
        bdBox[i].fluidVelZ = 0.0;
        bdBox[i].fluidVolF = 0.0;
        bdBox[i].dragFX = 0.0;
        bdBox[i].dragFY = 0.0;
        bdBox[i].dragFZ = 0.0;
    }
    
    demPart = allocatePar(parArraySize);
    for(int i=0; i<parArraySize; i++){
        demPart[i].dt = 0.0;
        demPart[i].currentTime = 0.0;
        demPart[i].nrmDisp = 0.0;
        demPart[i].noOfNeigh = 0;
        demPart[i].pos = allocateDoubleArray(DIM);
        demPart[i].angVel = allocateDoubleArray(DIM);
        demPart[i].vel = allocateDoubleArray(DIM);
        demPart[i].hisDisp = allocateDoubleArray(DIM);

        demPart[i].hisDisp[0]  = 0.0;
        demPart[i].hisDisp[1]  = 0.0;
        demPart[i].hisDisp[2]  = 0.0;   

        demPart[i].angVel[0] = 0.0;
        demPart[i].angVel[1] = 0.0;
        demPart[i].angVel[2] = 0.0;

        demPart[i].force = allocateDoubleArray(DIM);
        demPart[i].momentum = allocateDoubleArray(DIM);
        demPart[i].momentum[0] = 0.0;
        demPart[i].momentum[1] = 0.0;
        demPart[i].momentum[2] = 0.0;
        demPart[i].faceNode1 = allocateDoubleArray(DIM);
        demPart[i].faceNode2 = allocateDoubleArray(DIM);
        demPart[i].faceNode3 = allocateDoubleArray(DIM);
        demPart[i].surfNorm = allocateDoubleArray(DIM);
        demPart[i].displacement = 0.0;
        demPart[i].insertable = 1;
        //demPart[i].noOfCntF = 0;
        
    }


    //walls = allocateIntArray(noOfWalls);
    //printf();
    //dpmList = (Tracked_Particle)malloc(2*sizeof(Tracked_Particle));

}

/* Generate cells for the problem domain which is used by particle and fluent cells. 
These cells are used to reduce the number of searches when particle and fluid cells exchange 
information 
*/
void buildDEMCFDCellMap(){
    printf("buildDEMCFDCellMap\n");
}



/*
Insert particle to BdBox
param:
pI - particle index
cI - cell index
*/
/*void insertToBdBox(int p, int cI){
    bdBox[cI].parts[bdBox[cI].noOfParticles] = p;
    bdBox[cI].noOfParticles++;
    demPart[p].prevCellIndex = cI;
    if(bdBox[cI].noOfParticles > NO_OF_PARTICLES_IN_BDCELL){
        printf("bdBox[cI].noOfParticles > NO_OF_PARTICLES_IN_BDCELL\n");
    }
}*/

/*
Delete particle from bdBox
param:
pI - particle index
cI = cell index
*/
/*
void deleteParticle(int p, int cI){
    for(int i=0; i<bdBox[cI].noOfParticles; i++){
        int np = bdBox[cI].parts[i];
        if(p == np){
            bdBox[cI].parts[i] = bdBox[cI].parts[bdBox[cI].noOfParticles-1];
            bdBox[cI].noOfParticles--;
            break;
        }
    }
}*/


/*
Initially neighbourlist array is empty. Fill neighbourlist. 
*/
void initialize(double *nbList, int *parIndex, int *cellSE, int np,
    double *pos, double *parDia)
{
 
}

/*
Assign scale factors
*/
void setReduceUnits()
{
    //Scale factors for reduced units
	refLength = largestParDia;
	refDensity = largestParDensity;
	lengthFactor = 1.0/refLength;
	volumeFactor = pow(lengthFactor,3);
	massFactor = 6.0/(PI*pow(refLength,3)*refDensity);
	timeFactor = sqrt(gravity/refLength);
	densityFactor = 6.0/(PI*refDensity);
	forceFactor = 6.0/(gravity*PI*pow(refLength,3)*refDensity);
	pressureFactor = 6.0/(gravity*PI*refLength*refDensity);
	StressFactor = pressureFactor;
	energyFactor = 6.0/(gravity*PI*pow(refLength,4)*refDensity);
	momentFactor = energyFactor;
	powerFactor = 6.0/(pow(gravity,1.5)*PI*pow((double)refLength,3.5)*refDensity);
	velocityFactor = 1.0/sqrt(refLength*gravity);
	accFactor = 1.0/gravity;
	angVelFactor = sqrt(refLength/gravity);
	angAccFactor = refLength/gravity;
	freqFactor = sqrt(refLength/gravity);
	inertiaFactor = 6.0/(PI*pow(refLength,5)*refDensity);

    cutGap = 1.2*largestParDia*conversion*lengthFactor;

    dsmaxCff = sfc*(2.0-pois)/(2.0*(1.0-pois));
    //writeLogNum("logfile3.log"," DS MAX",dsmaxCff);
    dti = 0.0;
    dd = 0.0;
    dsmax = 0.0;

    // scale particle properties
    ymod = ymod*pressureFactor;
    cyldia = cyldia*lengthFactor;
	timeStep = timeStep*timeFactor;
	maxTime	= maxTime*timeFactor;

    elasticMod = ymod/(1.0-pow(pois,2));

    double haa = 0.;
    for (int i=0; i<np; i++){
        demPart[i].pos[0] = demPart[i].pos[0]*lengthFactor;
        demPart[i].pos[1] = demPart[i].pos[1]*lengthFactor;
        demPart[i].pos[2] = demPart[i].pos[2]*lengthFactor;
        demPart[i].dt = timeStep;
        demPart[i].dia = demPart[i].dia*lengthFactor;
        demPart[i].mass = (4.0/3.0)*PI*pow((0.5*demPart[i].dia),3.0)*dens*densityFactor;
        demPart[i].inert = 2.0*demPart[i].mass*pow(0.5*demPart[i].dia,2)/5.0;
        demPart[i].displacement = 2.0*allowedDisp;
        haa = 6.5E-20;
        //demPart[i].ha = haa*forceFactor*lengthFactor;       
    }


    //Find allowed displacement for neighbourlist update
    rIn = largestParDia*conversion*lengthFactor;
    rOut = 1.55*rIn; //By definition
    allowedDisp = 0.5*(rOut-rIn);

    //Adjust boundary limits
    xmin = xmin*lengthFactor;
    ymin = ymin*lengthFactor;
    zmin = zmin*lengthFactor;
    xmax = xmax*lengthFactor;
    ymax = ymax*lengthFactor;
    zmax = zmax*lengthFactor;

    domainDx = domainDx*lengthFactor;
    domainDy = domainDy*lengthFactor;
    domainDz = domainDz*lengthFactor;

    cellRadius = 0.5*sqrt(domainDx*domainDx+domainDy*domainDy);
    //writeLogNum("logfile3.log","CEll R ",cellRadius/lengthFactor);
    printf("Hmarker Constant  %lf\n ",haa);
}

/*
void allocateMat(struct MatType *m){
    printf("allocation.c\n");
    //struct MatType *mt = (struct MatType *)malloc(num_of_mats*sizeof(struct MatType));
    //m = mt;    
}*/

/*
Allocate an integer type array
return: int* 
*/
int *allocateIntArray(int size)
{
    int *val = (int*)malloc(size*sizeof(int));
    memset(val,0,size*sizeof(int));
    return val;
}

/*
Allocate a double* type array
return: double* 
*/
double *allocateDoubleArray(int size)
{
    double *val = (double*)malloc(size*sizeof(double));
    memset(val,0.0,size*sizeof(double));
    return val;
}

/*
Allocate a char* type array
return: char* 
*/
char *allocateCharArray(int size)
{
    char *val = (char*)malloc(size*sizeof(char));
    return val;
}

/*
Allocate bounding box type array
return: BdBox*
*/
struct BdBox *allocateBdBoxArray(int size)
{
    //printf("BD BOX SIZE %d\n",size);
    struct BdBox *bdB = (struct BdBox*)malloc(size*sizeof(struct BdBox));
    return  bdB;
}

/*
Allocate particle array
*/
struct demParticle *allocatePar(int np)
{
    struct demParticle *par = (struct demParticle*)malloc(np*sizeof(struct demParticle));
    return par;
}

