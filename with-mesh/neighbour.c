#include "common.h"

/*
Check whether contact calculation has been completed
*/
int contactNotDone(int ip, int jp){
    for(int i=0; i<demPart[ip].contacts; i++){
        if(demPart[ip].contDone[i] == jp){
           return 0;
        }
    }
    return 1;
}

/*
Check for existing neighbour
*/
int insertable(int ip, int jp){
    for(int i=0; i<demPart[ip].noOfNeigh; i++){
        if(demPart[ip].neigh[i] == jp){
            return 0;
        }
    }
    return 1;
}

/*
Add a new neighbour to the neighbour list
*/
void addNeighbour(int ip, int jp){
    demPart[ip].neigh[demPart[ip].noOfNeigh] = jp;
    demPart[ip].noOfNeigh++;
}

/*
Delete neighbour
*/
void deleteNeighbour(int ip, int jp){
    for (int i=0; i<demPart[ip].noOfNeigh; i++){
        int neigh = demPart[ip].neigh[i];
        if(neigh == jp){
            demPart[ip].neigh[i] = demPart[ip].neigh[demPart[ip].noOfNeigh-1];
            demPart[ip].noOfNeigh--;
            break;
        }
    }

}

/*
Insert particle to BdBox
param:
pI - particle index
cI - cell index
*/
void insertToBdBox(int p, int cI){

    bdBox[cI].parts[bdBox[cI].noOfParticles] = p;
    bdBox[cI].noOfParticles++;
    demPart[p].prevCellIndex = cI;
    if(bdBox[cI].noOfParticles > NO_OF_PARTICLES_IN_BDCELL){
        printf("bdBox[cI].noOfParticles > NO_OF_PARTICLES_IN_BDCELL\n");
    }
}

/*
Insert faces to DEM cells
Done only at the begining of the programm
*/
void addFacesAndEdges(){
    double max_no_of_faces = 0.0;
     //printf("wall min max %lf,%lf,%lf\n",wallMxmin/lengthFactor,wallMymin/lengthFactor,wallMzmin/lengthFactor);
     //printf("wall min max %lf,%lf,%lf\n",wallMxmax/lengthFactor,wallMymax/lengthFactor,wallMzmax/lengthFactor);
     //printf("wall min max %lf,%lf,%lf\n",wallMxmin/lengthFactor,wallMymin/lengthFactor,wallMzmin/lengthFactor);
     //printf("wall dx %lf,%lf,%lf\n",wallMdx/lengthFactor,wallMdy/lengthFactor,wallMdz/lengthFactor);

    /*****  Insert faces into wallMesh  *******/
    for(int i=0; i<noOfWallSurfaces; i++){
      double centX = face[i].centroid[0];
      double centY = face[i].centroid[1];
      double centZ = face[i].centroid[2];
      //printf("wall centroid %lf,%lf,%lf\n",centX/lengthFactor,centY/lengthFactor,centZ/lengthFactor);
     
      int iIndex = ceil((centX-wallMxmin)/wallMdx);
      int jIndex = ceil((centY-wallMymin)/wallMdy);
      int kIndex = ceil((centZ-wallMzmin)/wallMdz);
      int cI = iIndex + jIndex*wmxDiv + kIndex*wmxDiv*wmyDiv;

    //printf("NO of surfaces %d\n", noOfWallSurfaces);
      for(int r=kIndex-1; r<kIndex+2; r++){
         for(int q=jIndex-1; q<jIndex+2; q++){
             for(int p=iIndex-1; p<iIndex+2; p++){
                int cI = p + q*wmxDiv + r*wmxDiv*wmyDiv;
                double cX = p*wallMdx - wallMdx*0.5 + wallMxmin;
                double cY = q*wallMdy - wallMdy*0.5 + wallMymin;
                double cZ = r*wallMdz - wallMdz*0.5 + wallMzmin;
                
                double dist = sqrt((centX-cX)*(centX-cX)+(centY-cY)*(centY-cY)+(centZ-cZ)*(centZ-cZ));
                
                if(dist < (wallMdx+wallMdy+wallMdz)/3.0 && faceExist(cI, i) == 0){
                    //printf("DIST %lf\n", 1e3*dist/lengthFactor);
                    wallMeshBox[cI].wallFaceId[wallMeshBox[cI].noOfFaces] = i;
                    wallMeshBox[cI].noOfFaces = wallMeshBox[cI].noOfFaces + 1;
                    wallMeshBox[cI].isBoundary = 1;
                    max_no_of_faces = fmax(max_no_of_faces, wallMeshBox[cI].noOfFaces);
                    //printf("Max no of faces %d\n", wallMeshBox[cI].noOfFaces);
                }
             }
         }
       }
      
    }

    //printf("no of surfaces %d\n",noOfWallSurfaces);
    for(int i=0; i<wmxDiv*wmyDiv*wmzDiv; i++){
        if(wallMeshBox[i].noOfFaces > 0){
            printf("Max no of faces %d\n", wallMeshBox[i].noOfFaces);
        }
    }



    /*****  Insert edges into wallMesh  *******/
    for(int i=0; i<0.5*noOfVertices; i++){
      double centX = edge[i].center[0];
      double centY = edge[i].center[1];
      double centZ = edge[i].center[2];
      printf("edge center %lf,%lf,%lf\n",centX/lengthFactor,centY/lengthFactor,centZ/lengthFactor);
     
      int iIndex = ceil((centX-wallMxmin)/wallMdx);
      int jIndex = ceil((centY-wallMymin)/wallMdy);
      int kIndex = ceil((centZ-wallMzmin)/wallMdz);
      //int cI = iIndex + jIndex*wmxDiv + kIndex*wmxDiv*wmyDiv;

    //printf("NO of surfaces %d\n", noOfWallSurfaces);
      for(int r=kIndex-1; r<kIndex+2; r++){
         for(int q=jIndex-1; q<jIndex+2; q++){
             for(int p=iIndex-1; p<iIndex+2; p++){
                int cI = p + q*wmxDiv + r*wmxDiv*wmyDiv;
                double cX = p*wallMdx - wallMdx*0.5 + wallMxmin;
                double cY = q*wallMdy - wallMdy*0.5 + wallMymin;
                double cZ = r*wallMdz - wallMdz*0.5 + wallMzmin;
                
                double dist = sqrt((centX-cX)*(centX-cX)+(centY-cY)*(centY-cY)+(centZ-cZ)*(centZ-cZ));
                
                if(dist < (wallMdx+wallMdy+wallMdz)/3.0 && edgeExist(cI, i) == 0){
                    //printf("DIST %lf\n", 1e3*dist/lengthFactor);
                    wallMeshBox[cI].edgeId[wallMeshBox[cI].noOfEdges] = i;
                    wallMeshBox[cI].noOfEdges = wallMeshBox[cI].noOfEdges + 1;
                    wallMeshBox[cI].isEdgeBoundary = 1;
                    max_no_of_faces = fmax(max_no_of_faces, wallMeshBox[cI].noOfEdges);
                    //printf("Max no of faces %d\n", wallMeshBox[cI].noOfFaces);
                }
             }
         }
       }
    }
    //printf("no of surfaces %d\n",noOfWallSurfaces);
    for(int i=0; i<wmxDiv*wmyDiv*wmzDiv; i++){
        if(wallMeshBox[i].noOfEdges > 0){
            printf("Max no of edges %d\n", wallMeshBox[i].noOfEdges);
        }
    }
}

/*
Check for existing face
*/
int faceExist(int cI, int fId){
    for (int i=0; i<wallMeshBox[cI].noOfFaces; i++){
        if(wallMeshBox[cI].wallFaceId[i] == fId){
            return 1;
        }
    }
    return 0;
}
/*
Check for existing edge
*/
int edgeExist(int cI, int eId){
    for (int i=0; i<wallMeshBox[cI].noOfEdges; i++){
        if(wallMeshBox[cI].edgeId[i] == eId){
            return 1;
        }
    }
    return 0;
}

/*
Delete particle from bdBox
param:
pI - particle index
cI = cell index
*/
void deleteParticle(int p, int cI){
    for(int i=0; i<bdBox[cI].noOfParticles; i++){
        int np = bdBox[cI].parts[i];
        if(p == np){
            bdBox[cI].parts[i] = bdBox[cI].parts[bdBox[cI].noOfParticles-1];
            bdBox[cI].noOfParticles -= 1;
            break;
        }
    }
}

/*
Update neighbour list for a given particle
For a given particle scan through all neighbour cells and fetch particles within
the neighbour region
param:
pI - particle index
*/
void updateNeighbourList(int ip){
    if(demPart[ip].displacement > allowedDisp){
        //delete particle from previous cell
        deleteParticle(ip, demPart[ip].prevCellIndex);
        int iIndex = ceil((demPart[ip].pos[0]-xmin)/domainDx);
        int jIndex = ceil((demPart[ip].pos[1]-ymin)/domainDy);
        int kIndex = ceil((demPart[ip].pos[2]-zmin)/domainDz);
        int cI = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;
        //insert to new cell
        insertToBdBox(ip, cI);
        demPart[ip].displacement = 0.0;

        for(int i=0; i<demPart[ip].noOfNeigh; i++){
            int jp = demPart[ip].neigh[i];
            double xyz = getCenterDist(ip,jp);
            if (xyz > rOut){
                deleteNeighbour(jp,ip); //delete ip from jp
            }
        }

        //Update new neighbourlist for ip
        demPart[ip].noOfNeigh = 0;
        iIndex = ceil((demPart[ip].pos[0]-xmin)/domainDx);
        jIndex = ceil((demPart[ip].pos[1]-ymin)/domainDy);
        kIndex = ceil((demPart[ip].pos[2]-zmin)/domainDz);

        if(iIndex*jIndex*kIndex > xDiv*yDiv*zDiv){
            writeLogNum("logfile2.log","iIndex*jIndex*kIndex > xDiv*yDiv*zDiv ",iIndex*jIndex*kIndex);
        }
    
        for(int r=kIndex-1; r<kIndex+2; r++){
            for(int q=jIndex-1; q<jIndex+2; q++){
                for(int p=iIndex-1; p<iIndex+2; p++){
                    int neighCellIndex = p + q*xDiv + r*xDiv*yDiv;
                    if(neighCellIndex>xDiv*yDiv*zDiv){
                        writeLogNum("logfile2.log","neighCellIndex>xDiv*yDiv*zDiv ",neighCellIndex);
                    }
                    if(iIndex*jIndex*kIndex < 0){
                        writeLogNum("logfile2.log","iIndex*jIndex*kIndex < 0 ",iIndex*jIndex*kIndex);
                    }

                    for(int j=0; j<bdBox[neighCellIndex].noOfParticles; j++){
                        int jp = bdBox[neighCellIndex].parts[j];
                        if(getCenterDist(ip,jp) < rOut && ip != jp){
                            addNeighbour(ip,jp);
                        }
                    }
                }
            }
        }
    }//end of demPart[ip].displacement > allowedDisp
}


