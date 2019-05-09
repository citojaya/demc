#include "common.h"


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
    //delete inactive particles from the domain
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



