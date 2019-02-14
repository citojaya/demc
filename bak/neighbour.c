#include "common.h"



/*
Add a new neighbour to the neighbour list
*/
void addNeighbour(int ip, int jp){
    if(demPart[ip].noOfNeigh < NBSIZE){
        demPart[ip].neigh[demPart[ip].noOfNeigh] = jp;
        demPart[ip].noOfNeigh++;
    }
    else{
        //writeLogLine("demPart[ip].noOfNeigh > NBSIZE");
        printf("demPart[ip].noOfNeigh > NBSIZE");
    }
}

/*
Delete neighbour
*/
void deleteNeighbour(int ip, int jp){
    for (int i=0; i<demPart[ip].noOfNeigh; i++){
        int neigh = demPart[ip].neigh[i];
        if(neigh == jp){
            demPart[ip].neigh[i] = demPart[ip].neigh[demPart[i].noOfNeigh-1];
            demPart[ip].noOfNeigh--;
            break;
        }
    }

}

/*
Add particles to bounding box
*/
void addToBdBox(){
    //printf("NP %d\n", np);
    // //Reset to zero
    for(int i=0; i<xDiv*yDiv*zDiv; i++){
        bdBox[i].noOfParticles = 0;
    }

    // Injection *I;
    // Injection *Ilist = Get_dpm_injections();
    // //Reset neighbourlsit to zero
     for(int i=0; i<np; i++)
     {
    //     Particle *p;
    //     loop(p,I->p)
    //     {
    // //Assign to bDBox cells and update neighbour list
    
        int iIndex = ceil((demPart[i].pos[0]-xmin)/domainDx);
        int jIndex = ceil((demPart[i].pos[1]-ymin)/domainDy);
        int kIndex = ceil((demPart[i].pos[2]-zmin)/domainDz);
        int cellIndex = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;
        
        //printf("ADD TO CELL INDEX %d\n",cellIndex);
    //     //writeLog("logfile2.log","ADD TO BD BOX",(double)cellIndex);
    //     //insert particle to cell

         insertToBdBox(i,cellIndex);
    //     }
     }
    //     //printf("Cell index %d, %d, %d\n",iIndex,jIndex,kIndex);
    
}


/*
Update neighbour list for a given particle
For a given particle scan through all neighbour cells and fetch particles within
the neighbour region
param:
pI - particle index
*/
void updateNeighbourList(int ip){
    demPart[ip].noOfNeigh = 0;
    int iIndex = ceil((demPart[ip].pos[0]-xmin)/domainDx);
    int jIndex = ceil((demPart[ip].pos[1]-ymin)/domainDy);
    int kIndex = ceil((demPart[ip].pos[2]-zmin)/domainDz);

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
                    //printf("NEIGH CENT DIST %d %d %lf\n",ip,jp, getCenterDist(ip,jp)*1e3/lengthFactor);
                    if(getCenterDist(ip,jp)-0.5*(demPart[ip].dia+demPart[jp].dia) < cutGap && ip != jp){
                        //writeLogNum("logfile2.log","NEIGH ADDED ",jp);
                        //printf("NEIGH ADDED %d %d\n", ip, jp);
                        addNeighbour(ip,jp);
                    }
                }
            }
        }
    }
}



