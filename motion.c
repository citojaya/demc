#include "common.h"

/*Update particle position and insert them into cells*/
void updatePosition(int p){
    demPart[p].insertable = 1;
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
    // if(demPart[p].pos[0] > XMAX_BOUND*lengthFactor){
    //     demPart[p].active = 0;
    // }

    //Update DPM particle position
    //if(demPart[p].active == 1){
        //Insert to cell
        //ceil gives upper bound, but we need lower bound. Therefore use -1 for index
        int iIndex = ceil((demPart[p].pos[0]-xmin)/domainDx);
        int jIndex = ceil((demPart[p].pos[1]-ymin)/domainDy);
        int kIndex = ceil((demPart[p].pos[2]-zmin)/domainDz);
        int cellIndex = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;  
        //if(demPart[p].active == 1){
       
    //}
    // Thread *tc = P_CELL_THREAD(p);
    // cell_t c = P_CELL(p);
    // C_UDMI(c,tc,0) = 0.0;  
}



