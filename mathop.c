#include "common.h"

double getCenterDist(int ip, int jp){
    return sqrt((parPosX[ip]-parPosX[jp])*(parPosX[ip]-parPosX[jp])+
        (parPosY[ip]-parPosY[jp])*(parPosY[ip]-parPosY[jp])+
        (parPosZ[ip]-parPosZ[jp])*(parPosZ[ip]-parPosZ[jp]));
 }

