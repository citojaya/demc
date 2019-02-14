
#include "common.h"
#include "surf.h"
#include "mem.h"

#define C2 100.0

/* Initialize dem information
Executes only once at the begining of program*/
DEFINE_INIT(my_init, d)
{ 
  xmax = 0.0, ymax=0.0, zmax=0.0;
  xmin = 0.0, ymin=0.0, zmin=0.0;
 
  cell_t c;
  real x[ND_ND];
  Thread *t = Lookup_Thread(d,1);
  begin_c_loop(c,t)
    C_CENTROID(x,c,t);
    if(x[0] > xmax) xmax = x[0];
    if(x[1] > ymax) ymax = x[1];
    if(x[2] > zmax) zmax = x[2];
    if(x[0] < xmin) xmin = x[0];
    if(x[1] < ymin) ymin = x[1];
    if(x[2] < zmin) zmin = x[2];
  end_c_loop(c,t)

  xmin = xmin - offset*conversion;
  xmax = xmax + offset*conversion;
  ymin = ymin - offset*conversion;
  ymax = ymax + offset*conversion;
  zmin = zmin - offset*conversion;
  zmax = zmax + offset*conversion;
  int xD = floor((xmax-xmin)/(largestParDia*conversion*multif3));
  int yD = floor((ymax-ymin)/(largestParDia*conversion*multif3));
  int zD = floor((zmax-zmin)/(largestParDia*conversion*multif3));

  xDiv = xD;
  yDiv = yD;
  zDiv = zD;
  domainDx = (xmax-xmin)/xD;
  domainDy = (ymax-ymin)/yD;
  domainDz = (zmax-zmin)/zD;

  // printf("xmin, xmax %lf,%lf\n",xmin,xmax);
  // printf("ymin, ymax %lf,%lf\n",ymin,ymax);
  // printf("zmin, zmax %lf,%lf\n",zmin,zmax);

  // printf("dx,dy,dz %lf,%lf,%lf\n",domainDx,domainDy,domainDz);
  // printf("xD,yD,zD %d,%d,%d\n",xD,yD,zD);
  //Particle *p;
  demInit();
  fflush(stdout);
  //exit(0);
}

DEFINE_DPM_INJECTION_INIT(solid_paritcles, I)
{
  printf("TRACKED PARTICLE\n");
  Injection *I2;
  Injection *Ilist = Get_dpm_injections();
  //Get_dmp_injections();
  loop(I2,Ilist)
  {
    printf("INJECTIONS INITIAL\n");

    Particle *p;
    loop(p,I2->p_init)
    {
      //Set particle mass according to DEM density input
      demPart[np].cfdp = p;
      P_MASS(p) = (4.0/3.0)*PI*pow((0.5*P_DIAM(p)),3.0)*dens;
      printf("INITIAL PAR POS %lf,%lf,%lf\n",P_POS(p)[0],P_POS(p)[1],P_POS(p)[2]);
      np++;
    }
  }

  //Setup DEM scaling 
  //setScaleFactors(); 
  //Assign to bDBox cells
  addToBdBox();

  //Set DEM gravitational force
  //assignGravity();
  
  //Update neighbour list
  updateNeighbourList();
}

/*Execute at the begining of each iteration before conservation equations are solved
  Porosity and drag force information from DEM is transferred to CFD*/
DEFINE_ADJUST(my_adjust, d)
{
  printf("DEFINE_ADJUST **********************************************\n");
  cell_t cell;
  Thread *thread;
  //reset no of fluid cells in bounding box cells

  Injection *I;
  Injection *Ilist = Get_dpm_injections();
  
  // Update FLUENT particle postion and velocity 
  int ip = 0; 
  loop(I,Ilist)
  {
    Particle *p;
    loop(p,I->p)
    {
      // P_POS(p)[0] = demPart[ip].pos[0]/lengthFactor;
      // P_POS(p)[1] = demPart[ip].pos[1]/lengthFactor;
      // P_POS(p)[2] = demPart[ip].pos[2]/lengthFactor;
      // P_VEL(p)[0] = demPart[ip].vel[0]/velocityFactor;
      // P_VEL(p)[1] = demPart[ip].vel[1]/velocityFactor;
      // P_VEL(p)[2] = demPart[ip].vel[2]/velocityFactor;
      ip++;
    }
  }

  fflush(stdout);
}


/* Execute at end of each CFD timestep and update background cells which is required for DEM
    NOTE: Each timestep may have number of iterations
*/
DEFINE_EXECUTE_AT_END(execute_at_end)
{
  printf("DEFINE_EXECUTE_AT_END **********************************************\n");
  Domain *d;
  Thread *t;
  cell_t c;
  real x[ND_ND]; //Fluent real array for storing centroid coordinates

  d = Get_Domain(1);
  t = Lookup_Thread(d,1);

  Injection *I;
  Injection *Ilist = Get_dpm_injections();
  
  //Update DEM particle position and velocity
  int ip = 0; 
  loop(I,Ilist)
  {
    Particle *p;
    loop(p,I->p)
    {
      // demPart[ip].pos[0] = P_POS(p)[0]*lengthFactor;
      // demPart[ip].pos[1] = P_POS(p)[1]*lengthFactor;
      // demPart[ip].pos[2] = P_POS(p)[2]*lengthFactor;
      // demPart[ip].vel[0] = P_VEL(p)[0]*velocityFactor;
      // demPart[ip].vel[1] = P_VEL(p)[1]*velocityFactor;
      // demPart[ip].vel[2] = P_VEL(p)[2]*velocityFactor;
      ip++;
    }
  }

  //start DEM simulation
  //demLoop();
  //save DEM output data 
  //demSave();
  fflush(stdout);
}

DEFINE_DPM_BC(bc_reflect, p, t, f, f_normal, dim)
{
  printf("DEFINE_DPM_BC **********************************************\n");
  Injection *I;
  Injection *Ilist = Get_dpm_injections();
  
  //Update DEM particle position and velocity
  int ip = 0; 
  loop(I,Ilist)
  {
    Particle *p;
    int ip = 0;
    loop(p,I->p)//Loop over particles
    {
      //printf("F_NORMAL %lf,%lf,%lf\n",f_normal[0],f_normal[1],f_normal[2]);
      Node *node;
      int n;
      int i=0;
          
      demPart[ip].surfNorm[0] = -f_normal[0]; 
      demPart[ip].surfNorm[1] = -f_normal[1];
      demPart[ip].surfNorm[2] = -f_normal[2];
      f_node_loop(f,t,n)//Loop over face nodes
      {
        node = F_NODE(f,t,n);
        if(i==0){
          demPart[ip].faceNode1[0] = NODE_X(node);
          demPart[ip].faceNode1[1] = NODE_Y(node);
          demPart[ip].faceNode1[2] = NODE_Z(node);
        }
        if(i==1){
          demPart[ip].faceNode2[0] = NODE_X(node);
          demPart[ip].faceNode2[1] = NODE_Y(node);
          demPart[ip].faceNode2[2] = NODE_Z(node);
        }
        if(i==2){
          demPart[ip].faceNode3[0] = NODE_X(node);
          demPart[ip].faceNode3[1] = NODE_Y(node);
          demPart[ip].faceNode3[2] = NODE_Z(node);
        }
        i++;
      }
      //Calculate surface contact force
      boundaryContactForce(ip, demPart[ip].faceNode1, demPart[ip].faceNode2,demPart[ip].faceNode3, demPart[ip].surfNorm);
      ip++;
   
      //printf("NO OF NODES %d\n", F_NNODES(f,t)); 
    }
  }
  return PATH_ACTIVE;
}
// DEFINE_DPM_BC(bc_reflect,p,t,f,f_normal,dim)
// {
//   printf("DPM_BC\n");
//   return PATH_ABBORT;
// }

DEFINE_SOURCE(xmom_source,c,t,dS,eqn)
{
  real x[ND_ND];
  real con, source;
  C_CENTROID(x,c,t);
  con = C2*0.5*C_R(c,t)*x[1];
  source = -con*fabs(C_U(c, t))*C_U(c,t);
  dS[eqn] = -2.*con*fabs(C_U(c,t));
  return source;
} 

DEFINE_SOURCE(ymom_source,c,t,dS,eqn)
{
  real x[ND_ND];
  real con, source;
  C_CENTROID(x,c,t);
  con = C2*0.5*C_R(c,t)*x[1];
  source = -con*fabs(C_U(c, t))*C_U(c,t);
  dS[eqn] = -2.*con*fabs(C_U(c,t));
  return source;
} 

DEFINE_SOURCE(zmom_source,c,t,dS,eqn)
{
  real x[ND_ND];
  real con, source;
  C_CENTROID(x,c,t);
  con = C2*0.5*C_R(c,t)*x[1];
  source = -con*fabs(C_U(c, t))*C_U(c,t);
  dS[eqn] = -2.*con*fabs(C_U(c,t));
  return source;
} 

DEFINE_EXECUTE_AT_EXIT(execute_at_exit)
{
  // Delete dynamic memeory
  free(uVec);
  free(ipRVec);
  free(jpRVec);

  free(bdBox);
    
  free(ijVec);
    //free(tempVec);
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
  printf("EXECUTE AT EXIT\n");
  fflush(stdout);
}

/* Rotate the capsule for a given speed*/
// DEFINE_CG_MOTION(piston,dt,vel,omega,time,dtime)
//  {
//     omega[2]=50;
//  } 
