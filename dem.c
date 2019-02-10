
#include "common.h"
#include "surf.h"
#include "mem.h"

#define C2 100.0

/* Initialize dem information
Executes only once at the begining of FLUEMT*/
DEFINE_INIT(my_init, d)
{ 
  xmax = 0.0, ymax=0.0, zmax=0.0;
  xmin = 0.0, ymin=0.0, zmin=0.0;
  largestParDia = 0.0;
 
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

  xmin = 0.042;
  xmax = 0.057;
  ymin = -0.001;
  ymax = 0.001;
  zmin = -0.007;
  zmax = 0.014;
  xmin = xmin - 4.0*largestParDia*conversion;
  xmax = xmax + 4.0*largestParDia*conversion;
  ymin = ymin - 4.0*largestParDia*conversion;
  ymax = ymax + 4.0*largestParDia*conversion;
  zmin = zmin - 4.0*largestParDia*conversion;
  zmax = zmax + 4.0*largestParDia*conversion;
  xDiv = floor((xmax-xmin)/(largestParDia*conversion*multif3));
  yDiv = floor((ymax-ymin)/(largestParDia*conversion*multif3));
  zDiv = floor((zmax-zmin)/(largestParDia*conversion*multif3));

  domainDx = (xmax-xmin)/xDiv;
  domainDy = (ymax-ymin)/yDiv;
  domainDz = (zmax-zmin)/zDiv;

  demInit();
  fflush(stdout);
}

/*
Initialize DPM particle information
Executes only once on injected particles
*/
DEFINE_DPM_INJECTION_INIT(solid_paritcles, I)
{
  Injection *I2;
  Injection *Ilist = Get_dpm_injections();
  //Get_dmp_injections();
  np = 0;
  loop(I2,Ilist)
  {
    Particle *p;
    loop(p,I2->p_init)
    {
      //Set particle mass according to DEM density input
      P_MASS(p) = (4.0/3.0)*PI*pow((0.5*P_DIAM(p)),3.0)*dens;
      np++;
    }
  }

  //Setup DEM scaling 
  setReduceUnits();



  //Insert particles to cell 
  loop(I2,Ilist)
  {
    Particle *p;
    loop(p,I2->p_init)
    {
      //Insert to cell
      int iIndex = ceil((P_POS(p)[0]*lengthFactor-xmin)/domainDx);
      int jIndex = ceil((P_POS(p)[1]*lengthFactor-ymin)/domainDy);
      int kIndex = ceil((P_POS(p)[2]*lengthFactor-zmin)/domainDz);
      int cellIndex = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;   
      insertToBdBox(p->part_id, cellIndex);//add particle to the cell

      // writeLogLine("logfile2.log","**************\n");
      // writeLogNum("logfile2.log","BD BOX  PARTICLE ",p->part_id);
      // writeLogNum("logfile2.log","CELL INDEX ",cellIndex);
      //writeLogNum("logfile2.log","BD BOX  PARTICLE Y ",P_POS(p)[1]);
      //writeLogNum("logfile2.log","BD BOX  PARTICLE Z ",P_POS(p)[2]);
      //writeLogLine("logfile2.log","**************\n");
      demPart[p->part_id].pos[0] =  P_POS(p)[0]*lengthFactor;
      demPart[p->part_id].pos[1] =  P_POS(p)[1]*lengthFactor;
      demPart[p->part_id].pos[2] =  P_POS(p)[2]*lengthFactor;
      demPart[p->part_id].dia = P_DIAM(p)*lengthFactor;
      demPart[p->part_id].mass = P_MASS(p)*massFactor;
      demPart[p->part_id].inert = 2.0*demPart[p->part_id].mass*pow(0.5*demPart[p->part_id].dia,2)/5.0; 
      demPart[p->part_id].vel[0] = 0.0;
      demPart[p->part_id].vel[1] = 0.0;
      demPart[p->part_id].vel[2] = 0.0;

      if(P_DIAM(p)*lengthFactor > largestParDia){
        largestParDia = P_DIAM(p)*lengthFactor;
      }
    }
  }

  //Update neighbour list
  loop(I2,Ilist)
  {
    Particle *p;
    loop(p,I2->p_init)
    {
      updateNeighbourList(p->part_id);
    }
  }

  writeLogNum("logfile2.log","XMIN ",xmin/lengthFactor);
  writeLogNum("logfile2.log","YMIN ",ymin/lengthFactor);
  writeLogNum("logfile2.log","ZMIN ",zmin/lengthFactor);
  writeLogLine("logfile2.log","----------------\n");
  writeLogNum("logfile2.log","XMAX ",xmax/lengthFactor);
  writeLogNum("logfile2.log","YMAX ",ymax/lengthFactor);
  writeLogNum("logfile2.log","ZMAX ",zmax/lengthFactor);
  writeLogLine("logfile2.log","----------------\n");
  writeLogNum("logfile2.log","X LENGTH ",(xmax-xmin)/lengthFactor);
  writeLogNum("logfile2.log","Y LENGTH ",(ymax-ymin)/lengthFactor);
  writeLogNum("logfile2.log","Z LENGTH ",(zmax-zmin)/lengthFactor);
  writeLogLine("logfile2.log","----------------\n");
  writeLogNum("logfile2.log","DOMAIN DX ",domainDx/lengthFactor);
  writeLogNum("logfile2.log","DOMAIN DY ",domainDy/lengthFactor);
  writeLogNum("logfile2.log","DOMAIN DZ ",domainDz/lengthFactor);
  writeLogLine("logfile2.log","----------------\n");
  writeLogNum("logfile2.log","DOMAIN XDIV ",xDiv);
  writeLogNum("logfile2.log","DOMAIN YDIV ",yDiv);
  writeLogNum("logfile2.log","DOMAIN ZDIV ",zDiv);
  writeLogLine("logfile2.log","----------------\n");

  //fflush(stdout);
}

// DEFINE_DPM_SCALAR_UPDATE(scalar_update, c, t, initialize, p)
// {
//   int ip=p->part_id;
//   demPart[ip].force[0] = 0.0;
//   demPart[ip].force[1] = -demPart[ip].mass; //gravitational force
//   demPart[ip].force[2] = 0.0;


//   demPart[ip].pos[0] =  P_POS(p)[0]*lengthFactor;
//   demPart[ip].pos[1] =  P_POS(p)[1]*lengthFactor;
//   demPart[ip].pos[2] =  P_POS(p)[2]*lengthFactor;

//   demPart[ip].vel[0] = P_VEL(p)[0]*velocityFactor;
//   demPart[ip].vel[1] = P_VEL(p)[1]*velocityFactor;
//   demPart[ip].vel[2] = P_VEL(p)[2]*velocityFactor;

//   real gap  =lengthFactor*(P_POS(p)[1] - P_DIAM(p)*0.5);
//   uVec[0] = 0.0;
//   uVec[1] = 1.0;
//   uVec[2] = 0.0;
 
//   if(gap < 0) //If contact exists calculate contact force
//   {
//           //if(fabs(gap) > P_DIAM(p)*lengthFactor*overlapLimit){gap = -P_DIAM(p)*lengthFactor*overlapLimit;}
//     surfaceContactForce(ip, -gap, uVec);
//   }
//   moveDPM(p);

//   fflush(stdout);
// }

/*
Executes at every DPM time step
This macro is used to call DEM functions which calculate particle-particle contact forces 
*/
//DEFINE_DPM_SOURCE(particle_source, c, t, S, strength, p)
// //

// DEFINE_DPM_SCALAR_UPDATE(scalar_update, c, t, initialize, p)
// {
//   writeLog("logfile2.log","DEFINE_DPM_SCALAR_UPDATE ",(real)p->part_id);
//   demLoop(p);// Update DEM current time
//   move(p); //move to the new position

//   //find cell index
//   int iIndex = ceil((P_POS(p)[0]*lengthFactor-xmin)/domainDx);
//   int jIndex = ceil((P_POS(p)[1]*lengthFactor-ymin)/domainDy);
//   int kIndex = ceil((P_POS(p)[2]*lengthFactor-zmin)/domainDz);
//   int cellIndex = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;
//   deleteParticle(p, cellIndex);//delete particle from the existing cell 
//   insertToBdBox(p, cellIndex);//add particle to the cell
// }


// DEFINE_DPM_BODY_FORCE(particle_body_force,p,i)
// {
//      writeLog("logfile2.log","PART ",p->part_id);
//      writeLog("logfile2.log","BODY FORCE ",0.0;
// }
/*
When particles are in contact with surface this macro is executed by FLUENT
Loop through all particles and update contact surface normal vector and contact surface node positions
Only applicable when contact surface is triangular
*/
DEFINE_DPM_BC(bc_reflect, p, t, f, f_normal, dim)
{
  long int ip = p->part_id;

  //fflush(stdout);
  //return PATH_END;
  return PATH_ABORT;
  //return PATH_ACTIVE;
}

// DEFINE_CONTACT(contact_props, dt, contacts)
// {
//   writeLogLine("logfile2.log", "DEFINE_CONTACT\n");
//   Objp *o;
//   face_t face;
//   Thread *thread;
//   loop(o, contacts)
//   {
//     face = O_F(o);
//     thread = O_F_THREAD(o);
//     writeLogNum("logfile2.log", "THREAD ID", THREAD_ID(thread));
//   }

// }

DEFINE_DPM_EROSION(wall_erosion, p, t, f, normal, alpha, Vmag, mdot)
{
  writeLogNum("logfile2.log", "PART ",p->part_id);
  writeLogNum("logfile2.log", "THREAD ID", THREAD_ID(t));
}

DEFINE_EXECUTE_AT_END(execute_at_end)
{
  time_count++;
 
  int iter = CURRENT_TIMESTEP/(timeStep/timeFactor);
  Injection *I;
  Injection *Ilist = Get_dpm_injections();

  // loop(I,Ilist)
  // {
  //   Particle *p;
  //   loop(p,I->p)
  //   { 
  //     demPart[p->part_id].vel[0] = P_VEL(p)[0]*velocityFactor;
  //     demPart[p->part_id].vel[1] = P_VEL(p)[1]*velocityFactor;
  //     demPart[p->part_id].vel[2] = P_VEL(p)[2]*velocityFactor;
  //     // demPart[p->part_id].pos[0] = P_POS(p)[0]*velocityFactor;
  //     // demPart[p->part_id].pos[1] = P_POS(p)[1]*velocityFactor;
  //     // demPart[p->part_id].pos[2] = P_POS(p)[2]*velocityFactor;
  //   }
  // }
  for(int i=0; i<iter; i++)
  {
    loop(I,Ilist)
    {
      Particle *p;
      loop(p,I->p)
      { 
        forceCalculation(p);
      }
    }

    loop(I,Ilist)
    {
      Particle *p;
      loop(p,I->p)
      { 
        updatePosition(p);
      }
    }
    
    //Reset boundary box
    for(int i=0; i<xDiv*yDiv*zDiv; i++)
    {
      bdBox[i].noOfParticles = 0;
    }

    loop(I,Ilist)
    {
      Particle *p;
      loop(p,I->p)
      {
        //Update DPM velocity
        if(updateDPM == 0){
          P_POS(p)[0] = demPart[p->part_id].pos[0]/lengthFactor;
          P_POS(p)[1] = demPart[p->part_id].pos[1]/lengthFactor; 
          P_POS(p)[2] = demPart[p->part_id].pos[2]/lengthFactor;
        }
        //Insert to cell
        // int iIndex = ceil((P_POS(p)[0]*lengthFactor-xmin)/domainDx);
        // int jIndex = ceil((P_POS(p)[1]*lengthFactor-ymin)/domainDy);
        // int kIndex = ceil((P_POS(p)[2]*lengthFactor-zmin)/domainDz);
        int iIndex = ceil((demPart[p->part_id].pos[0]-xmin)/domainDx);
        int jIndex = ceil((demPart[p->part_id].pos[1]-ymin)/domainDy);
        int kIndex = ceil((demPart[p->part_id].pos[2]-zmin)/domainDz);
        int cellIndex = iIndex + jIndex*xDiv + kIndex*xDiv*yDiv;     
        insertToBdBox(p->part_id, cellIndex);//add particle to the cell
      }
    }

    loop(I,Ilist)
    {
      Particle *p;
      loop(p,I->p)
      {
        if(demPart[p->part_id].displacement > allowedDisp)
        {
          updateNeighbourList(p->part_id);
        }
      }
    }
  }//end if iteration

  if(time_count > 2)
  {
    demSave();
    time_count = 0;
  }
}

/*DEFINE_EXECUTE_AT_END(execute_at_end)
{
  int i;
  int n;
  Node *node;
 
  Thread *tf;
  Domain *d = Get_Domain(1);
  Thread *t = Lookup_Thread(d,1);
  face_t f;
  cell_t c;

  //thread_loop_c(t,d)
///  {
    Injection *I;
    Injection *Ilist = Get_dpm_injections();

    loop(I,Ilist)
    {
      Particle *p;
      loop(p,I->p)
      {
        int ip = p->part_id;
        demPart[ip].force[0] = 0.0;
        demPart[ip].force[1] = 0.0;//-demPart[ip].mass; //gravitational force
        demPart[ip].force[2] = 0.0;


        // demPart[ip].pos[0] =  P_POS(p)[0]*lengthFactor;
        // demPart[ip].pos[1] =  P_POS(p)[1]*lengthFactor;
        // demPart[ip].pos[2] =  P_POS(p)[2]*lengthFactor;

        // demPart[ip].vel[0] = P_VEL(p)[0]*velocityFactor;
        // demPart[ip].vel[1] = P_VEL(p)[1]*velocityFactor;
        // demPart[ip].vel[2] = P_VEL(p)[2]*velocityFactor;


        cell_t c = P_CELL(p);
        c_face_loop(c,t,i)//Loop over faces
        {
          int t_id = THREAD_ID(C_FACE_THREAD(c,t,i));
          //for(int n=0; n<noOfWalls; n++)
          //{
          // if(t_id == walls[n])
           /* if(t_id == 5)
            {
              f = C_FACE(c,t,i);
              tf = C_FACE_THREAD(c,t,i);

              //printf("FACES -------------%d\n",i);
              //int face_t_id = THREAD_ID(C_FACE_THREAD(cc,t,i));
              int j=0;
              f_node_loop(f,tf,n)//Loop over face nodes
              {
                node = F_NODE(f,tf,n);
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
            }

              real *uVec = allocateDoubleArray(DIM);
              getUnitVector(demPart[ip].faceNode1,demPart[ip].faceNode2,demPart[ip].faceNode3,uVec);
              real *parPos = allocateDoubleArray(DIM);
              parPos[0] = P_POS(p)[0]*lengthFactor;
              parPos[1] = P_POS(p)[1]*lengthFactor;
              parPos[2] = P_POS(p)[2]*lengthFactor;
              //real gap = getOverlap(parPos,P_DIAM(p)*lengthFactor, demPart[ip].faceNode1,  
                //                    demPart[ip].faceNode2, demPart[ip].faceNode3,  uVec);

              real gap  =lengthFactor*(P_POS(p)[1] - P_DIAM(p)*0.5);
              uVec[0] = 0.0;
              uVec[1] = 1.0;
              uVec[2] = 0.0;
              writeLog("logfile2.log","PART ",p->part_id);
              writeLog("logfile2.log","PAR TIME ",P_TIME(p));
              if(gap < 0) //If contact exists calculate contact force
              {
                //if(fabs(gap) > P_DIAM(p)*lengthFactor*overlapLimit){gap = -P_DIAM(p)*lengthFactor*overlapLimit;}
                //surfaceContactForce(p, -gap, uVec);
              }
              free(uVec);
              free(parPos);
            }
          }
        }
       //writeLog("logfile2.log","PART ",
       // moveDPM(p);
      }//end of particle loop
    }
 // }

  //cell_t c = P_CELL(p);
 

  //writeLog("logfile2.log","EXECUTE AT END\n",CURRENT_TIME);
  demSave();

}*/

DEFINE_ADJUST(define_adjust, d)
{
    // if(updateDPM == 1)
    // {
    //writeLog("logfile2.log","ADJUST ",0.0);=[;\]


      Injection *I;
      Injection *Ilist = Get_dpm_injections();

      loop(I,Ilist)
      {
        Particle *p;
        loop(p,I->p)
        { 
          //writeLog("logfile2.log","DPM PARTICLE ID ",p->part_id);
          //writeLog("logfile2.log","DPM TIME ",P_TIME(p));/
          
          // P_POS(p)[0] = demPart[p->part_id].pos[0]/lengthFactor;
          // P_POS(p)[1] = demPart[p->part_id].pos[1]/lengthFactor;
          // P_POS(p)[2] = demPart[p->part_id].pos[2]/lengthFactor;

          // P_VEL(p)[0] = 0.0;//demPart[p->part_id].vel[0]/velocityFactor;
          // P_VEL(p)[1] = 0.0;//demPart[p->part_id].vel[1]/velocityFactor;
          // P_VEL(p)[2] = 0.0;//demPart[p->part_id].vel[2]/velocityFactor;
        
        }
      }
    //}
    updateDPM = 0;
  
  //writeLog("logfile2.log","DEFINE_ADJUST\n",CURRENT_TIME);
  //demLoop();
}

DEFINE_ON_DEMAND(demand)
{
  writeInjectionFile("initial.inj");
}

DEFINE_EXECUTE_AT_EXIT(execute_at_exit)
{
  // Delete dynamic memeory
  //free(dpmList);
  fclose(LogFile);
  free(walls);
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
