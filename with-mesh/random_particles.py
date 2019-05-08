#!/usr/bin/env python

#import vtk
import random
import os
import sys
import re
import glob
import math
import numpy as np
from optparse import OptionParser, OptionValueError

arrSize = 1500

# parse command line
p = OptionParser(usage="""usage: %prog [options] <casename>

Generate random particles within a given domain

python random_particles.py <diameter>

""")

def orderedList():
    minX = minY = minZ = 100000
    maxX = maxY = maxZ = -100000

    parX = []
    parY = []
    parZ = []

    if opts.verbose:print ("Generating random particles\n")

    xDiv = 35 #max 49
    yDiv = 6  #max 8
    zDiv = 20 #35

    xShift = 102
    yShift = -4
    zShift = -47 #47

    no_of_part = 0

    gap = float(dia)*0.12
    print(xDiv,yDiv,zDiv)
    for k in range(zDiv):
        z = k*float(dia)+gap
        if(z < minZ):
            minZ = z
        if(z > maxZ):
            maxZ = z
        parZ.append(z)
        for j in range(yDiv):
            y = j*float(dia)+gap
            if(y < minY):
                minY = y
            if(y > maxY):
                maxY = y
            parY.append(y)
            for i in range(xDiv):
                x = i*float(dia)+gap
                if(x < minX):
                    minX = x
                if(x > maxX):
                    maxX = x
                parX.append(x)
                no_of_part += 1   
    print("MIN MAX ",minX+xShift,maxX+xShift,minY+yShift,maxY+yShift,minZ+zShift,maxZ+zShift)           
    print("Number of particles generated ",no_of_part,"\n")

    if opts.verbose:print ("Writing to injection file initial.inj")
    f = open("random.inj","w")
    #f.write(str(len(vertices_data))+"\n")

    for k in range(zDiv):
        for j in range(yDiv):
            for i in range(xDiv):
                sx = str(round((parX[i]+xShift)*1e-3,6))
                sy = str(round((parY[j]+yShift)*1e-3,6))
                sz = str(round((parZ[k]+zShift)*1e-3,6))
                sd = str(round(float(dia)*1e-3,4))
                # s = " ".join(str(round(parX[i]*1e-3,4)))
                # if(parZ[k]+zShift > 25 and parZ[k]+zShift < 45):
                f.write(sx+" "+sy+" "+sz+" "+sd+"\n")
                # no_of_part  += 1
                # f.write("(("+str(round(parX[i]*1e-3,4))+" "+str(parY[j])+" "+str(parZ[k])+" 0.0 0.0 0.0 "+str(float(dia)*1e-3)+" 0.0 1.0))\n")
                

    f.close()

def distance(px,py,pz,px2,py2,pz2):
    px = float(px)
    py = float(py)
    pz = float(pz)
    px2 = float(px2)
    py2 = float(py2)
    pz2 = float(pz2)

    dist = np.sqrt(np.power((px-px2),2)+np.power((py-py2),2)+np.power((pz-pz2),2))
    return dist

def insertable(px,py,pz,dia,particle):
    for x in range(len(particle)):
        (px2,py2,pz2) = particle[x]
        if(distance(px,py,pz,px2,py2,pz2) < float(dia)*1.02): 
            return False
        # if(float(px) == float(px2) and float(py) == float(py2) and float(pz) == float(pz2)):
        #     return False
        
    
    return True
    


def randomList(parts,xmin,xmax,ymin,ymax,zmin,zmax,centDist,dia):
    xMin = 1000
    yMin = 1000
    zMin = 1000
    xMax = -1000
    yMax = -1000
    zMax = -1000

    particle = []
    noOfParts = 0
    xdiv = (int)((xmax-xmin)/centDist)
    ydiv = (int)((ymax-ymin)/centDist)
    zdiv = (int)((zmax-zmin)/centDist)
    print("xdiv, ydiv, zdiv",xdiv,ydiv,zdiv)

    if opts.verbose:print ("Writing to injection file initial.inj")
    f = open("initial.inj","w")

    f.write(parts+"\n")
    while(noOfParts < int(parts)):
        px = random.randint(xmin,xmax)*0.1
        py = random.randint(ymin,ymax)*0.1
        pz = random.randint(zmin,zmax)*0.1
        if(insertable(px,py,pz,dia,particle)):
            sx = str(round((px)*1e-3,6))
            sy = str(round((py)*1e-3,6))
            sz = str(round((pz)*1e-3,6))
            size = str(round(0.1*float(dia)*1e-3,6))
            xMin = min(xMin,px)
            yMin = min(yMin,py)
            zMin = min(zMin,pz)
            xMax = max(xMax,px)
            yMax = max(yMax,py)
            zMax = max(zMax,pz)

            f.write(sx+" "+sy+" "+sz+" "+size+"\n")
            particle.append((px,py,pz))
            noOfParts += 1
            if(noOfParts%100 == 0):
                print("No of particles",noOfParts)

    # print("No of particles",noOfParts)
    print("Min Max",xMin,xMax,yMin,yMax,zMin,zMax)
    f.close()
    # print(particle)

def copyInjection(infile, outfile):
    zshift = 0.05
    f1 = open(infile, "r")
    f2 = open(outfile, "w")
    for line in f1:
        tuple = line.split()
        xx = tuple[0]
        yy = tuple[1]
        zz = float(tuple[2])
        zz += zshift
        zz = str(round(zz,6))

        f2.write(line)
        f2.write(xx+" "+yy+" "+zz+" 0.0 0.0 0.0 "+tuple[6]+" 0.0 1.0))\n")
        #print (xx[2:],yy,str(zz))

    f1.close()
    f2.close()

p.add_option("-v", action="store_true", dest="verbose",  help="Verbose")
#p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, VTK output")
(opts, args) = p.parse_args()

if len(args) < 2:
   print ("check input values")
   sys.exit(0)
# tuple = args.split()
(dia,parts) = args

# orderedList()

centerDist = float(dia)+0.02*float(dia)
xmin = 310
xmax = 390
ymin = -30
ymax = 30
zmin = -40
zmax = 0
randomList(parts,xmin,xmax,ymin,ymax,zmin,zmax,centerDist,dia)
#copyInjection("initial.inj", "combined.inj")

print ("DONE")
sys.exit(3)
