#!/usr/bin/env python

#import vtk
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


p.add_option("-v", action="store_true", dest="verbose",  help="Verbose")
#p.add_option("-a", action="store_true", dest="ascii",  help="ASCII, instead of binary, VTK output")
(opts, args) = p.parse_args()

if len(args) < 1:
   print ("check input values")
   sys.exit(0)
(dia) = args[0]



#mesh_files = glob.glob("*.mesh")
#padit = int(math.log10(len(mesh_files))+1)

# space = 1.0
# xDiv = math.floor(float(xmax)/(float(dia)+space))
# yDiv = math.floor(float(ymax)/(float(dia)+space))
# zDiv = math.floor(float(zmax)/(float(dia)+space))

# print("BEFORE xDiv, yDiv, zDiv ",xDiv,yDiv,zDiv)


minX = minY = minZ = 100000
maxX = maxY = maxZ = -100000
# zDiv -= 25

# print("AFTER xDiv, yDiv, zDiv ",xDiv,yDiv,zDiv)

parX = []
parY = []
parZ = []

if opts.verbose:print ("Generating random particles\n")
xDiv = 3
yDiv = 40
zDiv = 25


xDiv = 10
yDiv = 2
zDiv = 10


yShift = -0.4
zShift = -4.8
xShift = 46

no_of_part = 0

gap = float(dia)*0.5
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
f = open("initial.inj","w")
#f.write(str(len(vertices_data))+"\n")

f.write(str(no_of_part)+"\n");
for k in range(zDiv):
   for j in range(yDiv):
        for i in range(xDiv):
            sx = str(round((parX[i]+xShift)*1e-3,6))
            sy = str(round((parY[j]+yShift)*1e-3,6))
            sz = str(round((parZ[k]+zShift)*1e-3,6))
            sd = str(round(float(dia)*1e-3,6))
            # s = " ".join(str(round(parX[i]*1e-3,4)))
            # if(parZ[k]+zShift > 25 and parZ[k]+zShift < 45):
            f.write(sx+" "+sy+" "+sz+" "+sd+"\n")
            # no_of_part  += 1
            # f.write("(("+str(round(parX[i]*1e-3,4))+" "+str(parY[j])+" "+str(parZ[k])+" 0.0 0.0 0.0 "+str(float(dia)*1e-3)+" 0.0 1.0))\n")
            

f.close()
# Write veritces 


# count = 0
# path = "particle.dat"
# if opts.verbose:print (" "+path)
# fin = open(path,"r")
# for line in fin:
#   tuple = line.split()
#   if(count >0):
#     if(tuple[0] == "TIME"):
#       f.write('\nZONE T = particle\n')
#     else:
#       f.write(" ".join(tuple)+" 0.0\n")
#     # f.write(tuple[0]+" "+tuple[1]+" "+tuple[2]+" "+tuple[3]+" "+tuple[4]+" "+tuple[5]+" "+tuple[6]+" 0.0\n")
#   count += 1

# fin.close()

# f.close()

# print('written to "',file_name[:-5]+'.dat"')
print ("DONE")
sys.exit(3)
