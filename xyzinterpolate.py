#!/usr/bin/python
import sys,string,math,operator
from optparse import OptionParser
import argparse

parser = argparse.ArgumentParser(description="make a (linear) interpolation between 2 xyz-files")

parser.add_argument("-n", dest="steps",  type=int, default=9, help="number of interpolation steps")
parser.add_argument("-e", dest="extra",  type=int, default=0, help="extend with extrapolation steps")
parser.add_argument("-o", dest="output", type=int, default=0, help="type of output (0=xyz, 1=gNEB, 2=LAMMPS)")
parser.add_argument("-c", dest="cellut", type=int, default=0, help="type of output (0=xyz, 1=gNEB, 2=LAMMPS)")
parser.add_argument("-p", dest="cell",   type=float,nargs=3,help='cell size  px py pz',default=[100,100,100])
parser.add_argument("-d", dest="shift",  type=float,nargs=3,help='cell shift dx dy dz',default=[0,0,0])
parser.add_argument('file1', type=argparse.FileType('r'),help='xyz-file 1')
parser.add_argument('file2', type=argparse.FileType('r'),help='xyz-file 2')

def main():
  args = parser.parse_args()

  coords1,types1=loadxyz(args.file1)
  coords2,types2=loadxyz(args.file2)
  N=len(coords1)
  steps=args.steps
  extra=args.extra

  # coordinate differences
  delta=[]; rmsd=0.0; atomdisp=[]
  for i in range(N):      # for all coordinate lines
      deltai = []
      atomdisp.append(0)
      for j in range(3):  # and all dimensions
          d = (coords2[i][j]-coords1[i][j])/(steps-1)
          deltai.append(d)   #find the difference
          rmsd += d*d
          atomdisp[i] += d*d
      atomdisp[i] = math.sqrt(atomdisp[i])
      delta.append(deltai)
  maxd = max(atomdisp); maxi=atomdisp.index(maxd)
  sys.stderr.write("%s -> %s in %d frames, displacement per step = %.3g / %.3g (%d) = total / max\n"%\
          (args.file1.name,args.file2.name,steps+extra,math.sqrt(rmsd),maxd,maxi))
  
  distance_check=1
  # now we want to create a number of intermediate steps (coordinate files)
  for s in range(steps+extra):
    if distance_check: 
      rlist = []
      for i in range(N):
        xi = coords1[i]
        for j in range(i+1,N):
          xj = coords1[j]
          dr = 0.
          for k in range(3): 
              dx = xi[k]-xj[k]
              if(args.cell): dx -= round(dx/args.cell[k])*args.cell[k]
              dr += dx**2
          rlist.append((dr,i,j))
      m = min(rlist, key=lambda x: x[0])
      sys.stderr.write("%d minimum distance = %.2f (i,j = %d,%d)\n"%(s,m[0],m[1],m[2]))
    if args.output==0: printxyz(coords1,types1,s)
    elif args.output==1: printgneb(coords1,types1,s+1)
    elif args.output==2: printdata(coords1,types1,s+1,args.cell,args.shift)
    else: sys.exit("output type unknown")
    for i in range(N): #change coordinates 
      for j in range(3): coords1[i][j] += delta[i][j]


def printxyz(coords,types,s):
  N = len(coords)
  print(str(N)+'\nframe '+str(s))
  for i in range(N): print("%3s %16.10f %16.10f %16.10f" % ( types[i], coords[i][0], coords[i][1], coords[i][2] ))
  
def printgneb(coords,types,s):
  a0 = 1. / .5291772108
  N = len(coords)
  print("Image : "+str(s))
  for i in range(N): print("%3s %16.10f %16.10f %16.10f %d" % ( types[i], coords[i][0]*a0, coords[i][1]*a0, coords[i][2]*a0, i))
  
def printdata(coords,types,s,cell,shift):
  N = len(coords)
  numerictypes = list(sorted(set(types)))
  atomids = [numerictypes.index(t)+1 for t in types]
  print("frame %s"%s)
  print("%d atoms"%N)
  print("%d atom types\n"%len(numerictypes))
  xlo=shift[0];xhi=cell[0]+shift[0]
  ylo=shift[1];yhi=cell[1]+shift[1]
  zlo=shift[2];zhi=cell[2]+shift[2]
  print("%g %g xlo xhi"%(xlo,xhi))
  print("%g %g ylo yhi"%(ylo,yhi))
  print("%g %g zlo zhi\n\nAtoms\n"%(zlo,zhi))
  for i in range(N): print("%3d %d %6.3f %16.10f %16.10f %16.10f" % ( i+1, atomids[i], 0, coords[i][0], coords[i][1], coords[i][2]))
  
def loadxyz(f):
    f.readline                          # header1 (number of atoms)
    size = int(f.readline().strip())
    f.readline()                        # skip header2 
  
    # read coordinates
    coords = []; types = []
    for i in range(size):
        items = f.readline().split()
        types.append(items[0])
        pos = [float(items[1]), float(items[2]), float(items[3])]
        coords.append(pos)
    return coords,types

if __name__ == "__main__":
    sys.exit(main())
