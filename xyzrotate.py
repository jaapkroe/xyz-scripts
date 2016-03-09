#!/usr/bin/python2

import sys, string, argparse
from math import sqrt
import numpy as np

parser = argparse.ArgumentParser(description="rotate xyz file")#, usage=usage)
parser.add_argument("-m", "--method", type=int, default=2, help="method to define rotation matrix 1 (euler angles) or 2 (vector and angle) (default=2)")
parser.add_argument("-t", "--theta", type=float, default=0, help="rotation angle corresponding to vector (in degrees)")
parser.add_argument('-u', '--vector', type=float, nargs=3, default = [1, 0, 0],  help="rotation vector (e.g. -u '1 0 0')")
parser.add_argument('-s', '--scale', type=float, nargs=3, default = [1, 1, 1],  help="scaling vector (e.g. -u '1 1 1')")
parser.add_argument('-e', '--euler', type=float, nargs=3, default = [90, 0, 0], help="euler angles  (e.g. -e '90 45 0')")
parser.add_argument('-o', '--out', type=argparse.FileType('wb', 0), default="out.xyz", help="output file name")
parser.add_argument('-c', '--center', action="store_false", help="center molecule")
parser.add_argument('-v', '--verbose', action='count', default=0, help="be verbose")
parser.add_argument('xyzfile', type=file, help='xyz-file')

def main():
    args = parser.parse_args()
    cell, coords, types, cm = loadxyz(args.xyzfile, args.center)
    # rotate and print
    if args.method==1: 
      if args.verbose: print 'rotating over euler angles', np.array(args.euler)
      rotmat = rotation1(np.array(args.euler))
    elif args.method==2:
      if args.verbose: print 'rotating %s degrees over vector %s' %( args.theta, np.array(args.vector) )
      rotmat = rotation2(args.theta, np.array(args.vector))
    else: 
      raise ValueError, 'method %s does not exist' % args.method
    scalemat = np.diag ( np.array(args.scale) )
    if args.verbose: print 'rotation matrix :\n', rotmat, 'determinant = ',np.linalg.det(rotmat)
    if args.verbose: print 'scale matrix :\n', scalemat, 'determinant = ',np.linalg.det(scalemat)
    T = np.dot(scalemat,rotmat)
    if args.verbose: print 'transformation matrix :\n', T, 'determinant = ',np.linalg.det(T)
    cell=np.diagonal(T*cell)
    for i, c in enumerate(coords):
      coords[i]=np.dot(T,c) 
      if args.verbose>1: print "%s %12.8f %12.8f %12.8f -> %12.8f %12.8f %12.8f" % (types[i], c[0], c[1], c[2], coords[i][0], coords[i][1], coords[i][2])
      for j in range(3):
          coords[i][j] += cm[j]
    writexyz(args.out, cell, coords, types)

def loadxyz(file, center):
    # header1 (number of atoms)
    size = int(file.readline().strip())
    # header2 (cell size)
    items = string.split(file.readline())
    try: cell = [float(x) for x in items]
    except: cell = [0,0,0]
    if len(cell)!=3: cell = [0,0,0]
    # content (atom coordinates)
    coords = []; types = []
    for i in xrange(size):
        items = string.split(file.readline())
        coords.append([float(x) for x in items[1:4]])
        types.append(items[0])
    if center:
      # shift to cm
      cm = [ sum(x)/float(len(coords)) for x in zip(*coords) ]
      for c in coords: 
        for i in xrange(3): c[i] -= cm[i]
    else: cm = [0,0,0]
    return cell, coords, types, cm

def writexyz(file, cell, coords, types, mode='w'):
    file.write(str(len(coords))+'\n')
    file.write("%12.8f %12.8f %12.8f\n" % (cell[0], cell[1], cell[2]))
    for i, l in enumerate(coords):
      file.write("%s %12.8f %12.8f %12.8f\n" % (types[i], l[0], l[1], l[2]))  # xyz format

def rotation1(theta, R = np.zeros((3,3))):
    # function to return rotation matrix from three angles in degrees
    theta = np.array(theta)*np.pi/180.
    cx, cy, cz = np.cos(theta)
    sx, sy, sz = np.sin(theta)
    R.flat = (cx*cz, -cy*sz+sy*sx*cz,  sy*sz+cy*sx*cy, 
              cx*sz,  cy*cz+sy*sx*sz, -sy*cz+cy*sx*sy, 
              -sx  ,  sy*cx         ,  cy*cx         )
    return R

def rotation2(theta, u, R = np.zeros((3,3))):
    # function to return rotation matrix from one angles and a vector
    theta *= np.pi/180
    c=np.cos(theta); s=np.sin(theta)
    ux, uy, uz = u/sqrt(u[0]**2+u[1]**2+u[2]**2) # normalize
    R.flat = (c+ux**2*(1-c)   , ux*uy*(1-c)-uz*s, ux*uz*(1-c)+uy*s,
              uy*ux*(1-c)+uz*s, c+uy**2*(1-c)   , uy*uz*(1-c)-ux*s,
              uz*ux*(1-c)-uy*s, uz*uy*(1-c)+ux*s, c+uz**2*(1-c)   )
    return R

if __name__ == "__main__":
    sys.exit(main())
