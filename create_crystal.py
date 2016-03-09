#!/usr/bin/env python
from numpy import *
import argparse,sys

examples = "EXAMPLES:\n\
create_crystal.py -s1   -a 1.50       -n 4 4 4  # simple cubic,\n\
create_crystal.py -s2   -a 2.20       -n 3 3 3  # BCC,\n\
create_crystal.py -s3   -a 2.20       -n 4 4 4  # FCC / rock-salt,\n\
create_crystal.py -s4   -a 1.42 3.30  -n 4 4 2  # graphite,\n\
create_crystal.py -s4   -a 1.42       -n 5 5 1  # graphene,\n\
create_crystal.py -s5   -a 3.46       -n 2 2 2  # diamond / zinc-blende,\n\
create_crystal.py -s6   -a 2.50 4.00  -n 4 4 4  # wurtzite"
parser = argparse.ArgumentParser(description="create a crystal structure file",epilog=examples,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-s", dest="structure",type=int, default=1, help="choose structure (1=sc, 2=bcc, 3=fcc, 4=graphite, 5=diamond, 6=wurtzite)")
parser.add_argument('-n', metavar='i', dest="ncells", type=int, nargs=3, default = [1, 1, 1],  help="number of unit cells in x, y and z direction (e.g. -n 3 2 4)")
parser.add_argument('-a', metavar='x', dest="lattice_constants", type=float, nargs='+', default = [1], help="supply the necessary lattice constants")
parser.add_argument('--scaled', action='store_true', help="use scaled coordinates (orthorhombic cells only)")
parser.add_argument('-t', dest="atomtype",type=str, default="C", help="primary atom type")
parser.add_argument('-d', dest="diatomic",type=str, help="second atom type (creates diatomic structure)")
parser.add_argument('-b', dest="basename", default="base", help="basename for output file")
parser.add_argument('-r', dest="orthorhombic", action='store_false', help="add atoms to create an orthorhomic cell")
parser.add_argument("-o", dest="output",type=int, default=1, help="type of output (1=xyz,2=lammps,3=for,4=(x)bs),5=obj")

def main():
    args      = parser.parse_args()
    lc        = args.lattice_constants
    basename  = args.basename

    if args.structure is 1:
      print("SIMPLE CUBIC"); a=lc[0]
      cell,coords = sc(a)
    elif args.structure is 2:
      print("BCC"); a=lc[0]
      cell,coords = bcc(a)
    elif args.structure is 3:
      print("FCC"); a=lc[0]
      cell,coords = fcc(a)
    elif args.structure is 4 and len(lc)==1:
        print("GRAPHENE"); rcc=lc[0]; 
        cell,coords = graphene(rcc)
        # expand to orthorhombic cell
        if args.orthorhombic:
            cell,coords = expand((1,2,1),cell,coords)
            cell=(array([2*cell[0][0],0,0]),array([0,cell[1][1],0]),array([0,0,cell[2][2]]))
    elif args.structure is 4 and len(lc)==2:
        print("GRAPHITE"); rcc = lc[0]; d = lc[1]
        cell,coords = graphite(rcc,d)
        # expand to orthorhombic cell
        if args.orthorhombic:
            cell,coords = expand((1,2,1),cell,coords)
            cell=(array([2*cell[0][0],0,0]),array([0,cell[1][1],0]),array([0,0,cell[2][2]]))
    elif args.structure is 5:
      print("DIAMOND"); a=lc[0]
      cell,coords = diamond(a,array([0,0,0]))
    elif args.structure is 6 and len(lc)==2:
      print("WURTZITE"); a=lc[0]; c=lc[1]
      cell,coords = wurtzite(a,c)
    else:
      print('structure not implemented or incorrect number of lattice constants'); return
  
    # expand crystal (n,m,k) times
    cell,coords = expand(args.ncells,cell,coords)

    # define atom types
    natoms = len(coords)
    types = [args.atomtype]*natoms
    if args.diatomic: 
      for i in range(1,natoms,2): types[i] = args.diatomic

    # compute PBCs
    (a1,a2,a3) = cell
    a = sqrt(dot(a1,a1)); 
    b = sqrt(dot(a2,a2)); 
    c = sqrt(dot(a3,a3)); 
    if(a==0): a=1
    if(b==0): b=1
    if(c==0): c=1
    #print(a,b,c)
    alpha = arccos(dot(a2,a3)/(b*c))*180/pi
    beta  = arccos(dot(a1,a3)/(a*c))*180/pi
    gamma = arccos(dot(a1,a2)/(a*b))*180/pi
    #print(dot(a1,a2),(a*b))
    #print(dot(a1,a2)/(a*b))
    #print(arccos(dot(a1,a2)/(a*b)))
    #print(arccos(dot(a1,a2)/(a*b))*180/pi)
    print('-----------------------------------------')
    print('created ',len(coords),'atom sample')
    print('(super)cell vectors')
    print("a1 = array([% 10.6f % 10.6f % 10.6f])" % (a1[0],a1[1],a1[2]))
    print("a2 = array([% 10.6f % 10.6f % 10.6f])" % (a2[0],a2[1],a2[2]))
    print("a3 = array([% 10.6f % 10.6f % 10.6f])" % (a3[0],a3[1],a3[2]))
    print('-----------------------------------------')
    print("cell for VMD : pbc set { %.3f %.3f %.3f   %.2f %.2f %.2f }; pbc box"%(a,b,c,alpha,beta,gamma))
    print("if the cell is not orthorhombic the sample must rotated first to align with the x-axis")

    if args.scaled:
      for c in coords:
        for d in range(3):
          c[d]/=cell[d][d]

    if   args.output is 1: writexyz(basename+".xyz",cell,coords,types)
    elif args.output is 2: writelammps(basename+".data",cell,coords)
    elif args.output is 3: writefor(basename+".for",basename+".sys",cell,coords,types)
    elif args.output is 4: writebs(basename+".bs",cell,coords,types)
    elif args.output is 5: writeobj(basename+".obj",coords)
    elif args.output is 6: writexsf(basename+".xsf",cell,coords,types)
    else: print('output method not implemented'); return


######### crystal structure functions ##########
def sc(a,r=array([0,0,0])): # simple cubic
    r0 = r
    cell = a*identity(3)
    return cell,[r0]

def bcc(a,r=array([0,0,0])): # body centered cubic
    r0 = r
    r1 = r + dot((a/2.),[1,1,1])
    cell = a*identity(3)
    return cell,[r0,r1]

def fcc(a,r=array([0,0,0])): # face centered cubic
    r0 = r
    r1 = r + dot((a/2.),[1,1,0])
    r2 = r + dot((a/2.),[1,0,1])
    r3 = r + dot((a/2.),[0,1,1])
    cell = a*identity(3)
    return cell,[r0,r1,r2,r3]

def graphene(rcc,r=array([0,0,0])): # rcc = CC distance in graphene
    # translation vectors for generating graphene:
    a1 = dot(rcc,[3/2.,-sqrt(3)/2,0])
    a2 = dot(rcc,[3/2.,+sqrt(3)/2,0])
    a3 = array([0,0,1])
    cell = (a1,a2,a3)

    # two atoms per unit cell
    t = array([rcc,0,0])
    r0 = r 
    r1 = r + t
    coords = [r0,r1]
    return cell,coords

def graphite(rcc,d,r=array([0,0,0])): # rcc = CC distance in graphene, d = interplanar distance
    # translation vectors for generating graphite:
    a1 = dot(rcc,[3/2.,-sqrt(3)/2,0])
    a2 = dot(rcc,[3/2.,+sqrt(3)/2,0])
    a3 = dot(d,[0,0,2])
    cell = (a1,a2,a3)

    # four atoms per unit cell
    t0 = dot(rcc,[1,0,0])
    t1 = dot(d,[0,0,1])
    r0 = r
    r1 = r0 + t0
    r2 = r1 + t1
    r3 = r2 + t0
    coords = [r0,r1,r2,r3] 
    return cell,coords

def diamond(a,r=array([0,0,0])):
    # diamond (or zinc-blende, = diamond with alternating atoms)   r_CC = sqrt(3)/4 * a        (3.46)
    r0 = r 
    r1 = r + dot((a/4.),[1,1,1])
    r2 = r + dot((a/4.),[2,2,0])
    r3 = r + dot((a/4.),[3,3,1])
    r4 = r + dot((a/4.),[0,2,2])
    r5 = r + dot((a/4.),[3,1,3])
    r6 = r + dot((a/4.),[2,0,2])
    r7 = r + dot((a/4.),[1,3,3])
    coords = [r0,r1,r2,r3,r4,r5,r6,r7]
    cell = a*identity(3)
    return cell,coords

def wurtzite(a,c,r=array([0,0,0])):
    # from http://cst-www.nrl.navy.mil/lattice/struk/b4.html
    u=3./8. 

    # primitive vectors
    a1 = dot(a/2.,[1,-sqrt(3),0])
    a2 = dot(a/2.,[1, sqrt(3),0])
    a3 = dot(c   ,[0,0,1])
    cell = (a1,a2,a3)

    # basis vectors
    b1 = dot(1/3.,a1) + dot(2/3.,a2)
    b2 = dot(1/3.,a1) + dot(2/3.,a2) + dot(u,a3)
    b3 = dot(2/3.,a1) + dot(1/3.,a2) + dot(1/2.,a3)
    b4 = dot(2/3.,a1) + dot(1/3.,a2) + dot(1/2.+u,a3)

    coords = [b1,b2,b3,b4]
    return cell,coords

######## Misc. Functions #########
def expand(n,cell,coords): 
    # expand current cell and coordinates by a1,a2,a3
    (n1,n2,n3) = n
    newcoords = []
    (a1,a2,a3) = cell
    for x in range(n1):
        for y in range(n2):
            for z in range(n3):
                r = dot(x,a1) + dot(y,a2) + dot(z,a3)
                for c in coords:
                    newcoords.append(c + r)

    newcell = (dot(n1,a1),dot(n2,a2),dot(n3,a3))
    return newcell,newcoords

def periodic(cell,coords):
    # move atoms into periodic cell box
    for c in coords:
        for i in range(3):
            c[i] = c[i] % cell[i]

def translate(r,coords):
    coords = [r+c for c in coords]

# Input/Output Functions
def writexyz(file,cell,coords,types):
    # xyz
    f=open(file,'w')
    f.write(str(len(coords))+'\n')
    f.write("%12.8f %12.8f %12.8f\n" % (cell[0][0], cell[1][1], cell[2][2]))
    for i,l in enumerate(coords):
        f.write("%-2s %12.8f %12.8f %12.8f\n" % (types[i],l[0], l[1], l[2]))  # xyz format

def writelammps(file,cell,coords): 
    # LAMMPS
    f=open(file,'w')
    f.write('info: \n')
    f.write('\n')
    f.write(str(len(coords))+' atoms\n')
    f.write('1 atom types\n')
    f.write('\n')
    f.write("0.0 %12.8f xlo xhi" % cell[0] + '\n')
    f.write("0.0 %12.8f ylo yhi" % cell[1] + '\n')
    f.write("0.0 %12.8f zlo zhi" % cell[2] + '\n')
    f.write('\n')
    f.write('Atoms\n')
    f.write('\n')
    for i,l in enumerate(coords):
        f.write("%3d 1 %12.8f %12.8f %12.8f\n" % (i+1, l[0], l[1], l[2])) # lammps format
    f.close()

def writefor(file,sysfile,cell,coords,types):
    # tersoff
    b=0.5291772083
    g=open(sysfile,'w')
    g.write("&SYSTEM\n")
    g.write("SYMMETRY\n")
    g.write("1\n")
    g.write("CELL\n")
    g.write("%7.5f %7.5f %7.5f %7.5f %7.5f %7.5f\n" % (cell[0]/b, cell[1]/cell[0], cell[2]/cell[0],cos(cell[3]),cos(cell[4]),cos(cell[5])))
    g.write("CUTOFF\n")
    g.write(" 220.0\n")
    g.write("&END\n")
    g.close()

    f=open(file,'w')
    f.write("#BEGIN ID=1 FRAME=1 SERIES=f ENERGY=0\n")
    f.write("#SETPBC=%s\n"%sysfile)
    for i,c in enumerate(coords):
        f.write("%4d %2s%8.4f%8.4f%8.4f  %10.3E %10.3E %10.3E\n" %(i+1,types[i],c[0]/b,c[1]/b,c[2]/b,cos(cell[3]),cos(cell[4]),cos(cell[5])))
    f.write("#END\n")
    f.close()

def writebs(file,Acell,coords,types):
    # xbs format 
    f=open(file,'w')
    for i,l in enumerate(coords):
        f.write("atom %s %12.8f %12.8f %12.8f\n" % (types[i],l[0], l[1], l[2]))  # xyz format
    f.write('\n')
    f.write('line 0.0 0.0 0.0  %10.8f %10.8f %10.8f\n' % (Acell[0][0],Acell[0][1],Acell[0][2]))
    f.write('line 0.0 0.0 0.0  %10.8f %10.8f %10.8f\n' % (Acell[1][0],Acell[1][1],Acell[1][2]))
    f.write('line 0.0 0.0 0.0  %10.8f %10.8f %10.8f\n' % (Acell[2][0],Acell[2][1],Acell[2][2]))
    f.write('\n')
    f.write('spec O  0.30 red\n')
    f.write('spec C  0.25 black\n')
    f.write('spec B  0.25 0.5\n')
    f.write('spec N  0.2  blue\n')
    f.write('spec P  0.35 green\n')
    f.write('spec H  0.15 yellow\n')
    f.write('\n')
    f.write('bonds C O 0 1.6 0.04 black\n')
    f.write('bonds C H 0 1.6 0.04 black\n')
    f.write('bonds O H 0 1.6 0.04 black\n')
    f.write('bonds O O 0 1.7 0.04 black\n')
    f.write('bonds H H 0 1.3 0.04 black\n')
    f.write('bonds C C 0 1.7 0.04 black\n')
    f.write('bonds B C 0 1.7 0.04 black\n')
    f.write('bonds B B 0 1.7 0.04 black\n')
    f.write('bonds N C 0 1.7 0.04 black\n')
    f.close()

def writeobj(file,coords):
    f=open(file,'w')
    f.write("o crystal\n")
    for c in coords:
      f.write("v %f %f %f\n" % (c[0],c[1],c[2]))
    for i,u in enumerate(coords):
      for j,v in enumerate(coords):
        if i!=j: 
          dx=sqrt(sum((u[k]-v[k])**2 for k in range(3)))
          if dx<1.7: f.write("l %d %d\n" % (i+1,j+1))
    f.close()

def writexsf(file,cell,coords,types):
    f = open(file,'w')
    natoms = len(coords)
    f.write("ANIMSTEP 1\nCRYSTAL\nPRIMVEC\n")
    #v1 = cell[:,0]; v2 = cell[:,1]; v3 = cell[:,2]
    f.write("  %12.8f %12.8f %12.8f\n"%(cell[0][0],cell[0][1],cell[0][2]))
    f.write("  %12.8f %12.8f %12.8f\n"%(cell[1][0],cell[1][1],cell[1][2]))
    f.write("  %12.8f %12.8f %12.8f\n"%(cell[2][0],cell[2][1],cell[2][2]))
    f.write("PRIMCOORD %d\n%d 1\n"%(1,natoms))
    for j in range(natoms):
      c = coords[j]
      f.write("%d  % 8.4f % 8.4f % 8.4f  % 8.4f % 8.4f % 8.4f\n"%\
        (atomnumber[types[j]],c[0],c[1],c[2],0,0,0))

# element names
elements = ["X","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po", "At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg"]
# masses in a.m.u.
masses = [0.00000, 1.00794, 4.00260, 6.941, 9.012182, 10.811, 12.0107, 14.0067, 15.9994, 18.9984032, 20.1797, 22.989770,24.3050, 26.981538, 28.0855, 30.973761, 32.065, 35.453,39.948, 39.0983, 40.078, 44.955910, 47.867, 50.9415,51.9961, 54.938049, 55.845, 58.9332, 58.6934, 63.546,65.409, 69.723, 72.64, 74.92160, 78.96, 79.904, 83.798,85.4678, 87.62, 88.90585, 91.224, 92.90638, 95.94, 98.0,101.07, 102.90550, 106.42, 107.8682, 112.411, 114.818,118.710, 121.760, 127.60, 126.90447, 131.293, 132.90545,137.327, 138.9055, 140.116, 140.90765, 144.24, 145.0,150.36, 151.964, 157.25, 158.92534, 162.500, 164.93032,167.259, 168.93421, 173.04, 174.967, 178.49, 180.9479,183.84, 186.207, 190.23, 192.217, 195.078, 196.96655,200.59, 204.3833, 207.2, 208.98038, 209.0, 210.0, 222.0,223.0, 226.0, 227.0, 232.0381, 231.03588, 238.02891,237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0,258.0, 259.0, 262.0, 261.0, 262.0, 266.0, 264.0, 269.0,268.0, 271.0, 272.0]
# atomic radii
radii = [1.5, 1.2, 1.4, 1.82, 2.0, 2.0, 1.7, 1.55, 1.52,1.47, 1.54, 1.36, 1.18, 2.0, 2.1, 1.8, 1.8, 2.27, 1.88, 1.76,1.37, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.63, 1.4, 1.39, 1.07,2.0, 1.85, 1.9, 1.85, 2.02, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,2.0, 2.0, 1.63, 1.72, 1.58, 1.93, 2.17, 2.0, 2.06, 1.98, 2.16,2.1, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 1.72, 1.66,1.55, 1.96, 2.02, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,1.86, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
mass = dict(zip(elements,masses))
radius = dict(zip(elements,radii))
atomnumber = dict(zip(elements,range(len(elements))))

if __name__ == "__main__":
    sys.exit(main())
