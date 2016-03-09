#!/bin/bash

rcut=1.7
angles=0
bonds=1
stat=0
acc=3
numbers=1
linebreak=0
verbose=0

usage="${0##*/} [options] file\n
SELECTION
  -r  float -- distance cutoff                              ($rcut)
  -p  list  -- define PBCs : \"Lx Ly Lz\"
  -i  int   -- from index [1..n] (may also be a list, e.g. \"4 5 9\")
  -j  int   -- to   index [1..n]
  -t  str   -- from type (H, C, ..)
  -u  str   -- to   type

ANALYSIS
  -q        -- compute angles                               ($angles)
  -b        -- compute bond lengths                         ($bonds)
  -m        -- minimal value                                ($bonds)
  -s        -- statistics (number of bonds per atom)        ($stat)

PRINTING
  -a  int   -- accuracy for r (2 digits less for angles)    ($acc)
  -n        -- toggle display of indices in output          ($numbers)
  -l        -- toggle linebreak for long output             ($linebreak)
  -d int    -- print distance projected on axis (1/2/3)"

while getopts "r:i:j:t:u:x:p:qmd:bsa:nlvh" Option; do
  case $Option in
    r) rcut="$OPTARG";;
    i) ind="$OPTARG";;
    j) jnd="$OPTARG";;
    t) typ1="$OPTARG";;
    u) typ2="$OPTARG";;
    x) typx="$OPTARG";;
    q) angles=$((1-angles));;
    m) minimal=1;;
    d) direction=$OPTARG;;
    b) bonds=$((1-bonds));;
    a) acc=$OPTARG;;
    n) numbers=$((1-numbers));;
    l) linebreak=$((1-linebreak));;
    s) stat=$((1-stat)); numbers=0;;
    p) pbc="$OPTARG";;
    v) ((verbose+=1));;
    *) echo -e "$usage" > "/dev/stderr"; exit;;
  esac
done
shift $(($OPTIND - 1))

if [ $# -eq 0 ]; then while read arg; do echo $arg; done
else 
  for f in $*; do
    if [[ $# -gt 1 && $verbose>0 ]]; then echo "reading: $f" > "/dev/stderr"; fi
    if [ "${f##*.}" == "gz" ]; then gunzip -c $f ; else cat $f; fi;
  done 
fi | awk -v rcut="$rcut" -v ind="$ind" -v jnd="$jnd" -v acc="$acc" -v typ1="$typ1" -v typ2="$typ2" -v typx="$typx" -v angles="$angles" -v bonds="$bonds" -v numbers="$numbers" -v minimal="$minimal" -v stat="$stat" -v pbc="$pbc" -v direction="$direction" -v linebreak="$linebreak" 'BEGIN{
  i=1; rcut2=rcut*rcut;
  if(pbc) split(pbc,L," ");
  getline; n=$1; 
  if(ind) split(ind,ilist," ")
  if(jnd) split(jnd,jlist," ")
}
function wrap(x,L) { return x= (x<0 ? x-int(x/L-0.5)*L : x-int(x/L+0.5)*L) } # awk lacks a proper NINT function
{
  prev=$0; getline
  while(prev!=$0) {
    for(i=1;i<=n;i++) {t[i]=$1; x[i]=$2; y[i]=$3; z[i]=$4; getline}; 
    printf "%-3d ",frame++
    if(!ind) { for(i=1;i<=n;i++) ilist[i] = i }
    for(ii in ilist) { i=ilist[ii];
      if(minimal) rmin=""
      if(stat) { nbonds=0; nangles=0; }
      if(typ1 && typ1!=t[i]) continue
      lowj=(jnd>=0?jnd:(typ1||ind>=0?1:i+1));
      higj=(jnd>=0?jnd:n);
      if(!jnd) { jmin=((typ1||ind)?1:i+1)
                 for(j=jmin;j<=n;j++) jlist[j] = j}
      for(jj in jlist) { j=jlist[jj]
        if(typ2 && typ2!=t[j]) continue
        if(typx && typx==t[j]) continue
        if(typ1 == t[j] && j<i) continue
        dx[1]=x[j]-x[i]
        dx[2]=y[j]-y[i]
        dx[3]=z[j]-z[i]
        if(pbc) { dx[1] = wrap(dx[1],L[1]) ; dx[2] = wrap(dx[2],L[2]) ; dx[3] = wrap(dx[3],L[3]) }
        r2 = dx[1]*dx[1] + dx[2]*dx[2] + dx[3]*dx[3]
        if(minimal) { if(rmin=="" || r2<rmin) rmin=r2 }
        if(r2<rcut2 && i!=j) {
          if(angles) {
            for(kk in jlist) { k=jlist[kk]
              if(k<j&&k!=i) {
                dq[1]=x[k]-x[i]
                dq[2]=y[k]-y[i]
                dq[3]=z[k]-z[i]
                if(pbc) { dq[1] = wrap(dq[1],L[1]) ; dq[2] = wrap(dq[2],L[2]) ; dq[3] = wrap(dq[3],L[3]) }
                r2_k = dq[1]*dq[1] + dq[2]*dq[2] + dq[3]*dq[3]
                if(r2_k<rcut2) {
                  theta =(x[j]-x[i])*(x[k]-x[i])
                  theta+=(y[j]-y[i])*(y[k]-y[i])
                  theta+=(z[j]-z[i])*(z[k]-z[i])
                  theta/=sqrt(r2*r2_k)
                  theta=atan2(sqrt(1.-theta*theta),theta)
                  theta*=180./3.14159265359
                  if(numbers) printf "%s%d-%s%d-%s%d ",t[j],j,t[i],i,t[k],k
                  if(stat) nangles+=1;
                  else printf "%.*f ",(acc<2?0:acc-2),theta
                }
              }
            }
          }
          if(bonds && !minimal) {
            if(numbers) printf "%s%d-%s%d ",t[i],i,t[j],j
            if(stat) nbonds+=1;
            else if(direction>=1 && direction<=3) printf "%*.*f ",acc+3,acc,dx[direction]
            else printf "%.*f ",acc,sqrt(r2)
          }
        }
      }
      if(minimal) printf "%.*f ",acc,sqrt(rmin)
      if(!jnd) delete jlist
      if(linebreak) printf "\n"
      if(stat) { printf "%d ",nbonds; if(angles) printf "%d ",nangles }
    }
    printf "\n"
    n=$1; getline; prev=$0; getline;
  }
}' 
