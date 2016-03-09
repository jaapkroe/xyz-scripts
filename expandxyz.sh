#!/bin/bash

nrep="1 1 1"

usage="${0##*/} [options] file\n
  -n  list  -- number of repetitions \"N[1] N[2] N[3]\"
  -p  list  -- define PBCs : \"Lx Ly Lz\"
"

while getopts "n:p:vh" Option; do
  case $Option in
    n) nrep="$OPTARG";;
    p) pbc="$OPTARG";;
    v) ((verbose+=1));;
    *) echo -e "$usage" > "/dev/stderr"; exit;;
  esac
done
shift $(($OPTIND - 1))

if [ $# -eq 0 ]; then while read arg; do echo $arg; done
else 
  for f in $*; do
    if [ $# -gt 1 ]; then echo "#reading# $f" > "/dev/stderr"; fi
    if [ "${f##*.}" == "gz" ]; then gunzip -c $f ; else cat $f; fi;
  done 
fi | awk -v nrep="$nrep" -v pbc="$pbc" 'BEGIN{
  split(nrep,N," ");
  if(pbc) split(pbc,L," ");
  } {
  if(NR==1) natom=$1
  else if(NR==2) {
    if(!pbc) { L[1]=$1;L[2]=$2;L[3]=$3 }
    if(L[1]*L[2]*L[3]==0) { print "error in PBCs"; exit }
    if(N[1]*N[2]*N[3]<1) { print "error in N (<1?)"; exit }
    print N[1]*N[2]*N[3]*natom; 
    print N[1]*L[1],N[2]*L[2],N[3]*L[3]; 
  }
  else { 
    for(nx=0;nx<N[1];nx++) { 
      for(ny=0;ny<N[2];ny++) { 
        for(nz=0;nz<N[3];nz++) { 
          printf "%s %14.10f %14.10f %14.10f\n",$1,$2+nx*L[1],$3+ny*L[2],$4+nz*L[3]
        }
      }  
    }
  }
}' 

