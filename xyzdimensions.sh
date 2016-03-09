#!/bin/sh
awk '
BEGIN{xmin="";xmax="";ymin="";ymax="";zmin="";zmax=""}
{ 
  if(NR>2) {
    if($2<xmin || xmin == "") xmin=$2;
    if($2>xmax || xmax == "") xmax=$2;
    if($3<ymin || ymin == "") ymin=$3;
    if($3>ymax || ymax == "") ymax=$3;
    if($4<zmin || zmin == "") zmin=$4;
    if($4>zmax || zmax == "") zmax=$4;
    cm[0]+=$2
    cm[1]+=$3
    cm[2]+=$4
    n+=1
  }
}
END{
  printf "dx,dy,dz = %.3f %.3f %.3f\n",xmax-xmin,ymax-ymin,zmax-zmin
  printf "CM       = %.3f %.3f %.3f\n",cm[0]/n,cm[1]/n,cm[2]/n}' $1
