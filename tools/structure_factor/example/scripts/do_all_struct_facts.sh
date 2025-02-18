#!/bin/bash

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
SCRIPTPATH2="$( echo $SCRIPTPATH | rev | cut -d'/' -f3- | rev )"

for i in */config_*.xyz; 
do 
  python3 ${SCRIPTPATH2}/custom_scattering_calc_sq_avg_2d.py $i; 
done;

#accumulate all avg.x.y1.z1
for i in */; 
do 
  cd $i;
  pwd
  mx=;
  for j in config_*.x.y1.z1.Sq_3D.dat; do mx+="${j}:"; done;
  mx=${mx::-1};
  echo $mx;
  python3 ${SCRIPTPATH2}/accum_spaced_avg_2d.py $mx avg.x.y1.z1.3D.dat;
  cd ../;
done;

#accumulate all avg.x1.y.z1
for i in */; 
do 
  cd $i;
  my=;
  for j in config_*.x1.y.z1.Sq_3D.dat; do my+="${j}:"; done;
  my=${my::-1};
  python3 ${SCRIPTPATH2}/accum_spaced_avg_2d.py $my avg.x1.y.z1.3D.dat;
  cd ../;
done;

#accumulate all avg.x1.y1.z
for i in */; 
do 
  cd $i;
  mz=;
  for j in config_*.x1.y1.z.Sq_3D.dat; do mz+="${j}:"; done;
  mz=${mz::-1};
  python3 ${SCRIPTPATH2}/accum_spaced_avg_2d.py $mz avg.x1.y1.z.3D.dat;
  cd ../;
done;

#plot all three structure factor patterns and combine into single pdf (tmp.pdf)
python3 ${SCRIPTPATH}/plotIMshow.py */avg.x.y1.z1.3D.dat */avg.x1.y.z1.3D.dat */avg.x1.y1.z.3D.dat
pdfjam */avg*.png --nup 3x1 --landscape --outfile tmp.pdf;

bash ../scripts/copy_stuff.sh 
python3 ../scripts/subfig-struct_facts.py
