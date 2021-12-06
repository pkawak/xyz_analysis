#!/bin/bash

~/data/crystal_play/struct_fact_stuff/plotIMshow.py ${1}.x.y1.z1.Sq_3D.dat
~/data/crystal_play/struct_fact_stuff/plotIMshow.py ${1}.x1.y.z1.Sq_3D.dat
~/data/crystal_play/struct_fact_stuff/plotIMshow.py ${1}.x1.y1.z.Sq_3D.dat
pdfjam ${1}.x.y1.z1.Sq_3D.dat.Sq_plot.png ${1}.x1.y.z1.Sq_3D.dat.Sq_plot.png ${1}.x1.y1.z.Sq_3D.dat.Sq_plot.png --nup 3x1 --landscape --outfile ${1}.pdf
