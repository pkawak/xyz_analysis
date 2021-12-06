#!/bin/bash

#converts tag xyz, params, etc. files to needed format for this repo

chem_file=chem.in
tag_file=xyz.in
param_file=params.in

num=$(cat ${chem_file} | wc -l)
sed 's/e+00//g' ${chem_file} > ${chem_file}.int
sed -i -r -e 's/([0-9]+)\.0+\b/\1/g' -e 's/([0-9]+\.[0-9]*[1-9])0+\b/\1/g' ${chem_file}.int 
lx=$(grep Lx ${param_file} | cut -d' ' -f4 | cut -d ',' -f1)
paste ${chem_file}.int ${tag_file} > ${tag_file}.ident
sed -i 's/\t/ /g' ${tag_file}.ident
sed -i "1s/^/\/\/ ${num} 1 ${lx}\n/" ${tag_file}.ident
sed -i "1s/^/${num}\n/" ${tag_file}.ident
