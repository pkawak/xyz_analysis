#!/bin/bash

# prep file in opwl-xyz format
line1=$(grep -n "Atoms" nvt.postcross.data | cut -d':' -f1)
line2=$(grep -n "Veloc" nvt.postcross.data | cut -d':' -f1)
sed -n "$((${line1}+2)),$((${line2}-2))p;$((${line2}-2))q" nvt.postcross.data > nvt.postcross.data.opwl
sort -n nvt.postcross.data.opwl > nvt.postcross.data.opwl2
sed "/^[0-9]* [0-9]* [6|8]/d" nvt.postcross.data.opwl2 > nvt.postcross.data.opwl
sed -E -i 's/^[0-9]* //' nvt.postcross.data.opwl 
sed -E -i 's/^([0-9]*) [0-9]* /\1 /' nvt.postcross.data.opwl 
sed -i 's/ [-|0-9]* [-|0-9]* [-|0-9]*$/ 0/g' nvt.postcross.data.opwl 
N=$(cat nvt.postcross.data.opwl | wc -l)

# get some params from orig file
Nb=$(grep "^1 " nvt.postcross.data.opwl | wc -l)
Nc=$(echo "$N/$Nb" | bc -l | cut -d'.' -f1)
low_lx=$(grep xlo nvt.postcross.data | cut -d' ' -f1)
hih_lx=$(grep xlo nvt.postcross.data | cut -d' ' -f2)
lx=$(python3 -c "print($hih_lx-$low_lx)")

# make coordinates' origin 0
mid_lx=$(python3 -c "print(($hih_lx+$low_lx)/2)")
awk -v mid_lx="$mid_lx" '{$2=$2-mid_lx ; $3=$3-mid_lx ; $4=$4-mid_lx ; print }' nvt.postcross.data.opwl > nvt.postcross.data.opwl2
cp nvt.postcross.data.opwl2 nvt.postcross.data.opwl

# add two comment lines
echo "//" > nvt.postcross.data.opwl2
echo "//" >> nvt.postcross.data.opwl2
cat nvt.postcross.data.opwl >> nvt.postcross.data.opwl2
cp nvt.postcross.data.opwl2 nvt.postcross.data.opwl
rm nvt.postcross.data.opwl2

# need to make a stand in params.json for opwl to run
# solid cut off
rs=1.0
# choose max number of cells that is also even
ncell=$(echo $(($(($(echo "${lx}/${rs}" | bc -l | cut -d'.' -f1)/2))*2)))
Nmax=$(echo "10000/22/22/22*20" | bc -l | cut -d'.' -f1)
printf "{
    \"Nc\": ${Nc},
    \"Nb\": ${Nb},
    \"Lx\": ${lx},
    \"rs\": ${rs},
    \"ncell\": ${ncell},
    \"Nmax\": ${Nmax}
}" > params.json
