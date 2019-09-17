#!/bin/bash

#
# nei's Distance
# D=-log(10,I), I=0,0.5,1
# ${i}.gb are GB genotypes of individual i
#

cd /home/myproject/01NJtree
for i in {3..407}
do
for j in {3..407}
do
paste ${i}.gb ${j}.gb > ${i}${j}.gb.bed
grep -v "m" ${i}${j}.gb.bed > ${i}${j}.gb.tmp
awk '{if(($1==$3 && $2==$4) || ($1==$4 && $2==$3)) print 1; else if($1!=$3 && $2!=$4 && $1!=$4 && $2!=$3) print 0; else print 0.5}' ${i}${j}.gb.tmp > nes.${i}${j}.nes
awk 'BEGIN{total=0}{total+=$1}END{print -log(total/NR)}' nes.${i}${j}.nes >> nes.tmp.all
rm ${i}${j}.gb.bed
rm ${i}${j}.gb.tmp
rm nes.${i}${j}.nes
done
done
