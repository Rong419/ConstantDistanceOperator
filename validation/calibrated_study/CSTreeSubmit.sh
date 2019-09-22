#!/bin/sh
TEMPLATE=Treetemplate.sh

for taxa in {20,120}
do
    for sim in {1..100}
    do
    sed "s/FILE/calibrated${taxa}taxa_${sim}/g" ${TEMPLATE} > ./Treerandom.sl 
    echo "summarise calibrated${taxa}taxa_${sim}.trees"
    sbatch Treerandom.sl 
    rm -f Treerandom.sl 
    sleep 5
    done
  done
#done
