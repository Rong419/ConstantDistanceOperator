#!/bin/sh
TEMPLATE=Treetemplate.sl

for ft in {1..3}
do  
for sim in {1..3}
  do 
  sed "s/FILE/${ft}_${sim}/g" ${TEMPLATE} > ./Treerandom.sl 
    echo "summarise primates${ft}_${sim}.trees"
    sbatch Treerandom.sl 
    rm -f Treerandom.sl 
    sleep 5
 done
done
