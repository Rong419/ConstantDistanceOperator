#!/bin/sh
TEMPLATE=Treetemplate.sl

for data in {anolis,Shankarappa,RSV2}
do
  for param in {1..2}
  do
    for sim in {1..3}
    do
    sed "s/FILE/${data}_${param}_${sim}/g" ${TEMPLATE} > ./Treerandom.sl 
    echo "summarise ${data}_${param}_${sim}.trees"
    sbatch Treerandom.sl 
    rm -f Treerandom.sl 
    sleep 5
    done
  done
done
