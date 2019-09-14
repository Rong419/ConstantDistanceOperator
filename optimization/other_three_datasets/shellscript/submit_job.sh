#!/bin/sh
TEMPLATE=template.sh

for data in {anolis,Shankarappa,RSV2}
do  
 for param in {1..2}
 do 
    for sim in {1..3}
    do
    sed "s/DATA/${data}/g; s/PARAM/${param}/g; s/SIM/${sim}/g" ${TEMPLATE} > ./temp.sl 
    echo "submit ${data}_${param}_${sim}.xml"
    sbatch temp.sl 
    rm -f temp.sl 
    sleep 5
    done
  done
done
