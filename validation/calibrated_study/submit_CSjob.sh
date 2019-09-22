#!/bin/sh
TEMPLATE=template.sh


#for  taxa in {20,120}
#do 
    taxa=20
    for sim in {1..100}
    do
    sed "s/TAXA/${taxa}/g; s/SIM/${sim}/g" ${TEMPLATE} > ./temp.sl 
    echo "submit Calibrated_${taxa}taxa_${sim}.xml"
    sbatch temp.sl 
    rm -f temp.sl 
    sleep 5
    done
#done

