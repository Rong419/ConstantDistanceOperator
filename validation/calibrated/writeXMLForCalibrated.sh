#!/bin/sh
# This shell script is used to generate 100 xml files
# by replacing values of some parameters 
# including tree, frequencies, rates, kappa, ucld

for taxa in {20,120}
do

   TEMPLATE=cal_val_${taxa}_template.xml

   for sim in {1..100}
   do

   tree=$( sed -n ${sim}p ./true/${taxa}taxa/trees.txt)

   freq=$( sed -n ${sim}p ./true/${taxa}taxa/Freq.txt)

   kappa=$( sed -n ${sim}p ./true/${taxa}taxa/Kap.txt)

   ucld=$( sed -n ${sim}p ./true/${taxa}taxa/Ucld.txt)

   sed "s/TREE/${tree}/g; s/FREQUENCIES/${freq}/g; s/KAPPA/${kappa}/g; s/UCLDSTD/${ucld}/g; s/SIM/${sim}/g" ./${TEMPLATE} > ./xml/Calibrated_${taxa}taxa_${sim}.xml

   done
done

