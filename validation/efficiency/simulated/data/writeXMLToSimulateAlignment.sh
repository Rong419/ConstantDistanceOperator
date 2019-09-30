#!/bin/sh
# This shell script is used to generate 100 xml files
# by replacing values of some parameters 
# including tree, frequencies, rates, kappa, ucld

for taxa in {20,120}
do

   TEMPLATE=SimulateAlignment${taxa}taxa_template.xml

   for length in {5000,10000,20000}
   do
   
   sed "s/LENGTH/${length}/g" ./${TEMPLATE} > ./xml/${taxa}taxaAlignmentLength_${length}.xml

   done
done

