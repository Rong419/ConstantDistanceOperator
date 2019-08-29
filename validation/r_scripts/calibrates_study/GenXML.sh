#!/bin/sh
# This shell script is used to generate 100 xml files
# by replacing values of some parameters 
# including tree, frequencies, rates, kappa, ucld

TEMPLATE=testSimulatedAlignment.xml

for file in {1..100}
do

tree=$( sed -n ${file}p /Users/rzha419/Desktop/tree.txt)

freq=$( sed -n ${file}p /Users/rzha419/Desktop/Freq.txt)

rates=$( sed -n ${file}p /Users/rzha419/Desktop/Rates.txt)

kappa=$( sed -n ${file}p /Users/rzha419/Desktop/Kap.txt)

ucld=$( sed -n ${file}p /Users/rzha419/Desktop/Ucld.txt)

ratm=$( sed -n ${file}p /Users/rzha419/Desktop/RatM.txt)

sed "s/TREE/${tree}/g" ./${TEMPLATE} > ./temp1.xml

sed "s/FREQUENCIES/${freq}/g" ./temp1.xml > ./temp2.xml

sed "s/KAPPA/${kappa}/g" ./temp2.xml > ./temp3.xml

sed "s/RATES/${rates}/g" ./temp3.xml > ./temp4.xml

sed "s/RATEMEAN/${ratm}/g" ./temp4.xml > ./temp5.xml

sed "s/UCLDSTD/${ucld}/g" ./temp5.xml > ./XML/Alignment${file}.xml

done

