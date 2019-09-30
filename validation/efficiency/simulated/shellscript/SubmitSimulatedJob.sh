#!/bin/sh
TEMPLATE=SimulatedTemplate.sh


for  taxa in {20,120}
do 
	for model in {Category,Cons}
	do
		for length in {Long,Medium,Short}
		do
			for sim in {1..10}
			do
			sed "s/FILE/${length}${model}${taxa}taxa_${sim}/g" ${TEMPLATE} > ./temp.sl 
			echo "submit ${length}${model}${taxa}taxa_${sim}.xml"
			sbatch temp.sl 
			rm -f temp.sl 
			sleep 5
			done
		done
	done
done
