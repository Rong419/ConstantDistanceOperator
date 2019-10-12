#!/bin/sh

for data in {anolis,Shankarappa,RSV2}
do

	for model in {Cons,Category}
	do
	/Users/rzha419/Applications/BEAST2.6.0/bin/loganalyser -oneline /Users/rzha419/Desktop/validation/efficiency/others/${data}${model}/logs/*.log >/Users/rzha419/Desktop/validation/efficiency/others/ess/ESS_${data}${model}.txt
	done
	
done

