#!/bin/sh

for data in {anolis,Shankarappa,RSV2}
do
	for model in {Cons,Category}
	do
	/Applications/BEAST2.6.0/bin/loganalyser -oneline ~/Desktop/validation/efficiency/others/${data}${model}/logs/*.log >~/Desktop/validation/efficiency/others/ess/ESS_${data}${model}.txt
	done	
done



for model in {Cons,Category}
do
	/Applications/BEAST2.6.0/bin/loganalyser -oneline ~/Desktop/validation/efficiency/simulated/20taxa/Short${model}/logs/*.log >~/Desktop/validation/efficiency/simulated/20taxa/ess/ESS_Short${model}20taxa.txt
	/Applications/BEAST2.6.0/bin/loganalyser -burnin 25 -oneline ~/Desktop/validation/efficiency/simulated/20taxa/Medium${model}/logs/*.log >~/Desktop/validation/efficiency/simulated/20taxa/ess/ESS_Medium${model}20taxa.txt
done



