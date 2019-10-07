echo "Times" > look.txt
for data in {anolis,Shankarappa,RSV2}
do
	for param in {1..2}
	do
		for sim in {1..3}
		do
		date +"        Start ${data}_${param}_${sim}.xml at %D %H:%M:%S"
		date +"Start ${data}_${param}_${sim}.xml at %D %H:%M:%S" >> look.txt
		{ time java -jar ../ConstantDistanceOperator.jar -threads 1 -beagle_SSE ./xml/${data}_${param}_${sim}.xml > ./Output/output_${data}_${param}_${sim}.txt; } 2>./Output/error_${data}_${param}_${sim}.txt
		date +"        End ${data}_${param}_${sim}.xml at %D %H:%M:%S"
		date +"End ${data}_${param}_${sim}.xml at %D %H:%M:%S" >> look.txt
		done
	done
done
