echo "Times" > look.txt
for taxa in {20,120}
do
	for model in {Cons,Category}
	do
		for length in {Long,Short,Medium}
		do
			for sim in {1..10}
			do
			date +"        Start xml at %D %H:%M:%S"
			date +"Start .xml at %D %H:%M:%S" >> look.txt
			{ time java -jar ../ConstantDistanceOperator.jar -threads 1 -beagle_SSE ./xml/ShortCons20taxa_${sim}.xml > ./Output/output_${length}${model}${taxa}taxa_${sim}.txt; } 2>./Output/error_${length}${model}${taxa}taxa_${sim}.txt
			date +"        End xml at %D %H:%M:%S"
			date +"End xml at %D %H:%M:%S" >> look.txt
			done
		done
	done
done