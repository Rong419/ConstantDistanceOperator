# consoperators: Validation of ConstantDistance Operator

This file will guide you through reproducing the simulations of validating the ConstantDistance Operator, as is detailed in the paper.

## 1. Sample from prior
In this step, we aim at get samples from prior distributions with the operators working internal nodes and root separately. Then, we demonstrate the correctness of the operators by comparing the sampled distributions with the numerical results.
### Run BEAST analysis
To test the operator working on internal nodes, we have designed two scenarios with different trees and initial rates. In each scenario, we have ran the simulations with two different MCMC chain lengths. The script "write_xml_test_internalnode.R" is used to generate the corresponding xml files.
```
cd /validation/sampleprior/internal_nodes/
Rscript write_xml_test_internalnode.R /validation/sampleprior/internal_nodes/test_internalnode_template.xml /validation/sampleprior/internal_nodes/xml/ /validation/sampleprior/internal_nodes/test_internalnode_trees.txt /validation/sampleprior/internal_nodes/test_internalnode_rates.txt 10000000 2
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/sampleprior/internal_nodes/xml/internalnode_S1_1.xml
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/sampleprior/internal_nodes/xml/internalnode_S1_2.xml
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/sampleprior/internal_nodes/xml/internalnode_S2_1.xml
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/sampleprior/internal_nodes/xml/internalnode_S2_2.xml
```

For operators working the root, we directly run the xmls files in folder "/validation/sample_prior".
```
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/sampleprior/test_simpledistance.xml
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/sampleprior/test_smallpulley.xml
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/sampleprior/test_bigpulley.xml
```

### Conduct numerical integration and make comparisons
The scripts in folder "/validation/r_scripts/sample_prior/" will help us get the theoretical distribution of the sampled parameters and produce the figures to compare with the sampled distributions in .log and .trees files that are produced by BEAST2.
```
cd /validation/sampleprior/
Rscript internal_node.R "/validation/sampleprior/internal_nodes/test_internalnode_trees.txt" "/validation/sampleprior/internal_nodes/test_internalnode_rates.txt" "/ConstantDistanceOperator/out/artifacts/ConstantDistanceOperator_jar/" "/validation/sampleprior/"
Rscript simple_distance.R "/ConstantDistanceOperator/out/artifacts/ConstantDistanceOperator_jar/" "/validation/sampleprior/"
Rscript small_pulley.R "/ConstantDistanceOperator/out/artifacts/ConstantDistanceOperator_jar/" "/validation/sampleprior/"
Rscript big_pulley.R "/ConstantDistanceOperator/out/artifacts/ConstantDistanceOperator_jar/" "/validation/sampleprior/"
```


## 2. Well-calibrated simulation study
We further verify the ConstantDistance Operator by a well-calibrated simulation study.
### (2.1) Get samples of parameters from prior distributions
```
cd /validation/calibrated
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/calibrated/getPriorSamplesFor20Taxa.xml
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/calibrated/getPriorSamplesFor120Taxa.xml
```
### (2.2) Write .xml files for simulation
Run R scripts to get true values of parameters (true tree) from prior samples.
```
cd /validation/calibrated 
Rscript getTrueValue.R "validation/calibrated/priors/" "validation/calibrated/true/" 
Rscript getTrueTree.R "validation/calibrated/priors/" "validation/calibrated/true/" 
```

Run shell script to write .xml files by using the templates.
```
cd /validation/calibrated 
./writeXMLForCalibrated.sh
```

### (2.3) Run BEAST analysis to sample the simulated data
After sending the .xml files to the cluster (NeSI), we used the following command line and shell scripts (CStemplate.sh and submit_CSjob.sh) to submit the simulations.
```
cd /path/to/cluster/directoty
./submit_CSjob.sh
```
### (2.4) Analyse the results of simulations
Run R script to compare mean of sampled distributions to the real values.
```
Rscript OperatorCalibratedPlot.R 100 20 "/validation/calibrated/logs/" "/validation/calibrated/true/" "/validation/calibrated/figures/"
Rscript OperatorCalibratedPlot.R 100 120 "/validation/calibrated/logs/" "/validation/calibrated/true/" "/validation/calibrated/figures/"
```

## 3. Efficiency comparison
## (3.1) Simulated data sets
### 20 taxa
### 120 taxa
Sequence length 5000, 10000, 20000
```
cd /validation/simulated_data/20taxa/data
cd /validation/simulated_data/120taxa/data
```

Using categories (Alignment1.xml), continuous rates with ConstantDistance operator (Alignment2.xml) and continuous rates without ConstanceDistance operator (Alignment3.xml) .
```
cd /validation/simulated_data/20taxa/xmls/categories
cd /validation/simulated_data/20taxa/xmls/cons
cd /validation/simulated_data/20taxa/xmls/nocons
```
```
cd /validation/simulated_data/120taxa/xmls/categories
cd /validation/simulated_data/120taxa/xmls/cons
cd /validation/simulated_data/120taxa/xmls/nocons
```
## (3.2) Real data sets
### primates data set
Using categories (primates1.xml), continuous rates with ConstantDistance operator (primates2.xml) and continuous rates without ConstanceDistance operator (primates3.xml).
```
cd /validation/primates_data/
```
Run LogAnalyser to get ESS.
```
/Applications/BEAST\ 2.4.7/bin/loganalyser -oneline /Users/rzha419/Desktop/efficiency/RSV2/logs/*.log
```

## 4. Correlation analysis
We use a data set with sequence of seven ratites to study the relationships among rates and node times.
```
cd /validation/correlation/
```
### (4.1) Run BEAST analysis using bModelTest as substitution model
```
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/correlation/ratites.xml
```

### (4.2) Extract rates and node times from MCMC samples
After BEAST2 finishes running the .xml, we will get "ratites.log" and "ratites.trees". Use TreeAnnotator to summarise the "ratites.trees", which provides us the most credibility tree in "summary_ratites.trees".

Then, according to the clades in "summary_ratites.trees", we use TreeStat to filter all the tree samples in "ratites.trees", so that we get the rates and node times in the trees that have the same clades as the most credibility tree.

Finally, the node times, internal rates and external rates are written in "IntRates.log" and "ExtRates.log".

### (4.3) Make figures of correlation
Run the R script and we will get two separate figures, i.e., "RatesAndTimeCorrelation.pdf" and "DistanceCorrelation.pdf".
```
Rscript CorrelationFigurePlot.R "/validation/correlation/" "~/Desktop/validation/correlation/figures/"
```

## 5. Fixed tree analysis
In this section, we try to test the ConstantDisance operator by using a fixed tree whose distances are given by a maximum likelihood method. 
### (5.1) Run online software 
We 

Unrooted tree to a rooted time tree.

### (5.2) Run BEAST analysis by sampling the fixed tree using ConstantDistance operator
```
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/fixedtree/FixedRatites.xml
```

### (5.3) Analyse trees
Output "FixedRatites.trees". TreeTraceAnalysis