# conoperators: Validation of ConstantDistance Operator

This file will guide you through reproducing the simulations of validating the ConstantDistance Operator, as is detailed in the paper.

## 1. Sample from prior
In this step, we aim at get samples from prior distributions with the operators working internal nodes and root separately. Then, we demonstrate the correctness of the operators by comparing the sampled distributions with the numerical results.
### Run BEAST analysis
To test the operator working on internal nodes, we have designed two scenarios with different trees and initial rates. In each scenario, we have ran the simulations with two different MCMC chain lengths. The script "write_xml_test_internalnode.R" is used to generate the corresponding xml files.
```
cd /validation/sample_prior/
Rscript write_xml_test_internalnode.R /validation/sample_prior/internal_nodes/test_internalnode_template.xml /validation/sample_prior/internal_nodes/xml/ /validation/sample_prior/internal_nodes/test_internalnode_trees.txt /validation/sample_prior/internal_nodes/test_internalnode_rates.txt 10000000 2
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/sample_prior/internal_nodes/xml/internalnode_S1_1.xml
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/sample_prior/internal_nodes/xml/internalnode_S1_2.xml
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/sample_prior/internal_nodes/xml/internalnode_S2_1.xml
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/sample_prior/internal_nodes/xml/internalnode_S2_2.xml
```

For operators working the root, we directly run the xmls files in folder "/validation/sample_prior".
```
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/sample_prior/test_simpledistance.xml
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/sample_prior/test_smallpulley.xml
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/sample_prior/test_bigpulley.xml
```

### Conduct numerical integration and make comparisons
The scripts in folder "/validation/r_scripts/sample_prior/" will help us get the theoretical distribution of the sampled parameters and produce the figures to compare with the sampled distributions in .log and .trees files that are produced by BEAST2.
```
cd /validation/r_scripts/sample_prior/
Rscript internal_node.R
Rscript simple_distance.R
Rscript small_pulley.R
Rscript big_pulley.R
```


## 2. Well-calibrated simulation study
We further verify the ConstantDistance Operator by a well-calibrated simulation study.
### (2.1) Get samples of parameters from prior distributions
```
cd /validation/calibrated_study/20taxa
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/calibrated_study/20taxa/getPriorSamples.xml
cd /validation/calibrated_study/120taxa
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/calibrated_study/120taxa/getPriorSamples.xml
```
### (2.2) Simulate sequence alignment using parametric samples
```
cd /validation/calibrated_study/20taxa
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/calibrated_study/20taxa/.xml
cd /validation/calibrated_study/120taxa
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/calibrated_study/20taxa/.xml
```
### (2.3) Run BEAST analysis to sample the simulated data
```
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/calibrated_study/20taxa/.xml

```
### (2.4) Compare mean of sampled distributions to the real values
```
Rscript
```

## 3. Correlation analysis
We use a data set with sequence of seven ratites.
```
cd /validation/ratites_data/
```
### (3.1) Run BEAST analysis using bModelTest as substitution model
```
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/ratites_data/ratites1.xml
```
### (3.2) Run BEAST analysis using GTR+I+gamma as substitution model
```
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/ratites_data/ratites2.xml
```
### (3.3) Run BEAST analysis by sampling the fixed tree using ConstantDistance operator
```
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/ratites_data/ratites3.xml
```



## 4. Efficiency comparison
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
