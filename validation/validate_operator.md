# conoperators: Validation of ConstantDistance Operator

This file will guide you through reproducing the simulations of validating the ConstantDistance Operator, as is detailed in the paper.

## 1. Sample from prior
The tree includes three taxa.
```
cd /validation/sample_prior/
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/sample_prior/test_internalnode.xml
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/sample_prior/test_simpledistance.xml
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/sample_prior/test_smallpulley.xml
java -jar /path/to/jar_file/ConstantDistanceOperator.jar /validation/sample_prior/test_bigpulley.xml
```

## 2. Well-calibrated simulation study

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
## 3. Efficiency comparison
## (3.1) Simulated data sets
### 20 taxa
### 120 taxa

## (3.2) Real data sets
### primates data set