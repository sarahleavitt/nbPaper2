# nbPaper2

This directory contains the code to produce the results for
"Estimation of the generation interval using pair-wise relative transmission probabilities" 
by Leavitt et al. It contains code to run the simulations assessing 
using naive Bayes transmission probabilities to estimate the generation and serial 
intervals and the application to TB surveillance data from the Massachusetts 
department of health (DPH). These programs use the nbTransmission R package available
on CRAN or GitHub (https://github.com/sarahleavitt/nbTransmission)
as well as the simulation programs in the nbSimulation directory
(https://github.com/sarahleavitt/nbSimulation).
 

## Simulation Programs

The following programs were used to run a simulation study to asses using 
the naive Bayes transmission probabilities to estimate the generation interval 
with different clustering methods.


### SimRunSIFull.R

This program contains a function "simRunSI" which performs one iteration of a 
simulation which simulates an outbreak with pathogen genetics and covariates,
estimates the relative transmission probability using naive Bayes, assesses the 
performance, and estimates the generation interval using various methods.


### PerformSimulationSIFull.R

This program is a wrapper for "simRunSI" which runs the simulation coded in that 
function nSim times. It saves the effect estimates, generation interval estimates, 
and performance for all iterations. It is meant to be run with many different 
outbreak scenarious and creates labels based on the inputs. This is run by 
"SimQsubSIFull.qsub" with 9 different inputs: sample size, reproductive number, 
generation interval (shape and scale), mutation rate, SNP thresholds (lower and upper), 
and initial generation interval parameters (shape and scale).


### SimQsubSIFull.qsub

This a qsub file to run PerformSimulationSIFull.R as a batch job. 
It runs t parallel jobs and therefore the total number of simulations is 
t*nSim times. It needs 9 inputs: sample size, reproductive number, 
generation interval (shape and scale), mutation rate, SNP thresholds 
(lower and upper), and initial generation interval parameters (shape and scale) 
and gets them from the shell script "SimulationSIFull.sh".


### SimulationSIFull.sh

This is a shell script that runs SimQsubSIFull.qsub which runs 
PerformSimulationSIFull.R with 9 different scenarios varying the 9 different input 
parameters.

### ClusteringExample.R

This program simulates a mini outbreak and pulls out two example cases - one with a top
cluster of infectors and one without a top cluster of infectors to create the clustering
example figure.

### Paper2Results.R

This program reads in all of the results from the the simulation scenarios and then 
calculates all results including values in the text, tables and figures of the manuscript.


***


## Massachusetts Analysis Programs

The following programs were used to prepare, analyze, and evaluate applying the 
naive Bayes transmission method to TB surveillance data from the Massachusetts 
Department of Public Health (DPH). The dataset used is not publically available, but the 
next section describes the simulated dataset with the same form that is included
which can be used to create similar results.


### MassPrep.R

This program reads in the demographic and contact, creates 
cleaned individual-level and pair-level datasets to be used for analysis.


### MassAnalysis.R

This program calculates relative transmission probabilities using naive Bayes, 
training with contact data. It then estimates the serial interval both including
all cases and excluding one-month co-prevelant cases with 95% confidence intervals 
and the monthly and overall average reproductive numbers also with confidence 
intervals.


### MassQsub.qsub

A qsub program to run MassAnalysis.R on a shared computing cluster because it takes
5 hours to run with the ~900 individual Massachusetts dataset and 1 hour to run
with the simulated Massachusetts-like dataset.


### MassEvaluate.R

This program reads in all of the results from the Massachusetts analysis and then 
calculates all results including values in the text, tables and figures of the manuscript.


*** 

## Simulated Massachusetts-like Data

Because the Massachusetts dataset could not be shared publically, we created a
simulated dataset with the same form as the real dataset. This simulated dataset
is called "orderedMassSim.rds" and it can be used to run MassAnalysis.R starting
at line 64 as indicated in the comments. The results produced by this program can then
be used in MassEvaluate.R in order to replicate the form of the results, tables, and
figures in the manuscript though the actual values will be different.

In order to allow for replication of the Massachusetts analysis, the results of 
selected analyses with this simulated dataset are located at the end of the 
Supplementary Material of the manuscript.
