#!/bin/bash

#Running the different scenarios

#Arguments are:
	#sample size
	#reproductive number
	#generation interval shape
	#generation interval scale
	#mutation rate
	#lower SNP threshold (link)
	#upper SNP threshold (nonlink)
	#initial shape parameter
	#initial scale parameter
	
qsub -N SI_Baseline SimQsubSIFull.qsub 300 1.5 2.25 0.0122 25 3 7 3 4

qsub -N SI_LowR0 SimQsubSIFull.qsub 300 1.2 2.25 0.0122 25 3 7 3 4
qsub -N SI_HighR0 SimQsubSIFull.qsub 300 2.0 2.25 0.0122 25 2 6 3 4

qsub -N SI_LowGIV SimQsubSIFull.qsub 300 1.5 4 0.00685 25 3 7 3 4
qsub -N SI_HighGIV SimQsubSIFull.qsub 300 1.5 1 0.0274 25 2 7 3 4

qsub -N SI_LowGIM SimQsubSIFull.qsub 300 1.5 2.25 0.00609 25 2 5 3 3
qsub -N SI_HighGIM SimQsubSIFull.qsub 300 1.5 2.25 0.0365 25 6 13 3 10

qsub -N SI_LowN SimQsubSIFull.qsub 100 1.5 2.25 0.0122 25 3 7 3 4
qsub -N SI_HighN SimQsubSIFull.qsub 500 1.5 2.25 0.0122 25 3 7 3 4

qsub -N SI_LowMR SimQsubSIFull.qsub 300 1.5 2.25 0.0122 5 1 4 3 4
qsub -N SI_HighMR SimQsubSIFull.qsub 300 1.5 2.25 0.0122 50 4 11 3 4

