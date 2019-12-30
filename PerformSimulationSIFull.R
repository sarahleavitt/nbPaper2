#Sarah V. Leavitt
#Boston University Dissertation
#Paper 2

#################################################################################
# This program runs the full set of generation simulations
# by calling SimRunSIFull.R
#################################################################################

rm(list = ls())
options(scipen=999)

library(TransPhylo)
library(phangorn)
library(reshape2)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)
library(pROC)
library(caret)
library(gtools)
library(devtools)
library(roxygen2)


#### Batch Mode ####

setwd("/project/sv-thesis/dissertation_code/")
#Getting sample size from the arguements
argv <- commandArgs(trailingOnly = TRUE)
#Finding the task number for the run
iTask <- as.numeric(Sys.getenv("SGE_TASK_ID"))
#The number of simulations per split
nSim <- 50

## Parameters to change ##
sampleSize <- as.numeric(argv[1])
off.r <- as.numeric(argv[2])
w.shape <- as.numeric(argv[3])
w.scale <- as.numeric(argv[4])
#snps/genome/year
mutationR <- as.numeric(argv[5])
#Method parameters
lowerT <- as.numeric(argv[6])
upperT <- as.numeric(argv[7])
iPar1 <- as.numeric(argv[8])
iPar2 <- as.numeric(argv[9])


#### Interactive Mode ####

# setwd("~/Boston University/Dissertation/dissertation_code")
# iTask <- 1
# nSim <- 2
# 
# ## Parameters to change ##
# sampleSize <- 300
# off.r <- 1.5
# w.shape <- 2.25
# w.scale <- 0.0122
# #snps/genome/year
# mutationR <- 25
# #Method parameters
# lowerT <- 3
# upperT <- 7
# iPar1 <- 3
# iPar2 <- 4

load_all("../nbTransmission")
source("SimOutbreak.R")
source("SimulateOutbreakS.R")
source("SimCovariates.R")
source("SimEvaluate.R")
source("SimRunSIFull.R")




################## Setting the oubreak parameters ################

#Setting label
if(off.r < 1.3){
  label <- "LowR"
} else if(off.r > 1.7){
  label <- "HighR"
} else if(w.shape > 3){
  label <- "LowGIV"
} else if(w.shape < 2){
  label <- "HighGIV"
} else if(w.scale < 0.01){
  label <- "LowGIM"
} else if(w.scale > 0.02){
  label <- "HighGIM"
} else if(sampleSize < 300){
  label <- "LowN"
} else if(sampleSize > 300){
  label <- "HighN"
} else if(mutationR < 25){
  label <- "LowMR"
} else if(mutationR > 25){
  label <- "HighMR"
} else{
  label <- "Baseline"
}

## Constant parameters ##
neg <- 0.25
pi <- 1
off.p <- 0.5
w.shift <- 0
ws.shape <- w.shape
ws.scale <- w.scale
ws.shift <- w.shift
multOutbreaks <- FALSE
rootseq <- NULL
length <- 3000
observationDate <- "infectionDate"


## Derived parameters ##

#This gives the equivalent of mutationRate mutations/genome/year
rate <- mutationR / length
time <- (log(sampleSize, off.r) + 1)*qgamma(0.5, shape = w.shape, scale = w.scale)
#time <- 10 * qgamma(0.5, shape = w.shape, scale = w.scale)
initialPars <- c(iPar1, iPar2)
if(w.shape*w.scale < 0.25){
  dateVar <- "observationDiff"
}else{
  dateVar <- "observationDiffY"
}
pTraining <- 0.6


cat(paste0("Beginning '", label, "' run with sample size = ", sampleSize,
           ", R0 = ", off.r, ", serial interval = (", w.shape, ", ", w.scale,
           "), mutation rate = ", mutationR, ", \ntime = ", round(time, 2),
           " years, SNP thresholds = ", lowerT, "/", upperT, ", pTraining = ",
           pTraining, ", and initial parameters = (", iPar1, ", ", iPar2,
           ") using ", dateVar, ".\n"))


#Parameters used to manually run functions
covariates <- c("Y1", "Y2", "Y3", "Y4")
goldStandard <- "snpClose"
observationDate <- "infectionDate"
truth <- "transmission"
pVar <- "pScaled"
labelVar <- "label"
nReps <- 10
indIDVar = "individualID"
timeDiffVar = "infectionDiffY"
pVar = "pScaled"
clustMethod = "hc_absolute"
cutoffs = 0.05
shift = 0
epsilon = 0.0001
bootSamples = 0
alpha = 0.05



####################### Running the Simulations ######################

si <- NULL
performance <- NULL
idenInfo <- NULL

#Running the simuliaton for the outbreaks
for (iteration in 1:nSim){
  
  #Setting the seed for each run so I can re-run specific iterations with errors
  set.seed(iTask * 1000 + iteration)
  
  tryCatch({
    
    #Running the simulation
    res <- SimRunSI(label, observationDate)
    
    runID <- paste(iTask, iteration, sampleSize, sep = "_")
    sTemp <- res[[1]] %>% mutate(runID = runID)
    pTemp <- res[[2]] %>% mutate(runID = runID)
    
    si <- bind_rows(si, sTemp)
    performance <- bind_rows(performance, pTemp)
    #Printing message that the run is finished
    print(paste0("Completed run ", iteration, " (seed = ", iTask * 1000 + iteration, ")"))
    
  }, error = function(e){
    #Printing message that the run is finished
    print(paste0("Error in run ", iteration, " (seed = ", iTask * 1000 + iteration, ")"))
    cat("ERROR :", conditionMessage(e), "\n")})
}

#Saving dataframes with summary of results for all of the runs
saveRDS(si, file=paste0("../Simulation_Results/si", label, "_", iTask, ".rds"))
saveRDS(performance, file=paste0("../Simulation_Results/performance", label, "_", iTask, ".rds"))


