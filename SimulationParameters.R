#Sarah V. Leavitt
#Boston University Dissertation
#Paper 2

################################################################################
# This program identifies the thresholds to use for the generation interval
# simulations.
################################################################################

#rm(list = ls())
options(scipen=999)
setwd("~/Boston University/Dissertation/dissertation_code")

library(TransPhylo)
library(phangorn)
library(reshape2)
library(dplyr)
library(tidyr)
library(lubridate)
library(devtools)
library(ggplot2)
library(gtools)
library(devtools)
library(roxygen2)

source("SimOutbreak.R")
source("SimulateOutbreakS.R")
source("SimCovariates.R")
#load("../Datasets/SimulationParameters.RData")


#Constant parameters
neg <- 0.25
pi <- 1
off.p <- 0.5
w.shift <- 0
multOutbreaks <- FALSE
rootseq <- NULL
length <- 3000


#### Function to find the thresholds ####

testThreshold <- function(off.r, w.scale, w.shape, mutationR, sampleSize, nSim){
  
  allT <- NULL
  for(j in 1:nSim){
    rate <- mutationR / length
    time <- (log(sampleSize, off.r) + 1)*qgamma(0.5, shape = w.shape, scale = w.scale)
    ws.shape <- w.shape
    ws.scale <- w.scale
    
    #Simulate outbreak  
    obk <- simOutbreak(neg = neg, pi = pi, off.r = off.r, off.p = off.p,
                       w.scale = w.scale, w.shape = w.shape, w.shift = w.shift,
                       ws.scale = ws.scale, ws.shape = ws.shape, ws.shift = w.shift,
                       sampleSize = sampleSize, time = time,
                       multOutbreaks = multOutbreaks, length = length,
                       rate = rate)
    indData <- obk[[1]]
    pairData <- obk[[2]]
    nrow(indData)
    
    #Finding the lowest threshold where the proportion of true links less than t is at least 65%
    i <- 1
    nT <- 0
    while(nT < 0.65){
      sensTab <- prop.table(table(pairData$transmission, pairData$snpDist < i), 1)
      if(dim(sensTab)[2] == 2){
        if(sensTab[2, 2] > 0.65) {
          lowerT <- bind_cols("metric" = "lowerT", "t" = i, "value" = sensTab[2, 2])
        }
        nT <- sensTab[2, 2]
      }
      i <- i + 1 
    }
    
    #Finding the lowest threshold where the proportion of true links greater than t is at most 1%
    i <- 1
    nT <- 1
    while(nT > 0.01){
      sensTab <- prop.table(table(pairData$transmission, pairData$snpDist < i), 1)
      if(sensTab[2, 1] < 0.01) {
        upperT <- bind_cols("metric" = "upperT", "t" = i, "value" = sensTab[2, 1])
      }
      nT <- sensTab[2, 1]
      i <- i + 1
    }
  timeRow <- bind_cols("metric" = "time", t = time, "value" = NA) 
  allT <- allT %>% bind_rows(lowerT, upperT, timeRow)
  }
  return(allT)
}

set.seed(103020)

#Baseline
s1 <- testThreshold(off.r = 1.5, w.shape = 2.25, w.scale = 0.0122,
              sampleSize = 300, mutationR = 25, nSim = 15)
tapply(s1$t, s1$metric, summary)
tapply(s1$value, s1$metric, summary) 


#Reproductive number
s2 <- testThreshold(off.r = 1.2, w.shape = 2.25, w.scale = 0.0122,
              sampleSize = 300, mutationR = 25, nSim = 15)
tapply(s2$t, s2$metric, summary) 
tapply(s2$value, s2$metric, summary) 

s3 <- testThreshold(off.r = 2.0, w.shape = 2.25, w.scale = 0.0122,
                    sampleSize = 300, mutationR = 25, nSim = 20)
tapply(s3$t, s3$metric, summary)
tapply(s3$value, s3$metric, summary) 


#Sample size
s4 <- testThreshold(off.r = 1.5, w.shape = 2.25, w.scale = 0.0122,
                    sampleSize = 100, mutationR = 25, nSim = 15)
tapply(s4$t, s4$metric, summary)
tapply(s4$value, s4$metric, summary) 

s5 <- testThreshold(off.r = 1.5, w.shape = 2.25, w.scale = 0.0122,
                    sampleSize = 500, mutationR = 25, nSim = 15)
tapply(s5$t, s5$metric, summary)
tapply(s5$value, s5$metric, summary) 


#Mutation rate
s6 <- testThreshold(off.r = 1.5, w.shape = 2.25, w.scale = 0.0122,
                    sampleSize = 300, mutationR = 5, nSim = 15)
tapply(s6$t, s6$metric, summary)
tapply(s6$value, s6$metric, summary) 

s7 <- testThreshold(off.r = 1.5, w.shape = 2.25, w.scale = 0.0122,
                    sampleSize = 300, mutationR = 50, nSim = 15)
tapply(s7$t, s7$metric, summary)
tapply(s7$value, s7$metric, summary) 


#Generation interval CV
s8 <- testThreshold(off.r = 1.5, w.shape = 1, w.scale = 0.0274,
                    sampleSize = 300, mutationR = 25, nSim = 15)
tapply(s8$t, s8$metric, summary)
tapply(s8$value, s8$metric, summary)

s9 <- testThreshold(off.r = 1.5, w.shape = 4, w.scale = 0.00685,
                    sampleSize = 300, mutationR = 25, nSim = 15)
tapply(s9$t, s9$metric, summary)
tapply(s9$value, s9$metric, summary) 


#Generation interval mean
s10 <- testThreshold(off.r = 1.5, w.shape = 2.25, w.scale = 0.00609,
                    sampleSize = 300, mutationR = 25, nSim = 15)
tapply(s10$t, s10$metric, summary)
tapply(s10$value, s10$metric, summary)

s11 <- testThreshold(off.r = 1.5, w.shape = 2.25, w.scale = 0.0365,
                    sampleSize = 300, mutationR = 25, nSim = 15)
tapply(s11$t, s11$metric, summary)
tapply(s11$value, s11$metric, summary) 


#Saving workspace of results
save(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11,
     file = "../Datasets/SimulationParameters.RData")



#### Exploring relationship between n, R0, SI, and time ####

off.r = 2
w.shape = 2
w.scale = 1
shift = 0
sampleSize = NA
time = 6
n <- NULL
for(i in 1:1000){
  simu <- simulateOutbreakS(neg = neg, pi = pi, off.r = off.r, off.p = off.p,
                            w.shape = w.shape, w.scale = w.scale, shift = shift,
                            ws.shape = ws.shape, ws.scale = ws.scale,
                            nSampled = sampleSize,
                            dateStartOutbreak = 1980, dateT = 1980 + time)
  indData <- as.data.frame(extractTTree(simu)$ttree)
  n <- c(n, nrow(indData))
}
mean(n)
m <- median(rgamma(1000, shape = w.shape, scale = w.scale))
off.r ^ (time / m)


