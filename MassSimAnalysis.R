
#setwd("~/Boston University/Dissertation/nbPaper2")
setwd("/project/sv-thesis/nbPaper2")
rm(list = ls())

library(dplyr)
library(tidyr)
library(devtools)
load_all("../nbTransmission")
set.seed(103020)


orderedMass <- readRDS("orderedMassSim.rds")

#Estimating the probabilities with time difference
covariates <- c("Sex", "Age", "CountryOfBirth", "County", "Smear", "AnyImmunoSup",
                "SharedResG", "GENType", "TimeCat")

resMass <- nbProbabilities(orderedPair = orderedMass, indIDVar = "StudyID", pairIDVar = "EdgeID",
                           goldStdVar = "ContactTrain", covariates = covariates,
                           label = "ContactTime", l = 0.5, n = 10, m = 1, nReps = 50,
                           progressBar = FALSE)

resMassCov <- (orderedMass
               %>% full_join(resMass$probabilities, by = "EdgeID")
               #Setting probabilities to 0 if infectee was a recent immigrant but not if it was a training link
               %>% mutate(pScaledI2 = ifelse(!is.na(RecentArrival2.2) &
                                               RecentArrival2.2 == TRUE &
                                               (is.na(ContactTrain) | ContactTrain != TRUE), 0, pScaled))
)
resMassCoeff <- resMass$estimates

print("Finished estimating probabilities with time difference")


#Estimating probabilities without time difference
covariates2 <- c("Sex", "Age", "CountryOfBirth", "County", "Smear", "AnyImmunoSup",
                 "SharedResG", "GENType")

resMass2 <- nbProbabilities(orderedPair = orderedMass, indIDVar = "StudyID", pairID = "EdgeID",
                            goldStdVar = "ContactTrain", covariates = covariates2,
                            label = "ContactNoTime", l = 0.5, n = 10, m = 1, nReps = 50,
                            progressBar = FALSE)

resMassCov2 <- (orderedMass
                %>% full_join(resMass2$probabilities, by = c("EdgeID"))
                #Setting probabilities to 0 if infectee was a recent immigrant but not if it was a training link
                %>% mutate(pScaledI2 = ifelse(!is.na(RecentArrival2.2) &
                                                RecentArrival2.2 == TRUE &
                                                (is.na(ContactTrain) | ContactTrain != TRUE), 0, pScaled))
)
resMassCoeff2 <- resMass2$estimates

print("Finished estimating probabilities without time difference")




###################### Serial Interval ########################


siHC1 <- estimateSI(df = resMassCov2, indIDVar = "StudyID",
                    timeDiffVar = "CombinedDiffYM", pVar = "pScaledI2",
                    clustMethod = "hc_absolute", cutoffs = seq(0.025, 0.25, 0.025),
                    initialPars = c(1.2, 2), shift = 0, bootSamples = 1000, progressBar = FALSE)
siHC1$label <- "HC: Excluding 1-month co-prevalent cases"
print("HC: Excluding 1-month co-prevalent cases")

siHC2 <- estimateSI(df = resMassCov2, indIDVar = "StudyID",
                    timeDiffVar = "CombinedDiffYM", pVar = "pScaledI2",
                    clustMethod = "hc_absolute", cutoffs = seq(0.025, 0.25, 0.025),
                    initialPars = c(1.2, 2), shift = 1/12, bootSamples = 1000, progressBar = FALSE)
siHC2$label <- "HC: Excluding 2-month co-prevalent cases"
print("HC: Excluding 2-month co-prevalent cases")

siHC3 <- estimateSI(df = resMassCov2, indIDVar = "StudyID",
                    timeDiffVar = "CombinedDiffYM", pVar = "pScaledI2",
                    clustMethod = "hc_absolute", cutoffs = seq(0.025, 0.25, 0.025),
                    initialPars = c(1.2, 2), shift = 2/12, bootSamples = 1000, progressBar = FALSE)
siHC3$label <- "HC: Excluding 3-month co-prevalent cases"
print("HC: Excluding 3-month co-prevalent cases")


siKD1 <- estimateSI(df = resMassCov2, indIDVar = "StudyID",
                    timeDiffVar = "CombinedDiffYM", pVar = "pScaledI2",
                    clustMethod = "kd", cutoffs = seq(0.01, 0.1, 0.01),
                    initialPars = c(1.2, 2), shift = 0, bootSamples = 1000, progressBar = FALSE)
siKD1$label <- "KD: Excluding 1-month co-prevalent cases"
print("KD: Excluding 1-month co-prevalent cases")

siKD2 <- estimateSI(df = resMassCov2, indIDVar = "StudyID",
                    timeDiffVar = "CombinedDiffYM", pVar = "pScaledI2",
                    clustMethod = "kd", cutoffs = seq(0.01, 0.1, 0.01),
                    initialPars = c(1.2, 2), shift = 1/12, bootSamples = 1000, progressBar = FALSE)
siKD2$label <- "KD: Excluding 2-month co-prevalent cases"
print("KD: Excluding 2-month co-prevalent cases")

siKD3 <- estimateSI(df = resMassCov2, indIDVar = "StudyID",
                    timeDiffVar = "CombinedDiffYM", pVar = "pScaledI2",
                    clustMethod = "kd", cutoffs = seq(0.01, 0.1, 0.01),
                    initialPars = c(1.2, 2), shift = 2/12, bootSamples = 1000, progressBar = FALSE)
siKD3$label <- "KD: Excluding 3-month co-prevalent cases"
print("KD: Excluding 3-month co-prevalent cases")


siAll <- bind_rows(siHC1, siHC2, siHC3, siKD1, siKD2, siKD3)

#Saving the serial interval dataset
saveRDS(siAll, "../Datasets/MassSimSI.rds")




####################### Reproductive Number ###########################

#Initially calculating reproductive number to decide cut points
rInitial <- estimateR(df = resMassCov, dateVar = "CombinedDt", indIDVar = "StudyID",
                      pVar = "pScaledI2", timeFrame = "months")
rInitial$RtAvgDf
rt <- rInitial$RtDf

#Cutting the outbreak
totalTime <- max(rt$timeRank) - min(rt$timeRank)
monthCut1 <- ceiling(0.1 * totalTime)
monthCut2 <- ceiling(0.8 * totalTime)


#Calculating the reproductive number using 1 year definition for recent immigration
rFinal2 <- estimateR(resMassCov, dateVar = "CombinedDt", indIDVar = "StudyID",
                     pVar = "pScaledI2", timeFrame = "months",
                     rangeForAvg = c(monthCut1, monthCut2),
                     bootSamples = 1000, alpha = 0.05, progressBar = FALSE)

rFinal2$RiDf$label <- "Recent Arrival = 2 Years"
rFinal2$RtDf$label <- "Recent Arrival = 2 Years"
rFinal2$RtAvgDf$label <- "Recent Arrival = 2 Years"

print("Recent Arrival = 2 Year2")

RiData <- rFinal2$RiDf
RtData <- rFinal2$RtDf
RtAvg <- rFinal2$RtAvgDf


#Saving the confidence interval datasets
saveRDS(RiData, "../Datasets/MassSimRi.rds")
saveRDS(RtData, "../Datasets/MassSimRtCI.rds")
saveRDS(RtAvg, "../Datasets/MassSimRtAvgCI.rds")


