#Sarah V. Leavitt
#Boston University Dissertation
#Paper 2

################################################################################
# This program estimates transmission probabilities and serial interval, and 
# reproductive number for the Mass DPH data using the NB transmission method.
# The data are cleaned in MassPrep.R.

# NOTE: Program takes around 3 hours to run completely. The best way is to run
# the first section and then the next two sections as local jobs in parallel.
################################################################################


setwd("~/Boston University/Dissertation")
rm(list = ls())


######################### Section 1: Estimating Probabilities ############################

library(dplyr)
library(tidyr)
library(devtools)
load_all("../nbTransmission")

#Reading in cleaned datasets from MassPrep.R
set.seed(103020)
massInd <- readRDS("Datasets/MassInd.rds")
massPair <- readRDS("Datasets/MassPair.rds")

#How many pairs are different lineages (52%)
sum(massPair$Lineage == "Different", na.rm = TRUE)
sum(massPair$Lineage == "Different", na.rm = TRUE) / nrow(massPair)

#Creating an ordered dataset that also removes pairs with different lineages
orderedMass <- (massPair
                %>% filter(CombinedDiff >= 0, Lineage == "Same" | is.na(Lineage))
                %>% select(EdgeID, StudyID.1, StudyID.2, ContactGroup, Lineage.1, Lineage.2,
                           CombinedDt.1, CombinedDt.2, RecentArrival.1, RecentArrival.2,
                           County, Sex, Age, Spoligotype, MIRUDiff, MIRUDiffG, GENType,
                           PCRType, Lineage, CountryOfBirth, Smear, SharedResG, AnyImmunoSup,
                           TimeCat, CombinedDiff, CombinedDiffY, ContactTrain)
                #Creating a gold standard based on the GenType
                %>% mutate(miruLink = ifelse(GENType == "Same" & County == "Same", TRUE,
                                             ifelse(MIRUDiffG == "4+", FALSE, NA)))
)

print(table(orderedMass$ContactTrain, useNA = "always"))
print(prop.table(table(orderedMass$ContactTrain, useNA = "always")))

#Looking at all contact pairs
contactPairs <- (orderedMass
                 %>% filter(ContactGroup == TRUE)
                 %>% select(EdgeID, RecentArrival.1, RecentArrival.2, CombinedDt.1, CombinedDt.2,
                            CombinedDiff, Lineage.1, Lineage.2, Spoligotype, GENType, PCRType,
                            MIRUDiffG, SharedResG, ContactTrain)
)



#### Estimating Probabilities ####

#Estimating the probabilities with time difference
covariates <- c("Sex", "Age", "CountryOfBirth", "County", "Smear", "AnyImmunoSup",
                "SharedResG", "GENType", "TimeCat")

resMass <- nbProbabilities(orderedPair = orderedMass, indIDVar = "StudyID", pairIDVar = "EdgeID",
                           goldStdVar = "ContactTrain", covariates = covariates,
                           label = "ContactTime", l = 0.5, n = 10, m = 1, nReps = 10)

resMassCov <- (orderedMass
               %>% full_join(resMass$probabilities, by = "EdgeID")
               #Setting probabilities to 0 if infectee was a recent immigrant but not if it was a training link
               %>% mutate(pScaledI = ifelse(!is.na(RecentArrival.2) & 
                                              RecentArrival.2 == TRUE &
                                              (is.na(ContactTrain) | ContactTrain != TRUE), 0, pScaled),
                          pAvgI = ifelse(!is.na(RecentArrival.2) & 
                                           RecentArrival.2 == TRUE &
                                           (is.na(ContactTrain) | ContactTrain != TRUE), 0, pScaled))
)
resMassCoeff <- resMass$estimates

print("Finished estimating probabilities with time difference")


#Estimating probabilities without time difference
covariates2 <- c("Sex", "Age", "CountryOfBirth", "County", "Smear", "AnyImmunoSup",
                 "SharedResG", "GENType")

resMass2 <- nbProbabilities(orderedPair = orderedMass, indIDVar = "StudyID", pairID = "EdgeID",
                            goldStdVar = "ContactTrain", covariates = covariates2,
                            label = "ContactNoTime", l = 0.5, n = 10, m = 1, nReps = 20)

resMassCov2 <- (orderedMass
                %>% full_join(resMass2$probabilities, by = c("EdgeID"))
                #Setting probabilities to 0 if infectee was a recent immigrant but not if it was a training link
                %>% mutate(pScaledI = ifelse(!is.na(RecentArrival.2) & 
                                               RecentArrival.2 == TRUE &
                                               (is.na(ContactTrain) | ContactTrain != TRUE), 0, pScaled),
                           pAvgI = ifelse(!is.na(RecentArrival.2) & 
                                            RecentArrival.2 == TRUE &
                                            (is.na(ContactTrain) | ContactTrain != TRUE), 0, pScaled))
)
resMassCoeff2 <- resMass2$estimates

print("Finished estimating probabilities without time difference")


#Saving the results
saveRDS(resMassCov, "Datasets/MassResults.rds")
saveRDS(resMassCoeff, "Datasets/MassORs.rds")
saveRDS(resMassCov2, "Datasets/MassResults_NoTime.rds")
saveRDS(resMassCoeff, "Datasets/MassORsNoTime.rds")




###################### Section 2: Serial Interval ########################

#Loading libraries (repeating to run in parallel)
library(dplyr)
library(tidyr)
library(devtools)
load_all("../nbTransmission")


set.seed(103020)
resMassCov2 <- readRDS("../Datasets/MassResults_NoTime.rds")

siHC <- estimateSI(df = resMassCov2, indIDVar = "StudyID",
                  timeDiffVar = "CombinedDiffY", pVar = "pScaledI",
                  clustMethod = "hc_absolute", cutoffs = seq(0.025, 0.25, 0.025),
                  initialPars = c(1.2, 2), shift = 0, bootSamples = 2)
siHC$label <- "HC: No exclusions"

siHC12 <- estimateSI(df = resMassCov2, indIDVar = "StudyID",
                   timeDiffVar = "CombinedDiffY", pVar = "pScaledI",
                   clustMethod = "hc_absolute", cutoffs = seq(0.025, 0.25, 0.025),
                   initialPars = c(1.2, 2), shift = 1/12, bootSamples = 2)
siHC12$label <- "HC: Excluding 1-month co-prevalent cases"


siKD <- estimateSI(df = resMassCov2, indIDVar = "StudyID",
                   timeDiffVar = "CombinedDiffY", pVar = "pScaledI",
                   clustMethod = "kd", cutoffs = seq(0.01, 0.1, 0.01),
                   initialPars = c(1.2, 2), shift = 0, bootSamples = 2)
siKD$label <- "KD: No exclusions"

siKD12 <- estimateSI(df = resMassCov2, indIDVar = "StudyID",
                    timeDiffVar = "CombinedDiffY", pVar = "pScaledI",
                    clustMethod = "kd", cutoffs = seq(0.01, 0.1, 0.01),
                    initialPars = c(1.2, 2), shift = 1/12, bootSamples = 2)
siKD12$label <- "KD: Excluding 1-month co-prevalent cases"


siAll <- bind_rows(siHC, siHC12, siKD, siKD12)

#Saving the serial interval dataset
saveRDS(siAll, "../Datasets/MassSI.rds")




####################### Section 3: Reproductive Number ###########################

#Loading libraries (repeating to run in parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(devtools)
load_all("../nbTransmission")

set.seed(103020)
resMassCov <- readRDS("../Datasets/MassResults.rds")

#Initially calculating reproductive number to decide cut points
rInitial <- estimateR(df = resMassCov, dateVar = "CombinedDt", indIDVar = "StudyID",
                      pVar = "pScaledI", timeFrame = "months")
rInitial$RtAvgDf
rt <- rInitial$RtDf

#Cutting the outbreak
totalTime <- max(rt$timeRank) - min(rt$timeRank)
monthCut1 <- ceiling(0.1 * totalTime)
monthCut2 <- ceiling(0.8 * totalTime)

#Plotting where to cut
ggplot(data = rt) +
  geom_line(aes(x = timeRank, y = Rt)) +
  scale_y_continuous(name = "Rt") +
  geom_vline(aes(xintercept = monthCut1), linetype = 2, size = 0.7, col = "blue") +
  geom_vline(aes(xintercept = monthCut2), linetype = 2, size = 0.7, col = "blue")


#Calculating the reproductive number accounting for importation
rFinal <- estimateR(resMassCov, dateVar = "CombinedDt", indIDVar = "StudyID",
                    pVar = "pScaledI", timeFrame = "months",
                    rangeForAvg = c(monthCut1, monthCut2),
                    bootSamples = 2, alpha = 0.05)

rFinal$RiDf$label <- "Accounting for Importation"
rFinal$RtDf$label <- "Accounting for Importation"
rFinal$RtAvgDf$label <- "Accounting for Importation"


#Saving the confidence interval datasets
saveRDS(rFinal$RiDf, "../Datasets/MassRi.rds")
saveRDS(rFinal$RtDf, "../Datasets/MassRtCI.rds")
saveRDS(rFinal$RtAvgDf, "../Datasets/MassRtAvgCI.rds")



