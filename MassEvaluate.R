#Sarah V. Leavitt
#Boston University Dissertation
#Paper 2

################################################################################
# This program creates the figures and tables to analyze the Mass DPH data
# The data are cleaned in MassPrep.R and analyzed in MassAnalysis.R
################################################################################

setwd("~/Boston University/Dissertation")
rm(list = ls())
options(scipen = 999)

library(dplyr)
library(tidyr)
library(lubridate)
library(tableone)
library(devtools)
library(ggplot2)
library(ggforce)
load_all("../nbTransmission")


#Reading in cleaned datasets from MassPrep.R and results from MassAnalysis.R
massInd <- readRDS("Datasets/MassInd.rds")
massPair <- readRDS("Datasets/MassPair.rds")
resMassCov <- readRDS("Datasets/MassResults.rds")
resMassCov2 <- readRDS("Datasets/MassResults_NoTime.rds")
siAll <- readRDS("Datasets/MassSI.rds")
RiData <- readRDS("Datasets/MassRi.rds")
RtData <- readRDS("Datasets/MassRtCI.rds")
RtAvg <- readRDS("Datasets/MassRtAvgCI.rds")

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



################## Covariate Tables #####################

massInd <- massInd %>% replace_na(list(HaveContInv = "No"))

## Individual Level ##
indCat <- c("Sex", "Age", "USBorn", "RecentArrival", "Smear", "AnyImmunoSup",
            "ISUSRIF", "ISUSINH",  "ISUSPZA", "ISUSEMB", "ISUSSM", "ISUSETH",
            "County", "Lineage", "HaveContact", "HaveContInv")

covarInd <- CreateTableOne(vars = indCat, factorVars = indCat, data = massInd)
covarInd <- as.data.frame(print(covarInd, showAllLevels = TRUE, missing = TRUE))

sum(is.na(massInd$GENType))
sum(is.na(massInd$Spoligotype))
sum(is.na(massInd$MIRUComb))


## Pair Level ##
pairCat <- c("Sex", "Age", "CountryOfBirth", "Smear", "AnyImmunoSup",
             "SharedResG", "County", "GENType", "TimeCat", "ContactTrain")

covarPair <- CreateTableOne(vars = pairCat, factorVars = pairCat, data = orderedMass)
covarPair <- as.data.frame(print(covarPair, showAllLevels = TRUE))

#Stratified by contact group
covarPairC <- CreateTableOne(vars = pairCat, factorVars = pairCat,
                            data = orderedMass, strata = "ContactTrain", test = FALSE)
covarPairC <- as.data.frame(print(covarPairC, showAllLevels = TRUE))

covarPairAll <- cbind.data.frame(covarPair, covarPairC)

table(orderedMass$Lineage, useNA = "always")
prop.table(table(orderedMass$Lineage, useNA = "always"))



################### Assessing Probabilities ###################

#One possible clustering method and cutoff
resMassCov2C <- clusterInfectors(df = resMassCov2, indIDVar = "StudyID", pVar = "pScaledI",
                                 clustMethod = "hc_absolute", cutoff = 0.1)

topClust <- resMassCov2C %>% filter(cluster == 1)
length(unique(topClust$StudyID.2))
length(unique(topClust$StudyID.2)) / length(unique(resMassCov2C$StudyID.2))


#### Figure: Plot of Probabilities Colored by Cluster ####
resMassCov2C <- resMassCov2C %>% mutate(clusterC = ifelse(cluster == 1, "Top Cluster",
                                                         "Bottom Cluster"))
ggplot(data = resMassCov2C) +
  geom_histogram(aes(x = pScaledI, fill = clusterC),
                 binwidth = 0.1, position = "dodge") +
  scale_y_continuous(name = "Number of Case Pairs") +
  scale_x_continuous(name = "Relative Transmission Probability") +
  scale_fill_grey(start = 0.6, end = 0.3) +
  facet_zoom(ylim = c(0, 300)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ggsave(file = "Figures/MassProbs.png",
         width = 8, height = 6, units = "in", dpi = 300)
  

#Exploring probabilities of contact pairs and high probability pairs
contactPairs <- (resMassCov
                 %>% filter(ContactTrain == TRUE)
                 %>% select(EdgeID, StudyID.1, StudyID.2, CombinedDt.1, CombinedDt.2,
                            pScaledI, pAvg, nSamples, Spoligotype, MIRUDiffG,
                            County, CountryOfBirth)
)

highProb <- resMassCov %>% filter(pScaledI > 0.3)
table(highProb$GENType, useNA = "ifany")
table(highProb$ContactTrain, useNA = "ifany")

pairCat2 <- c("Sex", "Age", "CountryOfBirth", "Smear", "AnyImmunoSup",
             "SharedResG", "County", "GENType", "TimeCat")

covarPair2 <- CreateTableOne(vars = pairCat2, factorVars = pairCat2, data = orderedMass)
covarPair2 <- as.data.frame(print(covarPair2, showAllLevels = TRUE))

covarPairH <- CreateTableOne(data = highProb, vars = pairCat2, factorVars = pairCat2)
covarPairH <- as.data.frame(print(covarPairH, showAllLevels = TRUE))
covarPairAll2 <- cbind.data.frame(covarPair2, covarPairH)




######################## Serial Interval #########################

#Function to take the SI results and format to a nice table
formatSITable <- function(siTable){
  siTable2 <- (siTable
               %>% mutate(npIncluded = paste0(nIndividuals, " (",
                                              100 * round(pCluster, 3), ")"),
                          mean = paste0(round(meanSI, 2), " (", round(meanCILB, 2),
                                        ", ", round(meanCIUB, 2), ")"),
                          median = paste0(round(medianSI, 2), " (", round(medianCILB, 2),
                                          ", ", round(medianCIUB, 2), ")"),
                          sd = paste0(round(sdSI, 2), " (", round(sdCILB, 2),
                                      ", ", round(sdCIUB, 2), ")"))
               %>% select(label, cutoff, npIncluded, mean, median, sd)
  )
  return(siTable2)
}

#Tables of pooled results for text
pooled <- formatSITable(siAll %>% filter(cutoff == "pooled"))

#### Supplementary Tables: Detailed Serial Interval Results ####
siHC <- formatSITable(siAll %>% filter(clustMethod == "hc_absolute"))
siKD <- formatSITable(siAll %>% filter(clustMethod == "kd"))


#Creating alternative label
siAll <- siAll %>% mutate(label2 = gsub("[A-Z]{2}: ", "", label))

## Creating long dataset ##
meanDf <- (siAll
           %>% select(label2, clustMethod, cutoff, shape, scale,
                      est = meanSI, cilb = meanCILB, ciub = meanCIUB)
           %>% mutate(Parameter = "Mean")
)
medianDf <- (siAll
             %>% select(label2, clustMethod, cutoff, shape, scale,
                        est = medianSI, cilb = medianCILB, ciub = medianCIUB)
             %>% mutate(Parameter = "Median")
)
sdDf <- (siAll
         %>% select(label2, clustMethod, cutoff, shape, scale,
                    est = sdSI, cilb = sdCILB, ciub = sdCIUB)
         %>% mutate(Parameter = "Standard Deviation")
)
siAllLong <- bind_rows(meanDf, medianDf, sdDf)

#Finding the width for the error bars
errorWidth <- (siAllLong
               %>% filter(cutoff != "pooled")
               %>% group_by(clustMethod)
               %>% summarize(range = max(as.numeric(cutoff)) - min(as.numeric(cutoff)),
                             width = range / 40)
               %>% select(-range)
)
siAllLong2 <- (siAllLong
               %>% full_join(errorWidth, by = "clustMethod")
               %>% mutate(clustMethod = ifelse(clustMethod == "hc_absolute",
                                               "Hiearchical Clustering",
                                               "Kernel Density Estimation"))
)

#### Figure: Plot of Serial Interval Estimates with CIs ####
ggplot(data = siAllLong2 %>% filter(cutoff != "pooled"),
       aes(x = as.numeric(cutoff), y = est, color = label2)) +
  geom_point() +
  geom_errorbar(aes(ymin = cilb, ymax = ciub, width = width)) +
  geom_hline(data = siAllLong2 %>% filter(cutoff == "pooled"),
             aes(yintercept = est, color = label2)) +
  geom_hline(data = siAllLong2 %>% filter(cutoff == "pooled"),
             aes(yintercept = cilb, color = label2), linetype = "dotted") +
  geom_hline(data = siAllLong2 %>% filter(cutoff == "pooled"),
             aes(yintercept = ciub, color = label2), linetype = "dotted") +
  facet_grid(Parameter~clustMethod, scales = "free") +
  scale_x_continuous(name = "Clustering Cutoff/Binwidth") +
  scale_y_continuous(name = "Estimate in Years") +
  theme_bw() +
  scale_color_grey(start = 0.6, end = 0.3) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                    size = 11),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                    size = 11)) +
  ggsave(file = "Figures/MassSI.png",
         width = 7, height = 7, units = "in", dpi = 300)



################### Reproductive Number ######################

#Average Rt accounting for importation
RtAvg

#Cutting the outbreak
totalTime <- max(RtData$timeRank) - min(RtData$timeRank)
monthCut1 <- ceiling(0.1 * totalTime)
monthCut2 <- ceiling(0.8 * totalTime)


#### Figure: Plot of Rt Estimates by Month with CIs ####
ggplot(data = RtData, aes(x = timeRank, y = Rt)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = ciLower, ymax = ciUpper), width = 0.7, color = "grey40") +
  scale_y_continuous(name = "Monthly Effective Reproductive Number") + 
  scale_x_continuous(name = "Year of Observation", breaks = seq(3, 89, 12),
                     labels = seq(2010, 2017, 1)) +
  geom_vline(aes(xintercept = monthCut1), linetype = "dotted", size = 0.7) +
  geom_vline(aes(xintercept = monthCut2), linetype = "dotted", size = 0.7) +
  geom_hline(data = RtAvg, aes(yintercept = RtAvg), size = 0.7) +
  theme_bw() +
  geom_hline(data = RtAvg, aes(yintercept = ciLower), linetype = "dashed",
             size = 0.5, color = "grey40") +
  geom_hline(data = RtAvg, aes(yintercept = ciUpper), linetype = "dashed",
             size = 0.5, color = "grey40") +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                    size = 12),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                    size = 12)) +
  ggsave(file = "Figures/MassRt.png",
         width = 8, height = 6, units = "in", dpi = 300)


