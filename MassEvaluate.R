#Sarah V. Leavitt
#Boston University Dissertation
#Paper 2

################################################################################
# This program creates the figures and tables to analyze the Mass DPH data
# The data are cleaned in MassPrep.R and analyzed in MassAnalysis.R
################################################################################

setwd("~/Boston University/Dissertation/nbPaper2")
#rm(list = ls())
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
massInd <- readRDS("../Datasets/MassInd.rds")
massPair <- readRDS("../Datasets/MassPair.rds")
resMassCov <- readRDS("../Datasets/MassResults.rds")
resMassCov2 <- readRDS("../Datasets/MassResults_NoTime.rds")
siAll <- readRDS("../Datasets/MassSI.rds")
RiData <- readRDS("../Datasets/MassRi.rds")
RtData <- readRDS("../Datasets/MassRtCI.rds")
RtAvg <- readRDS("../Datasets/MassRtAvgCI.rds")

#Creating an ordered dataset that also removes pairs with different lineages
orderedMass <- (massPair
                %>% filter(CombinedDiff >= 0, Lineage == "Same" | is.na(Lineage))
                %>% select(EdgeID, StudyID.1, StudyID.2, ContactGroup, Lineage.1, Lineage.2,
                           CombinedDt.1, CombinedDt.2, RecentArrival1.1, RecentArrival1.2,
                           RecentArrival2.1, RecentArrival2.2,
                           County, Sex, Age, Spoligotype, MIRUDiff, MIRUDiffG, GENType,
                           PCRType, Lineage, CountryOfBirth, Smear, SharedResG, AnyImmunoSup,
                           TimeCat, CombinedDiff, CombinedDiffY, ContactTrain)
                #Creating a gold standard based on the GenType
                %>% mutate(miruLink = ifelse(GENType == "Same" & County == "Same", TRUE,
                                             ifelse(MIRUDiffG == "4+", FALSE, NA)))
)


################## Creating Dataset for Maps #################

# countyPrev <- (massInd
#                %>% group_by(County)
#                %>% summarize(nCases = n())
#                %>% filter(!is.na(County))
# )
# 
# write.csv(countyPrev, "../MA_Map/countyPrev.csv", row.names = FALSE)



################## Covariate Tables #####################

massInd <- massInd %>% replace_na(list(HaveContInv = "No"))

## Individual Level ##
indCat <- c("Sex", "Age", "USBorn", "RecentArrival2", "Smear", "AnyImmunoSup",
            "ISUSRIF", "ISUSINH",  "ISUSPZA", "ISUSEMB", "ISUSSM", "ISUSETH",
            "County", "Lineage")

covarInd <- CreateTableOne(vars = indCat, factorVars = indCat, data = massInd)
covarInd <- as.data.frame(print(covarInd, showAllLevels = TRUE))

sum(is.na(massInd$GENType))
sum(is.na(massInd$Spoligotype))
sum(is.na(massInd$MIRUComb))

#Finding the amount of missing values for each variable
findMissingness <- function(data){
  
  #Calculating n(%) missing
  numMiss <- apply(data, 2, function(x)sum(is.na(x)))
  percMiss <- 100 * (numMiss / nrow(data))
  missData <- cbind.data.frame(numMiss, percMiss)
  #If a category is not completely missing but rounds to 100% that >99.99% missing is printed
  #Similarly, if it is mostly not missing but rounds to 0%, <0.01% is printed
  missData$percMissf <- ifelse(missData$numMiss != nrow(data) 
                               & round(missData$percMiss, 2) == 100, ">99.99",
                               ifelse(missData$numMiss != 0 & round(missData$percMiss, 2) == 0, "<0.01",
                                      sprintf("%.2f", round(missData$percMiss, 2))))
  missData$nPercMiss <- paste(missData$numMiss, " (", missData$percMissf, "%)", sep="")
  
  return(missData)
}
findMissingness(massInd[, indCat])


## Pair Level ##
pairCat <- c("Sex", "Age", "CountryOfBirth", "Smear", "AnyImmunoSup",
             "SharedResG", "County", "GENType", "TimeCat", "ContactTrain")

covarPair <- CreateTableOne(vars = pairCat, factorVars = pairCat, data = orderedMass)
covarPair <- as.data.frame(print(covarPair, showAllLevels = TRUE))

table(orderedMass$Lineage, useNA = "always")
prop.table(table(orderedMass$Lineage, useNA = "always"))



################### Assessing Probabilities ###################

# #One possible clustering method and cutoff
# resMassCov2C <- clusterInfectors(df = resMassCov2, indIDVar = "StudyID", pVar = "pScaledI2",
#                                  clustMethod = "hc_absolute", cutoff = 0.1)
# 
# topClust <- resMassCov2C %>% filter(cluster == 1)
# length(unique(topClust$StudyID.2))
# length(unique(topClust$StudyID.2)) / length(unique(resMassCov2C$StudyID.2))


#### Figure: Plot of Probabilities Colored by Cluster ####
resMassCov2C <- resMassCov2C %>% mutate(clusterC = ifelse(cluster == 1, "Top Cluster",
                                                         "Bottom Cluster"))
ggplot(data = resMassCov2C) +
  geom_histogram(aes(x = pScaledI2, fill = clusterC),
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
  ggsave(file = "../Figures/MassProbs.png",
         width = 8, height = 6, units = "in", dpi = 300)

## COLOR VERSION ##
ggplot(data = resMassCov2C) +
  geom_histogram(aes(x = pScaledI2, fill = clusterC),
                 binwidth = 0.1, position = "dodge") +
  scale_y_continuous(name = "Number of Case Pairs") +
  scale_x_continuous(name = "Relative Transmission Probability") +
  facet_zoom(ylim = c(0, 300)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ggsave(file = "../Figures/MassProbs_color.png",
         width = 8, height = 6, units = "in", dpi = 300)

## PRESENTATION VERSION ##
ggplot(data = resMassCov2C) +
  geom_histogram(aes(x = pScaledI2, fill = clusterC),
                 binwidth = 0.1, position = "dodge") +
  scale_y_continuous(name = "Number of Case Pairs") +
  scale_x_continuous(name = "Relative Transmission Probability") +
  facet_zoom(ylim = c(0, 300)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ggsave(file = "../Figures/MassProbs_pres.png",
         width = 8, height = 6, units = "in", dpi = 300)
  




######################## Serial Interval #########################

#Function to take the SI results and format to a nice table
formatSITable <- function(siTable){
  siTable2 <- (siTable
               %>% mutate(npIncluded = paste0(nIndividuals, " (",
                                              sprintf("%.1f", 100 * round(pCluster, 3)), ")"),
                          mean = paste0(sprintf("%.2f", round(meanSI, 2)), " (",
                                        sprintf("%.2f", round(meanCILB, 2)), ", ",
                                        sprintf("%.2f", round(meanCIUB, 2)), ")"),
                          median = paste0(sprintf("%.2f", round(medianSI, 2)), " (",
                                          sprintf("%.2f", round(medianCILB, 2)), ", ",
                                          sprintf("%.2f", round(medianCIUB, 2)), ")"),
                          sd = paste0(sprintf("%.2f", round(sdSI, 2)), " (",
                                      sprintf("%.2f", round(sdCILB, 2)), ", ",
                                      sprintf("%.2f", round(sdCIUB, 2)), ")"),
                          cutoff = ifelse(cutoff != "pooled", sprintf("%.3f", as.numeric(cutoff)),
                                          cutoff))
               %>% select(label, cutoff, npIncluded, mean, median, sd)
  )
  return(siTable2)
}

#Tables of pooled results for text
pooled <- formatSITable(siAll %>% filter(cutoff == "pooled"))
pooled %>% select(-npIncluded)

#### Supplementary Tables: Detailed Serial Interval Results ####
siHC <- formatSITable(siAll %>% filter(clustMethod == "hc_absolute",
                                       !grepl("Recent", label)))
siKD <- formatSITable(siAll %>% filter(clustMethod == "kd",
                                       !grepl("Recent", label)))


#Creating alternative label
siAll <- siAll %>% mutate(label2 = gsub("[A-Z]{2}: ", "", label),
                          label2 = factor(label2, levels = c("Excluding 3-month co-prevalent cases",
                                                             "Excluding 1-month co-prevalent cases",
                                                             "No exclusions",
                                                             "Recent Arrival = 1 Year")))

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
siAllLongPlot <- siAllLong2 %>% filter(cutoff != "pooled", !grepl("Recent", label2))
siAllLongPooled <- siAllLong2 %>% filter(cutoff == "pooled", !grepl("Recent", label2))
ggplot(data = siAllLongPlot, aes(x = as.numeric(cutoff), y = est, color = label2)) +
  geom_point() +
  geom_errorbar(aes(ymin = cilb, ymax = ciub, width = width)) +
  geom_hline(data = siAllLongPooled, aes(yintercept = est, color = label2)) +
  geom_hline(data = siAllLongPooled, aes(yintercept = cilb, color = label2), linetype = "dotted") +
  geom_hline(data = siAllLongPooled, aes(yintercept = ciub, color = label2), linetype = "dotted") +
  facet_grid(Parameter~clustMethod, scales = "free") +
  scale_x_continuous(name = "Clustering Cutoff/Binwidth") +
  scale_y_continuous(name = "Estimate in Years") +
  theme_bw() +
  #scale_color_grey(start = 0.7, end = 0.3) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                    size = 11),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                    size = 11)) +
  ggsave(file = "../Figures/MassSI.png",
         width = 7, height = 7, units = "in", dpi = 300)


## PRESENTATION VERSION ##
ggplot(data = siAllLongPlot, aes(x = as.numeric(cutoff), y = est, color = label2)) +
  geom_point() +
  geom_errorbar(aes(ymin = cilb, ymax = ciub, width = width)) +
  geom_hline(data = siAllLongPooled, aes(yintercept = est, color = label2)) +
  geom_hline(data = siAllLongPooled, aes(yintercept = cilb, color = label2), linetype = "dotted") +
  geom_hline(data = siAllLongPooled, aes(yintercept = ciub, color = label2), linetype = "dotted") +
  facet_grid(Parameter~clustMethod, scales = "free") +
  scale_x_continuous(name = "Clustering Cutoff/Binwidth") +
  scale_y_continuous(name = "Estimate in Years") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  guides(color = guide_legend(nrow = 2, byrow = FALSE)) +
  ggsave(file = "../Figures/MassSI_pres.png",
         width = 7, height = 7.5, units = "in", dpi = 300)



################### Reproductive Number ######################

#Average Rt accounting for importation
RtAvg

#Cutting the outbreak
totalTime <- max(RtData$timeRank) - min(RtData$timeRank)
monthCut1 <- ceiling(0.1 * totalTime)
monthCut2 <- ceiling(0.8 * totalTime)


#### Figure: Plot of Rt Estimates by Month with CIs ####
RtData2 <- RtData %>% filter(label == "Recent Arrival = 2 Years")
RtAvg2 <- RtAvg %>% filter(label == "Recent Arrival = 2 Years")
ggplot(data = RtData2, aes(x = timeRank, y = Rt)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = ciLower, ymax = ciUpper), width = 0.7, color = "grey40") +
  scale_y_continuous(name = "Monthly Effective Reproductive Number") + 
  scale_x_continuous(name = "Year of Observation", breaks = seq(3, 89, 12),
                     labels = seq(2010, 2017, 1)) +
  geom_vline(aes(xintercept = monthCut1), linetype = "dotted", size = 0.7) +
  geom_vline(aes(xintercept = monthCut2), linetype = "dotted", size = 0.7) +
  geom_hline(data = RtAvg2, aes(yintercept = RtAvg), size = 0.7) +
  theme_bw() +
  geom_hline(data = RtAvg2, aes(yintercept = ciLower), linetype = "dashed",
             size = 0.5, color = "grey40") +
  geom_hline(data = RtAvg2, aes(yintercept = ciUpper), linetype = "dashed",
             size = 0.5, color = "grey40") +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                    size = 12),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                    size = 12)) +
  ggsave(file = "../Figures/MassRt.png",
         width = 8, height = 6, units = "in", dpi = 300)


## PRESENTATION VERSION ##
ggplot(data = RtData2, aes(x = timeRank, y = Rt)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = ciLower, ymax = ciUpper), width = 0.7, color = "grey40") +
  scale_y_continuous(name = "Monthly Effective Reproductive Number") + 
  scale_x_continuous(name = "Year of Observation", breaks = seq(3, 89, 12),
                     labels = seq(2010, 2017, 1)) +
  geom_vline(aes(xintercept = monthCut1), linetype = "dotted", size = 0.7) +
  geom_vline(aes(xintercept = monthCut2), linetype = "dotted", size = 0.7) +
  geom_hline(data = RtAvg2, aes(yintercept = RtAvg), size = 0.7) +
  theme_bw(base_size = 16) +
  geom_hline(data = RtAvg2, aes(yintercept = ciLower), linetype = "dashed",
             size = 0.5, color = "grey40") +
  geom_hline(data = RtAvg2, aes(yintercept = ciUpper), linetype = "dashed",
             size = 0.5, color = "grey40") +
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ggsave(file = "../Figures/MassRt_pres.png",
         width = 8, height = 6, units = "in", dpi = 300)


  
################### Recent Arrival Sensitivity Analysis ######################


#### Supplementary Figure: Serial Interval ####

siSensLong <- (siAllLong2
               %>% filter(cutoff != "pooled",
                          label2 %in% c("No exclusions", "Recent Arrival = 1 Year"))
               %>% mutate(label2 = ifelse(label2 == "No exclusions", "Recent Arrival = 2 Years",
                                          as.character(label2)))
)
siSensLongPooled <- (siAllLong2
                     %>% filter(cutoff == "pooled",
                                label2 %in% c("No exclusions", "Recent Arrival = 1 Year"))
                     %>% mutate(label2 = ifelse(label2 == "No exclusions", "Recent Arrival = 2 Years",
                                                as.character(label2)))
)

ggplot(data = siSensLong, aes(x = as.numeric(cutoff), y = est, color = label2)) +
  geom_point() +
  geom_errorbar(aes(ymin = cilb, ymax = ciub, width = width)) +
  geom_hline(data = siSensLongPooled, aes(yintercept = est, color = label2)) +
  geom_hline(data = siSensLongPooled, aes(yintercept = cilb, color = label2), linetype = "dotted") +
  geom_hline(data = siSensLongPooled, aes(yintercept = ciub, color = label2), linetype = "dotted") +
  facet_grid(Parameter~clustMethod, scales = "free") +
  scale_x_continuous(name = "Clustering Cutoff/Binwidth") +
  scale_y_continuous(name = "Estimate in Years") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                    size = 11),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                    size = 11)) +
  ggsave(file = "../Figures/MassSISens.png",
         width = 7, height = 7, units = "in", dpi = 300)



#### Supplementary Figure: Reproductive Number ####

ggplot(data = RtData, aes(x = timeRank, y = Rt, color = label)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(name = "Monthly Effective Reproductive Number") + 
  scale_x_continuous(name = "Year of Observation", breaks = seq(3, 89, 12),
                     labels = seq(2010, 2017, 1)) +
  geom_vline(aes(xintercept = monthCut1), linetype = "dotted", size = 0.7) +
  geom_vline(aes(xintercept = monthCut2), linetype = "dotted", size = 0.7) +
  geom_hline(data = RtAvg, aes(yintercept = RtAvg, color = label), size = 0.7) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                    size = 12),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),
                                    size = 12)) +
ggsave(file = "../Figures/MassRtSens.png",
       width = 8, height = 6, units = "in", dpi = 300)



