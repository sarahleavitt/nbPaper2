#Sarah V. Leavitt
#Boston University Dissertation
#Paper 2

################################################################################
# This program makes tables and figures for the simulation results of the 
# generation interval estimation paper
################################################################################

rm(list = ls())
options(scipen=999)
setwd("~/Boston University/Dissertation/nbPaper2")

library(TransPhylo)
library(phangorn)
library(reshape2)
library(dplyr)
library(tidyr)
library(lubridate)
library(devtools)
library(ggplot2)
library(gridExtra)
library(dendextend)
library(purrr)

load_all("../nbTransmission")
source("../nbSimulation/SimOutbreak.R")
source("../nbSimulation/SimulateOutbreakS.R")
source("../nbSimulation/SimCovariates.R")
source("../nbSimulation/SimEvaluate.R")




##################### Clustering Example Plot ######################


#### Simulate Example Outbreak ####

#Parameters to change
sampleSize <- 50
off.r <- 1.5
w.shape <- 2.25
w.scale <- 0.0122
w.shift <- 0
mutationR <- 25 #snps/genome/year
lowerT <- 3
upperT <- 7
iPar1 <- 3
iPar2 <- 4
#Constant parameters
neg <- 0.25
pi <- 1
off.p <- 0.5
multOutbreaks <- FALSE
rootseq <- NULL
length <- 3000
#Derived parameters
rate <- mutationR / length
time <- (log(sampleSize, off.r) + 1)*qgamma(0.5, shape = w.shape, scale = w.scale)


#Simulate outbreak
set.seed(10001)
obk <- simOutbreak(neg = neg, pi = pi, off.r = off.r, off.p = off.p,
                   w.scale = w.scale, w.shape = w.shape, w.shift = w.shift,
                   ws.scale = w.scale, ws.shape = w.shape, ws.shift = w.shift,
                   sampleSize = sampleSize, time = time, startDate = 2010,
                   multOutbreaks = multOutbreaks, length = length, rate = rate)

indData <- obk[[1]]
pairData <- obk[[2]]
print(paste0("Simulated outbreak, n = ", nrow(indData)))

#Simulating covariates
covar <- simCovariates(indData, pairData)
covarPair <- covar[[1]]
covarInd <- covar[[2]]
print("Simulated covariates")

#Only using a proportion for training
pTraining = 1
trainingID <- (covarInd
               %>% filter(complete == TRUE, !is.na(sampleDate))
               %>% sample_frac(pTraining)
               %>% pull(individualID)
)

covarOrderedPair <- (covarPair
                     %>% filter(infectionDate.2 > infectionDate.1)
                     %>% mutate(snpClose = ifelse(snpDist < lowerT, TRUE,
                                           ifelse(snpDist > upperT, FALSE, NA)),
                                trainPair = individualID.1 %in% trainingID & individualID.2 %in% trainingID,
                                snpCloseGS = ifelse(trainPair == TRUE, snpClose, NA),
                                snpClose2 = snpDist < 3)
)


resGen <- nbProbabilities(orderedPair = covarOrderedPair, indIDVar = "individualID",
                          pairIDVar = "edgeID", goldStdVar = "snpCloseGS",
                          covariates = c("Y1", "Y2", "Y3", "Y4"), label = "NoTime",
                          n = 10, m = 1, nReps = 10)
nbResults <- resGen[[1]] %>% full_join(covarOrderedPair, by = "edgeID")
print("Completed SNP threshold gold standard analysis")

saveRDS(nbResults, "../Datasets/ClustExample.rds")



#### Clustering infectors ####

#nbResults <- readRDS("../Datasets/ClustExample.rds")

clustRes <- clusterInfectors(nbResults, indIDVar = "individualID", pVar = "pScaled",
                          clustMethod = "hc_absolute", cutoff = 0.05)

#Finding good example cases
ggplot(data = clustRes %>% filter(individualID.2 >= 10070, individualID.2 <= 10099),
       aes(x = pRank, y = pScaled, color = cluster, shape = transmission)) +
  geom_point() +
  facet_wrap(~individualID.2, scales = "free") +
  theme(legend.position = "none")

examples <- clustRes %>% filter(individualID.2 %in% c(10090, 10089))
ind1 <- examples %>% filter(individualID.2 == 10090) %>% arrange(pRank)
ind2 <- examples %>% filter(individualID.2 == 10089) %>% arrange(pRank)

ggplot(data = examples, aes(x = pRank, y = pScaled, color = cluster, shape = transmission)) +
  geom_jitter() +
  facet_wrap(~individualID.2) +
  xlab("Probability Rank") +
  ylab("Relative Probability") +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  theme_bw() +
  theme(legend.position = "none")


#### Plot of probabilities ####

#Note this gives the default blue and red c("#00BFC4", "#F8766D")
p1 <- ggplot(data = ind1, aes(x = pRank, y = pScaled, color = cluster)) +
  geom_jitter() +
  xlab("Probability Rank") +
  ylab("Relative Probability") +
  scale_color_manual(values = c("black", "darkgrey"), drop = FALSE) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Case A")

p2 <- ggplot(data = ind2, aes(x = pRank, y = pScaled, color = cluster)) +
  geom_jitter() +
  xlab("Probability Rank") +
  ylab("Relative Probability") +
  scale_color_manual(values = c("black", "darkgrey"), drop = FALSE) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Case B")

## COLOR VERSIONS ##
p1c <- ggplot(data = ind1, aes(x = pRank, y = pScaled, color = cluster)) +
  geom_jitter() +
  xlab("Probability Rank") +
  ylab("Relative Probability") +
  scale_color_manual(values = c("#00BFC4", "#F8766D"), drop = FALSE) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Case A")

p2c <- ggplot(data = ind2, aes(x = pRank, y = pScaled, color = cluster)) +
  geom_jitter() +
  xlab("Probability Rank") +
  ylab("Relative Probability") +
  scale_color_manual(values = c("#00BFC4", "#F8766D"), drop = FALSE) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Case B")


#### Plot of dendrograms ####

hclustS1 <- rev(hclust(dist(ind1$pScaled), method = "single"))
hclustD1 <- (hclustS1
            %>% as.dendrogram(.)
            %>% set("branches_k_color", k = 2)
            %>% set("hang_leaves", 0.02)
            %>% set("branches_lwd", 0.5)
            %>% set("labels", NA)
)
hclustS2 <- rev(hclust(dist(ind2$pScaled), method = "single"))
hclustD2 <- (hclustS2
             %>% as.dendrogram(.)
             %>% set("branches_k_color", k = 1)
             %>% set("hang_leaves", 0.02)
             %>% set("branches_lwd", 0.5)
             %>% set("labels", NA)
)

pHC1 <- ggplot(as.ggdend(hclustD1)) +
  theme_bw() +
  scale_color_manual(values = c("darkgrey", "black", "red")) +
  scale_x_continuous(limits = c(0, 40)) +
  xlab("Probability Rank") +
  ylab("Height")

pHC2 <- ggplot(as.ggdend(hclustD2)) +
  theme_bw() +
  scale_color_manual(values = c("darkgrey", "black", "red")) +
  scale_x_continuous(limits = c(0, 90)) +
  xlab("Probability Rank") +
  ylab("Height")

## COLOR VERSIONS ##
pHC1c <- ggplot(as.ggdend(hclustD1)) +
  theme_bw() +
  scale_color_manual(values = c("#F8766D", "#00BFC4", "red")) +
  scale_x_continuous(limits = c(0, 40)) +
  xlab("Probability Rank") +
  ylab("Height")

pHC2c <- ggplot(as.ggdend(hclustD2)) +
  theme_bw() +
  scale_color_manual(values = c("#F8766D", "#00BFC4", "red")) +
  scale_x_continuous(limits = c(0, 90)) +
  xlab("Probability Rank") +
  ylab("Height")



#### Plot of densities ####

pKD1 <- findClustersKD(df = ind1, pVar = "pScaled", cutoff = 0.01, plot = TRUE,
                       colors = c("black", "darkgrey"))
pKD2 <- findClustersKD(df = ind2, pVar = "pScaled", cutoff = 0.01, plot = TRUE,
                       colors = c("black", "darkgrey"))

## COLOR VERSIONS ##
pKD1c <- findClustersKD(df = ind1, pVar = "pScaled", cutoff = 0.01, plot = TRUE,
                        colors = c("#00BFC4", "#F8766D"))
pKD2c <- findClustersKD(df = ind2, pVar = "pScaled", cutoff = 0.01, plot = TRUE,
                        colors = c("#00BFC4", "#F8766D"))




#### Figure: Clustering Examples ####

lay <- rbind(c(1, 2), c(3, 4), c(5, 6))
grid.arrange(p1, p2, pHC1, pHC2, pKD1, pKD2, layout_matrix = lay)

pAll <- arrangeGrob(p1, p2, pHC1, pHC2, pKD1, pKD2, layout_matrix = lay)
ggsave(file = "../Figures/ClustExamples.png", plot = pAll,
       width = 6, height = 7, units = "in", dpi = 300)

## COLOR VERSION ##
grid.arrange(p1c, p2c, pHC1c, pHC2c, pKD1c, pKD2c, layout_matrix = lay)
pAllc <- arrangeGrob(p1c, p2c, pHC1c, pHC2c, pKD1c, pKD2c, layout_matrix = lay)
ggsave(file = "../Figures/ClustExamples_color.png", plot = pAllc,
       width = 6, height = 7, units = "in", dpi = 300)




###################### Simulation Results #########################

setwd("~/Boston University/Dissertation/Simulation_ResultsSI_1.21.20")

#Initializing dataframes
si <- NULL
perform <- NULL

#Reading in the results
for (file in list.files()){
  
  if(grepl("^si", file)){
    siTemp <- readRDS(file)
    si <- bind_rows(si, siTemp)
  }
  
  if(grepl("^perform", file)){
    pTemp <- readRDS(file)
    perform <- bind_rows(perform, pTemp) 
  }
}

si2 <- (si
        #Reorder label
        #Renaming incorrect columns
        %>% mutate(label = factor(label, levels = c("Baseline", "LowN", "HighN",
                                                    "LowR", "HighR", "LowMR",
                                                    "HighMR", "LowGV", "HighGV",
                                                    "LowGM", "HighGM"),
                                  labels = c("Baseline", "LowN", "HighN",
                                             "LowR", "HighR", "LowMR",
                                             "HighMR", "LowGV", "HighGV",
                                             "LowGM", "HighGM")),
                   obsCV = obsMean / obsSD,
                   cvSI = meanSI / sdSI,
                   relMeanDiff = meanDiff/obsMean,
                   relMedianDiff = medianDiff/obsMean,
                   relSDDiff = sdDiff / obsSD)
        %>% filter(!label %in% c("LowGM", "HighGM"))
)


#### Figure: Violin Plot of Pooled GI Estimates ####

plotData1 <- (si2
              %>% filter(prob %in% c("All", "Naive", "Npooled", "KDpooled", "HCpooled"))
              %>% mutate(probf = factor(prob, levels = c("Naive", "All", "Npooled",
                                                         "HCpooled", "KDpooled"),
                                        labels = c("SNP Distance", "PEM: All Pairs",
                                                   "PEM: Top N", "PEM: Hierarchical",
                                                   "PEM: Kernel Density")))
              %>% select(label, probf, Mean = meanDiff, Median = medianDiff, SD = sdDiff)
              %>% gather("Parameter", "absDiff", -label, -probf)
)

ggplot(data = plotData1, aes(y = absDiff, x = probf,
                             fill = Parameter, color = Parameter)) +
  facet_wrap(~label) +
  geom_violin(alpha = 0.7, draw_quantiles = 0.5) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(name = "Absolute Bias in Days", limits = c(-15, 15)) +
  scale_x_discrete(name = "Generation Interval Estimation Method") +
  scale_fill_grey(start = 0.3, end = 0.7) +
  scale_color_grey(start = 0.3, end = 0.7) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ggsave(file = "../Figures/GIResults.png",
         width = 7, height = 7, units = "in", dpi = 300)


#### Figure: MAPE Plot ####
errorL <- (si2
           %>% group_by(label, prob)
           %>% summarize(nRuns = sum(!is.na(meanDiff)),
                         avgNumInd = mean(nIndividuals, na.rm = TRUE),
                         Mean = 100 * mean(abs(meanDiff / obsMean), na.rm = TRUE),
                         Median = 100 * mean(abs(medianDiff / obsMedian), na.rm = TRUE),
                         SD = 100 * mean(abs(sdDiff / obsSD), na.rm = TRUE),
                         clustMethod = first(clustMethod),
                         cutoff = first(cutoff))
           %>% select(label, clustMethod, cutoff, prob, Mean, Median, SD)
           %>% gather("Parameter", "error", -label, -prob, -clustMethod, -cutoff)
)

pErrorHC <- ggplot(data = errorL %>% filter(grepl("^HC", prob),
                                prob != "HCpooled",
                                prob != "HC0"),
       aes(x = as.numeric(cutoff), y = error, color = Parameter)) +
  facet_wrap(~label) +
  geom_point() +
  scale_y_continuous(name = "Mean Absolute Percentage Error", limits = c(0, 45)) +
  scale_x_continuous(name = "Hierarchical Clustering Cutoff") +
  scale_color_grey(start = 0.3, end = 0.7) +
  geom_hline(data = errorL %>% filter(prob == "HCpooled"),
             aes(yintercept = error, color = Parameter)) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))


pErrorKD <- ggplot(data = errorL %>% filter(grepl("^KD", prob),
                                prob != "KDpooled"),
       aes(x = as.numeric(cutoff), y = error, color = Parameter)) +
  facet_wrap(~label) +
  geom_point() +
  scale_y_continuous(name = "Mean Absolute Percentage Error", limits = c(0, 45)) +
  scale_x_continuous(name = "Kernel Density Estimation Binwidth") +
  scale_color_grey(start = 0.3, end = 0.7) +
  geom_hline(data = errorL %>% filter(prob == "KDpooled"),
             aes(yintercept = error, color = Parameter)) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

grid.arrange(pErrorHC, pErrorKD, ncol = 2)
pError <- arrangeGrob(pErrorHC, pErrorKD, ncol = 2)
ggsave(file = "../Figures/GIError.png", plot = pError,
       width = 9, height = 7, units = "in", dpi = 300)



#### Figure: Performance Metrics ####
longData <- (perform
             %>% ungroup()
             %>% mutate(label = factor(label, levels = c("Baseline", "LowN", "HighN",
                                                         "LowR", "HighR", "LowMR",
                                                         "HighMR", "LowGIV", "HighGIV",
                                                         "LowGIM", "HighGIM"),
                                       labels = c("Baseline", "LowN", "HighN",
                                                  "LowR", "HighR", "LowMR",
                                                  "HighMR", "LowGV", "HighGV",
                                                  "LowGM", "HighGM")))
             %>% filter(!label %in% c("LowGM", "HighGM"))
             %>% select(runID, label, aucVal, pCorrect,
                        pTop5, pTop10, pTop25, pTop50)
             %>% gather(metric, value, -label, -runID)
             %>% mutate(metric = factor(metric, levels = c("aucVal", "pCorrect", "pTop5",
                                                           "pTop10", "pTop25", "pTop50"),
                                        labels = c("Area under the ROC",
                                                   "Proportion Correct",
                                                   "Proportion in Top 5%",
                                                   "Proportion in Top 10%",
                                                   "Proportion in Top 25%",
                                                   "Proportion in Top 50%")))
)

ggplot(data = longData, aes(x = label, y = value)) +
  geom_violin(alpha = 0.7, draw_quantiles = 0.5, fill = "grey") +
  scale_x_discrete(name = "Scenario") +
  scale_y_continuous(name = "Value") +
  facet_wrap(~ metric, nrow = 3, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.position = "none") +
  ggsave(file = "../Figures/GIMetrics.png",
         width = 8, height = 5, units = "in", dpi = 300)






####################### Supplementary Figures ##########################

suppPlots <- (si2
              %>% select(label, prob, Mean = meanDiff, Median = medianDiff, SD = sdDiff)
              %>% gather("Parameter", "absDiff", -label, -prob)
)

#### Supplementary Figure: Hierarchical Clustering Density ####
ggplot(data = suppPlots %>% filter(grepl("^HC", prob)),
       aes(y = absDiff, x = gsub("HC", "", prob), fill = Parameter, color = Parameter)) +
  facet_wrap(~label) +
  geom_violin(alpha = 0.7, draw_quantiles = 0.5) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(name = "Absolute Bias in Days", limits = c(-10.3, 10.3)) +
  scale_x_discrete(name = "Hierarchical Clustering Cutoff") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ggsave(file = "../Figures/GISuppHC.png",
         width = 9, height = 7, units = "in", dpi = 300)


#### Supplementary Figure: Kernel Density ####
ggplot(data = suppPlots %>% filter(grepl("^KD", prob)),
       aes(y = absDiff, x = gsub("KD", "", prob), fill = Parameter, color = Parameter)) +
  facet_wrap(~label) +
  geom_violin(alpha = 0.7, draw_quantiles = 0.5) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(name = "Absolute Bias in Days", limits = c(-10.3, 10.3)) +
  scale_x_discrete(name = "Kernel Density Estimation Binwidth") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ggsave(file = "../Figures/GISuppKD.png",
         width = 9, height = 7, units = "in", dpi = 300)


#### Supplementary Figure: Top N ####
suppPlotsN <- (suppPlots
               %>% filter(grepl("^N", prob), prob != "Naive")
               %>% mutate(prob = factor(prob, levels = c(paste0("N", 1:10), "Npooled"),
                                        labels = c(1:10, "pooled")))
)
ggplot(data = suppPlotsN, aes(y = absDiff, x = prob,
                              fill = Parameter, color = Parameter)) +
  facet_wrap(~label) +
  geom_violin(alpha = 0.7, draw_quantiles = 0.5) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(name = "Absolute Bias in Days", limits = c(-10, 15)) +
  scale_x_discrete(name = "Number of Infectors Included") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ggsave(file = "../Figures/GISuppN.png",
         width = 9, height = 7, units = "in", dpi = 300)




#### Histograms of Number of Missing ####
sumScenario <- (si2
                %>% group_by(label, prob)
                %>% summarize(nMiss = sum(is.na(meanSI)))
                #nIndividuals = mean(nIndividuals, na.rm = TRUE),
                #nInfectors = mean(nInfectors), na.rm = TRUE))
)

ggplot(data = sumScenario %>% filter(grepl("^HC", prob)),
       aes(y = 100 * nMiss/1000, x = gsub("HC", "", prob))) +
  facet_wrap(~label) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0) +
  scale_y_continuous(name = "Percentage of Missing Estimates") +
  scale_x_discrete(name = "Cutoff Value") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

ggplot(data = sumScenario %>% filter(grepl("^KD", prob)),
       aes(y = 100 * nMiss/1000, x = gsub("KD", "", prob))) +
  facet_wrap(~label) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0) +
  scale_y_continuous(name = "Percentage of Missing Estimates") +
  scale_x_discrete(name = "Cutoff Value") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
