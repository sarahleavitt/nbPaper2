#Sarah V. Leavitt
#Boston University Dissertation
#Paper 2

################################################################################
# This program makes the clustering examples figure for paper 2
################################################################################

#rm(list = ls())
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




######################### Simulate Example Outbreak ##########################

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
pTraining = 0.6
trainingID <- (covarInd
               %>% filter(complete == TRUE, !is.na(sampleDate))
               %>% sample_frac(pTraining)
               %>% pull(individualID)
)

orderedPair <- (covarPair
                %>% filter(infectionDate.2 > infectionDate.1)
                %>% mutate(snpClose = ifelse(snpDist < lowerT, TRUE,
                                             ifelse(snpDist > upperT, FALSE, NA)),
                           trainPair = individualID.1 %in% trainingID & individualID.2 %in% trainingID,
                           snpCloseGS = ifelse(trainPair == TRUE, snpClose, NA),
                           snpClose2 = snpDist < 3)
)


resGen <- nbProbabilities(orderedPair = orderedPair, indIDVar = "individualID",
                          pairIDVar = "edgeID", goldStdVar = "snpCloseGS",
                          covariates = c("Y1", "Y2", "Y3", "Y4"), label = "NoTime",
                          n = 10, m = 1, nReps = 10)
nbResults <- resGen[[1]] %>% full_join(orderedPair, by = "edgeID")
print("Completed SNP threshold gold standard analysis")





############################# Clustering infectors ################################

clustRes <- clusterInfectors(nbResults, indIDVar = "individualID", pVar = "pScaled",
                             clustMethod = "hc_absolute", cutoff = 0.05)

#Finding good example cases
ggplot(data = clustRes %>% filter(individualID.2 >= 10010, individualID.2 <= 10050),
       aes(x = pRank, y = pScaled, color = cluster, shape = transmission)) +
  geom_point() +
  facet_wrap(~individualID.2, scales = "free") +
  theme(legend.position = "none")

examples <- clustRes %>% filter(individualID.2 %in% c(10047, 10025))
ind1 <- examples %>% filter(individualID.2 == 10047) %>% arrange(pRank)
ind2 <- examples %>% filter(individualID.2 == 10025) %>% arrange(pRank)

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
  xlab("Probability Rank") +
  ylab("Height")

pHC2 <- ggplot(as.ggdend(hclustD2)) +
  theme_bw() +
  scale_color_manual(values = c("darkgrey", "black", "red")) +
  xlab("Probability Rank") +
  ylab("Height")




#### Plot of densities ####

pKD1 <- findClustersKD(df = ind1, pVar = "pScaled", cutoff = 0.02, plot = TRUE,
                       colors = c("black", "darkgrey"), size = 12)
pKD2 <- findClustersKD(df = ind2, pVar = "pScaled", cutoff = 0.02, plot = TRUE,
                       colors = c("black", "darkgrey"), size = 12)


#### Figure: Clustering Examples ####

lay <- rbind(c(1, 2), c(3, 4), c(5, 6))
grid.arrange(p1, p2, pHC1, pHC2, pKD1, pKD2, layout_matrix = lay)

pAll <- arrangeGrob(p1, p2, pHC1, pHC2, pKD1, pKD2, layout_matrix = lay)
ggsave(file = "../Figures/ClustExamples.pdf", plot = pAll,
       width = 6, height = 7, units = "in", dpi = 600)




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

pHC1c <- ggplot(as.ggdend(hclustD1)) +
  theme_bw() +
  scale_color_manual(values = c("#F8766D", "#00BFC4", "red")) +
  xlab("Probability Rank") +
  ylab("Height")

pHC2c <- ggplot(as.ggdend(hclustD2)) +
  theme_bw() +
  scale_color_manual(values = c("#F8766D", "#00BFC4", "red")) +
  xlab("Probability Rank") +
  ylab("Height")

pKD1c <- findClustersKD(df = ind1, pVar = "pScaled", cutoff = 0.02, plot = TRUE,
                        colors = c("#00BFC4", "#F8766D"), size = 12)
pKD2c <- findClustersKD(df = ind2, pVar = "pScaled", cutoff = 0.02, plot = TRUE,
                        colors = c("#00BFC4", "#F8766D"), size = 12)

grid.arrange(p1c, p2c, pHC1c, pHC2c, pKD1c, pKD2c, layout_matrix = lay)
pAllc <- arrangeGrob(p1c, p2c, pHC1c, pHC2c, pKD1c, pKD2c, layout_matrix = lay)
ggsave(file = "../Figures/ClustExamples_color.png", plot = pAllc,
       width = 6, height = 7, units = "in", dpi = 300)



## PRESENTATION VERSIONS ##
p1p <- ggplot(data = ind1, aes(x = pRank, y = pScaled, color = cluster)) +
  geom_jitter() +
  xlab("Probability Rank") +
  ylab("Relative Probability") +
  scale_color_manual(values = c("#00BFC4", "#F8766D"), drop = FALSE) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none") +
  ggtitle("Case A")

p2p <- ggplot(data = ind2, aes(x = pRank, y = pScaled, color = cluster)) +
  geom_jitter() +
  xlab("Probability Rank") +
  ylab("Relative Probability") +
  scale_color_manual(values = c("#00BFC4", "#F8766D"), drop = FALSE) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none") +
  ggtitle("Case B")

pHC1p <- ggplot(as.ggdend(hclustD1)) +
  theme_bw(base_size = 16) +
  scale_color_manual(values = c("#F8766D", "#00BFC4", "red")) +
  xlab("Probability Rank") +
  ylab("Height")

pHC2p <- ggplot(as.ggdend(hclustD2)) +
  theme_bw(base_size = 16) +
  scale_color_manual(values = c("#F8766D", "#00BFC4", "red")) +
  xlab("Probability Rank") +
  ylab("Height")

pKD1p <- findClustersKD(df = ind1, pVar = "pScaled", cutoff = 0.02, plot = TRUE,
                        colors = c("#00BFC4", "#F8766D"), size = 16)
pKD2p <- findClustersKD(df = ind2, pVar = "pScaled", cutoff = 0.02, plot = TRUE,
                        colors = c("#00BFC4", "#F8766D"), size = 16)

grid.arrange(p1p, p2p, pHC1p, pHC2p, pKD1p, pKD2p, layout_matrix = lay)
pAllp <- arrangeGrob(p1p, p2p, pHC1p, pHC2p, pKD1p, pKD2p, layout_matrix = lay)
ggsave(file = "../Figures/ClustExamples_pres.png", plot = pAllp,
       width = 6, height = 7, units = "in", dpi = 300)
