#Sarah V. Leavitt
#Boston University Dissertation
#Paper 2

################################################################################
# This program makes tables and figures for the simulation results of the 
# generation interval estimation paper
################################################################################

#rm(list = ls())
options(scipen=999)
setwd("~/Boston University/Dissertation/Simulation_ResultsSI_1.21.20")

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)


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



########################### Formatting the Results ###########################

si2 <- (si
        #Reorder label
        %>% mutate(label = factor(label, levels = c("Baseline", "LowN", "HighN",
                                                    "LowR", "HighR", "LowMR",
                                                    "HighMR", "LowGV", "HighGV",
                                                    "LowGM", "HighGM"),
                                  labels = c("Baseline", "LowN", "HighN",
                                             "LowR", "HighR", "LowMR",
                                             "HighMR", "LowGV", "HighGV",
                                             "LowGM", "HighGM")),
                   obsCV = obsMean / obsSD,
                   cvSI = meanSI / sdSI)
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
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ggsave(file = "../Figures/GIResults.eps",
         width = 7, height = 7, units = "in", dpi = 600)

## COLOR VERSION ##
ggplot(data = plotData1, aes(y = absDiff, x = probf,
                             fill = Parameter, color = Parameter)) +
  facet_wrap(~label) +
  geom_violin(alpha = 0.7, draw_quantiles = 0.5) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(name = "Absolute Bias in Days", limits = c(-15, 15)) +
  scale_x_discrete(name = "Generation Interval Estimation Method") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ggsave(file = "../Figures/GIResults_color.png",
         width = 7, height = 8, units = "in", dpi = 300)

## PRESENTATION VERSION ##
ggplot(data = plotData1, aes(y = absDiff, x = probf,
                             fill = Parameter, color = Parameter)) +
  facet_wrap(~label) +
  geom_violin(alpha = 0.7, draw_quantiles = 0.5) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(name = "Absolute Bias in Days", limits = c(-15, 15)) +
  scale_x_discrete(name = "Generation Interval Estimation Method") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ggsave(file = "../Figures/GIResults_pres.png",
         width = 7, height = 8, units = "in", dpi = 300)


## PRESENTATION VERSION - ABBREVIATED ##
ggplot(data = plotData1 %>% filter(label == "Baseline"), aes(y = absDiff, x = probf,
                             fill = Parameter, color = Parameter)) +
  geom_violin(alpha = 0.7, draw_quantiles = 0.5) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(name = "Absolute Bias in Days", limits = c(-15, 15)) +
  scale_x_discrete(name = "Generation Interval Estimation Method") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ggsave(file = "../Figures/GIResults_abb.png",
         width = 7, height = 7, units = "in", dpi = 300)




####################### Supplementary Figures ##########################


#### Supplementary Figure: MAPE Plot ####

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
        axis.text.x = element_text(angle = 45, hjust = 1),
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
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

grid.arrange(pErrorHC, pErrorKD, ncol = 2)
pError <- arrangeGrob(pErrorHC, pErrorKD, ncol = 2)
ggsave(file = "../Figures/GIError.png", plot = pError,
       width = 9, height = 7, units = "in", dpi = 300)


## COLOR VERSION ##
pErrorHCc <- ggplot(data = errorL %>% filter(grepl("^HC", prob),
                                             prob != "HCpooled",
                                             prob != "HC0"),
                    aes(x = as.numeric(cutoff), y = error, color = Parameter)) +
  facet_wrap(~label) +
  geom_point() +
  scale_y_continuous(name = "Mean Absolute Percentage Error", limits = c(0, 45)) +
  scale_x_continuous(name = "Hierarchical Clustering Cutoff") +
  geom_hline(data = errorL %>% filter(prob == "HCpooled"),
             aes(yintercept = error, color = Parameter)) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))


pErrorKDc <- ggplot(data = errorL %>% filter(grepl("^KD", prob),
                                             prob != "KDpooled"),
                    aes(x = as.numeric(cutoff), y = error, color = Parameter)) +
  facet_wrap(~label) +
  geom_point() +
  scale_y_continuous(name = "Mean Absolute Percentage Error", limits = c(0, 45)) +
  scale_x_continuous(name = "Kernel Density Estimation Binwidth") +
  geom_hline(data = errorL %>% filter(prob == "KDpooled"),
             aes(yintercept = error, color = Parameter)) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

grid.arrange(pErrorHCc, pErrorKDc, ncol = 2)
pErrorc <- arrangeGrob(pErrorHCc, pErrorKDc, ncol = 2)
ggsave(file = "../Figures/GIError_color.png", plot = pErrorc,
       width = 9, height = 7, units = "in", dpi = 300)


## PRESENTATION VERSION ##
pErrorHCp <- ggplot(data = errorL %>% filter(grepl("^HC", prob),
                                             prob != "HCpooled",
                                             prob != "HC0"),
                    aes(x = as.numeric(cutoff), y = error, color = Parameter)) +
  facet_wrap(~label) +
  geom_point() +
  scale_y_continuous(name = "Mean Absolute Percentage Error", limits = c(0, 45)) +
  scale_x_continuous(name = "Hierarchical Clustering Cutoff") +
  geom_hline(data = errorL %>% filter(prob == "HCpooled"),
             aes(yintercept = error, color = Parameter)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))


pErrorKDp <- ggplot(data = errorL %>% filter(grepl("^KD", prob),
                                             prob != "KDpooled"),
                    aes(x = as.numeric(cutoff), y = error, color = Parameter)) +
  facet_wrap(~label) +
  geom_point() +
  scale_y_continuous(name = "Mean Absolute Percentage Error", limits = c(0, 45)) +
  scale_x_continuous(name = "Kernel Density Estimation Binwidth") +
  geom_hline(data = errorL %>% filter(prob == "KDpooled"),
             aes(yintercept = error, color = Parameter)) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

grid.arrange(pErrorHCp, pErrorKDp, ncol = 2)
pErrorp <- arrangeGrob(pErrorHCp, pErrorKDp, ncol = 2)
ggsave(file = "../Figures/GIError_pres.png", plot = pErrorp,
       width = 9, height = 7, units = "in", dpi = 300)



#### Supplementary Figure: Performance Metrics ####
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





###################### Extra Results ############################

setwd("~/Boston University/Dissertation/Simulation_ResultsSI_6.30.20")


#Initializing dataframes
si_cov <- NULL

#Reading in the results
for (file in list.files()){
  
  if(grepl("^si", file)){
    siTemp <- readRDS(file)
    si_cov <- bind_rows(si_cov, siTemp)
  }
}

si_cov <- si_cov %>% mutate(meanInclude = obsMean > meanCILB & obsMean < meanCIUB,
                            medianInclude = obsMedian > medianCILB & obsMedian < medianCIUB,
                            sdInclude = obsSD > sdCILB & obsSD < sdCIUB)


#### Extra Figure: Time as a Covariate ####

plotTime <- (si_cov
             %>% group_by(prob)
             %>% filter(clustMethod %in% c("kd", "kdt"))
             %>% mutate(clustMethodf = factor(clustMethod, levels = c("kd", "kdt"),
                                              labels = c("Excluding time",
                                                         "Including time")))
             %>% summarize(nRuns = sum(!is.na(meanDiff)),
                           avgNumInd = mean(nIndividuals, na.rm = TRUE),
                           Mean = 100 * mean(abs(meanDiff / obsMean), na.rm = TRUE),
                           Median = 100 * mean(abs(medianDiff / obsMedian), na.rm = TRUE),
                           SD = 100 * mean(abs(sdDiff / obsSD), na.rm = TRUE),
                           clustMethodf = first(clustMethodf),
                           cutoff = first(cutoff),
                           .groups = "drop")
             %>% select(clustMethodf, cutoff, prob, Mean, Median, SD)
             %>% gather("Parameter", "error", -prob, -clustMethodf, -cutoff)
)

ggplot(data = plotTime %>% filter(cutoff != "pooled"),
       aes(x = as.numeric(cutoff), y = error, color = clustMethodf)) +
  facet_wrap(~Parameter) +
  geom_point() +
  scale_y_continuous(name = "Mean Absolute Percentage Error", limits = c(0, 25)) +
  scale_x_continuous(name = "Kernel Density Estimation Binwidth") +
  scale_color_grey(start = 0.3, end = 0.7) +
  geom_hline(data = plotTime %>% filter(cutoff == "pooled"),
             aes(yintercept = error, color = clustMethodf)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))



#### Supplementary Figure: Coverage plot ####

coverage <- (si_cov
             %>% filter(clustMethod != "kdt")
             %>% group_by(prob)
             %>% summarize(nRuns = sum(!is.na(meanDiff)),
                           Mean = 100 * sum(meanInclude) / nRuns,
                           Median = 100 * sum(medianInclude) / nRuns,
                           SD = 100 * sum(sdInclude) / nRuns,
                           clustMethod = first(clustMethod),
                           cutoff = first(cutoff),
                           .groups = "drop")
             %>% select(clustMethod, cutoff, prob, Mean, Median, SD)
             %>% gather("Parameter", "coverage", -prob, -clustMethod, -cutoff)
             %>% filter(!is.na(coverage))
             %>% mutate(clustMethodf = factor(clustMethod, levels = c("hc_absolute", "kd"),
                                              labels = c("Hierarchical Cluster",
                                                         "Kernel Density Estimation")))
)

ggplot(data = coverage %>% filter(cutoff != "pooled"),
       aes(x = as.numeric(cutoff), y = coverage, color = Parameter)) +
  geom_point() +
  facet_wrap(~clustMethodf, scales = "free_x") +
  scale_y_continuous(name = "Coverage", limits = c(70, 100)) +
  scale_x_continuous(name = "Clustering Cutoff or Binwdith") +
  geom_hline(data = coverage %>% filter(cutoff == "pooled"),
             aes(yintercept = coverage, color = Parameter)) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ggsave(file = "../Figures/GISuppCov.png",
         width = 7, height = 5, units = "in", dpi = 300)


