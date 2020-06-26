#Sarah V. Leavitt
#Boston University Dissertation
#Paper 2

################################################################################
# This program runs one iteration of a simulation to assess
# generation interval estimation
################################################################################


SimRunSI <- function(label, observationDate) {
  
  
  #### Simulate outbreak ####
  
  shift <- w.shift
  obk <- simOutbreak(neg = neg, pi = pi, off.r = off.r, off.p = off.p,
                     w.scale = w.scale, w.shape = w.shape, w.shift = w.shift,
                     ws.scale = ws.scale, ws.shape = ws.shape, ws.shift = ws.shift,
                     sampleSize = sampleSize, time = time, startDate = 2010,
                     multOutbreaks = multOutbreaks, length = length, rate = rate)
  
  indData <- obk[[1]]
  pairData <- obk[[2]]
  print(paste0("Simulated outbreak, n = ", nrow(indData)))
  
  #Simulating covariates
  covar <- simCovariates(indData, pairData, observationDate = observationDate)
  covarPair <- covar[[1]]
  covarInd <- covar[[2]]
  print("Simulated covariates")

  #Subseting to the pairs with the potential infector observed before the infectee
  covarOrderedPair <- (covarPair
                       %>% filter(infectionDate.2 > infectionDate.1)
                       %>% mutate(snpClose = ifelse(snpDist < lowerT, TRUE,
                                                    ifelse(snpDist > upperT, FALSE, NA)))
  )
  
  #Finding the biase for mean, median, and variance
  truePairs <- covarPair %>% filter(transmission == TRUE)
  obsMean <- mean(truePairs[, dateVar])
  obsMedian <- median(truePairs[, dateVar])
  obsSD <- sd(truePairs[, dateVar])
  
  
  #### Choosing pTraining cases for training set ####
  
  #Make sure that the training dataset has at least 3 true links (relevant for small sample sizes)
  nLinked <- 0
  nTries <- 0
  while(nLinked < 3){
    #Finding all pairs that can be included in the training dataset (no missing variables)
    trainingID <- (covarInd
                   %>% filter(complete == TRUE, !is.na(sampleDate))
                   %>% sample_frac(pTraining)
                   %>% pull(individualID)
    )
    
    covarOrderedPair <- (covarOrderedPair
                         %>% mutate(trainPair = ifelse(individualID.1 %in% trainingID &
                                                         individualID.2 %in% trainingID,
                                                       TRUE, FALSE))
    )
    nLinked <- sum(covarOrderedPair$trainPair == TRUE &
                     covarOrderedPair$snpClose == TRUE, na.rm = TRUE)
    nTries <- nTries + 1
  }
  if(nTries > 1){print(paste0("Number of tries: ", nTries))}
  
  
  
  #### Estimating probabilities ####
  
  covarOrderedPair <- covarOrderedPair %>% mutate(snpCloseGS = ifelse(trainPair == TRUE, snpClose, NA))
  
  covariates <- c("Y1", "Y2", "Y3", "Y4")
  resGen <- nbProbabilities(orderedPair = covarOrderedPair, indIDVar = "individualID",
                              pairIDVar = "edgeID", goldStdVar = "snpCloseGS",
                              covariates = covariates, label = label,
                              n = 10, m = 1, nReps = 10)
  allProbs <- resGen[[1]] %>% full_join(covarOrderedPair, by = "edgeID")
  
  print("Completed SNP threshold gold standard analysis")
  
  
  covariates2 <- c("Y1", "Y2", "Y3", "Y4", "timeCat")
  resGen2 <- nbProbabilities(orderedPair = covarOrderedPair, indIDVar = "individualID",
                            pairIDVar = "edgeID", goldStdVar = "snpCloseGS",
                            covariates = covariates2, label = label,
                            n = 10, m = 1, nReps = 10)
  allProbs2 <- resGen2[[1]] %>% full_join(covarOrderedPair, by = "edgeID")
  
  print("Completed SNP threshold gold standard analysis with time")
  
  
  
  #### Evaluating the performance ####

  pTemp <- (allProbs
            %>% group_by(label)
            %>% do(simEvaluate(.))
            %>% mutate(nTrainLinks = sum(covarOrderedPair$trainPair == TRUE))
  )
  
  
  #### Estimating the generation interval distribution ####
  
  siPars <- siMethods(allProbs, allProbs2, shift = shift,
                      initialPars = initialPars, dateVar = dateVar)
  
  
  #### Assessing the results ####
  
  #Finding the sensitivity of the gold standard - what percentage of the true links
  #are probable links
  sensitivity <- sum(allProbs$snpClose == TRUE & allProbs$transmission == TRUE, na.rm = TRUE) /
    sum(allProbs$transmission == TRUE, na.rm = TRUE)
  #Finding the ppv of the gold standard - what percentage of the probable links are true links
  ppv <- sum(allProbs$snpClose == TRUE & allProbs$transmission == TRUE, na.rm = TRUE) /
    sum(allProbs$snpClose == TRUE, na.rm = TRUE)
  
  #Finding the average number of links per case for those who have links in the training dataset
  trainLinks <- allProbs %>% filter(snpCloseGS == TRUE)
  linkData <- (allProbs
               %>% group_by(individualID.2)
               %>% summarize(nLinks = n())
  )
  posLinks <- mean(linkData$nLinks)
  
  #Finding the proportion of individuals included in training links
  pTrainingLinks <- length(unique(c(trainLinks$individualID.1,
                                    trainLinks$individualID.2))) / nrow(covarInd)
  
  
  #Finding the total time of the outbreak
  totalTime <- as.numeric(difftime(max(indData$infectionDate),
                                   min(indData$infectionDate),
                                   units = "days")) / 365
  
  
  siPars2 <- cbind(label = allProbs$label[1], nCases = nrow(indData),
                   nOutbreaks = length(unique(indData$outbreakID)),
                   pTraining, lowerT, pTrainingLinks,
                   nTrainLinks = nrow(trainLinks), totalTime,
                   obsMean, obsMedian, obsSD, sensitivity, ppv,
                   siPars, stringsAsFactors = FALSE)
  
  siPars2 <- siPars2 %>% mutate(meanDiff = meanSI - obsMean,
                                medianDiff = medianSI - obsMedian,
                                sdDiff = sdSI - obsSD)

  
  return(list(siPars2, pTemp))
} 



#### Function to calculate all of the different SI estimates ####

siMethods <- function(allProbs, allProbs2, shift, initialPars, dateVar){
  
  ## Constant Cutoff ##

  siParsN <- estimateSI(df = allProbs, indIDVar = "individualID", timeDiffVar = dateVar,
                        pVar = "pScaled",  clustMethod = "n", cutoffs = seq(1, 10, 1),
                        shift = shift, initialPars = initialPars)
  
  siParsN$prob <- ifelse(siParsN$clustMethod == "pooled", "Npooled",
                          paste0("N", siParsN$cutoff))

  print("Finished constant n")
  
  
  
  ## Hierarchical Clustering ##
  
  siParsHC <- estimateSI(df = allProbs, indIDVar = "individualID",
                         timeDiffVar = dateVar, pVar = "pScaled",
                         clustMethod = "hc_absolute", cutoffs = seq(0.025, 0.25, 0.025),
                         shift = shift, initialPars = initialPars,
                         bootSamples = bootSamples)
  
  siParsHC0 <- estimateSI(df = allProbs, indIDVar = "individualID",
                          timeDiffVar = dateVar, pVar = "pScaled",
                          clustMethod = "hc_absolute", cutoffs = 0,
                          shift = shift, initialPars = initialPars)
  
  siParsHC <- dplyr::bind_rows(siParsHC, siParsHC0)
  siParsHC$prob <- ifelse(siParsHC$clustMethod == "pooled", "HCpooled",
                          paste0("HC", siParsHC$cutoff))
  
  print("Finished hierarchical clustering")
  
  
  
  ## Kernel Density Clustering ##
  
  #Different binwidths and minimum gaps
  siParsKD <- estimateSI(df = allProbs, indIDVar = "individualID",
                         timeDiffVar = dateVar, pVar = "pScaled",
                         clustMethod = "kd", cutoffs = seq(0.01, 0.1, 0.01),
                         shift = shift, initialPars = initialPars,
                         bootSamples = bootSamples)
  
  siParsKD$prob <- ifelse(siParsKD$clustMethod == "pooled", "KDpooled",
                          paste0("KD", siParsKD$cutoff))
  
  #Including time
  siParsKDT <- estimateSI(df = allProbs2, indIDVar = "individualID",
                         timeDiffVar = dateVar, pVar = "pScaled",
                         clustMethod = "kd", cutoffs = seq(0.01, 0.1, 0.01),
                         shift = shift, initialPars = initialPars,
                         bootSamples = 0)
  
  siParsKDT$prob <- ifelse(siParsKD$clustMethod == "pooled", "KDTpooled",
                          paste0("KDT", siParsKD$cutoff))
  siParsKDT$clustMethod <- "kdt"
  
  print("Finished kernel density")
  

  
  ## Comparison Methods ##
  
  #Using just the training links
  trainLinks <- as.data.frame(allProbs %>% filter(snpCloseGS == TRUE))

  siParsT <- performPEM(df = trainLinks, indIDVar = "individualID",
                        timeDiffVar = dateVar, pVar = "pScaled",
                        initialPars = initialPars, shift = shift)
  
  siParsT <- siParsT %>% mutate(prob = "TrainLinks",
                                meanSI = shape * scale + shift,
                                medianSI = qgamma(0.5, shape = shape, scale = scale) + shift,
                                sdSI = sqrt(shape * scale^2),
                                nIndividuals = length(unique(trainLinks$individualID.2)))
  
  #Naive method
  naive <- optim(par = initialPars, logl_si, df = trainLinks, timeDiffVar = dateVar,
                 method = "L-BFGS-B", lower = c(0.0001, 0.0001), upper = c(Inf, Inf))
  siParsNa <- as.data.frame(t(naive$par))
  names(siParsNa) <- c("shape", "scale")
  siParsNa <- siParsNa %>% mutate(prob = "Naive",
                                  meanSI = shape * scale + shift,
                                  medianSI = qgamma(0.5, shape = shape, scale = scale) + shift,
                                  sdSI = sqrt(shape * scale^2),
                                  nIndividuals = length(unique(trainLinks$individualID.2)))
  
  print("Finished SNP distance method")

  #Using all pairs
  siParsA <- performPEM(df = allProbs, indIDVar = "individualID",
                        timeDiffVar = dateVar, pVar = "pScaled",
                        initialPars = initialPars, shift = shift)
  siParsA <- siParsA %>% mutate(prob = "All",
                                meanSI = shape * scale + shift,
                                medianSI = qgamma(0.5, shape = shape, scale = scale) + shift,
                                sdSI = sqrt(shape * scale^2),
                                nIndividuals = length(unique(allProbs$individualID.2)))
  print("Finished all pairs")
  
  #Combining methods
  siPars <- bind_rows(siParsA, siParsNa, siParsT, siParsN, siParsHC, siParsKD, siParsKDT)

  return(siPars)
}

