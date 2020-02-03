#Sarah V. Leavitt
#Boston University Dissertation
#Paper 2

#################################################################################
# This program reads in and formats the Mass DPH data to prepare for analysis
#################################################################################

setwd("~/Boston University/Dissertation/nbPaper2")
rm(list = ls())
set.seed(103020)

library(dplyr)
library(tidyr)
library(naniar)
library(lubridate)
library(reshape2)
library(readxl)
library(stringr)
library(tableone)
library(ggplot2)

#Reading in the original dataset
massInd <- read_excel("../Datasets/genotyping IRB data set-  jan 2010 to dec 2016 ver 2.5.xlsx",
                      na = c("", "#N/A"))

#Reading in my contact groupings
massContacts <- read.csv("../Datasets/DPH_ContactGroups.csv")

#Looking at resistance data
vars <- c("ISUSINH", "ISUSRIF", "ISUSPZA", "ISUSEMB", "ISUSSM", "ISUSETH", "ISUSKAN",
          "ISUSCYC", "ISUSCAP", "ISUSPAS", "ISUSAM", "ISUSRIB", "ISUSCIP", "ISUSOFL",
          "ISUSRPT", "ISUSLEVO", "ISUSMOXI", "ISUSQUIN")
CreateTableOne(vars = vars, data = massInd, factorVars = vars)
#Variables to use (enough resistance to be meaningful):
#ISUSPZA, ISUSRIF, ISUSPZA, ISUSEMB, ISUSSM, ISUSETH

#Looking at smear and culture data
table(massInd$SPSMEAR, massInd$Smearother)
table(massInd$SPCULT, massInd$Cultureother)

#Find out how many cases are involved in the contact tracing (36)
haveContact <- (massInd
                %>% filter(!is.na(`Linked Case1`))
                %>% select(StudyID, `Linked Case1`)
)

haveContactInvest <- (massInd
                      %>% filter(`Contact Investigation Done` == "Yes" |
                                   `Identified in Epi Investigation` == "Yes")
                      %>% select(StudyID, `Contact Investigation Done`,
                                 `Identified in Epi Investigation`, `Linked Case1`)
)



#Cleaning up the data
massInd2 <- (massInd
             %>% replace_with_na(list(ISUSINH = c("NOT", "UNK"), ISUSRIF = c("NOT", "UNK"),
                                      ISUSPZA = c("NOT", "UNK"), ISUSEMB = c("NOT", "UNK"),
                                      ISUSSM = c("NOT", "UNK"), ISUSETH = c("NOT", "UNK"),
                                      MIRU = c("NoResult", "NORESULT"),
                                      MIRU2 = c("NoResult", "NORESULT"),
                                      Spoligotype = c("NoResult", "NORESULT"),
                                      SEX = c("UNK"),
                                      `Genotyping Lineage` = c("antelope", "H37Rv/Ra")))
             %>% mutate(Smear = ifelse(SPSMEAR == "POS" | Smearother == "POS", "POS",
                                       ifelse(SPSMEAR == "NEG" | Smearother == "NEG", "NEG", NA)),
                        Culture = ifelse(SPCULT == "POS" | Cultureother == "POS", "POS",
                                         ifelse(SPCULT == "NEG" | Cultureother == "NEG", "NEG", NA)),
                        Age = ifelse(AGE3 %in% c("00-04", "05-14"), "<15", AGE3),
                        YearOfArrival = as.numeric(`Year of Arrival`),
                        USBorn = `Country of birth` == "UNITED STATES",
                        EpiLinkFound = ifelse(`Identified in Epi Investigation` == "0", "No",
                                              `Identified in Epi Investigation`))
             %>% extract(`Suspect Month-year`, into = c("SuspectM", "SuspectY"),
                         "([[:alnum:]]+)-([[:alnum:]]+)", remove = FALSE)
             %>% extract(`Counted Month-Year`, into = c("CountedM", "CountedY"),
                         "([[:alnum:]]+)-([[:alnum:]]+)", remove = FALSE)
             %>% mutate(SuspectDt = ymd(paste(SuspectY, SuspectM, "15", sep = "-")), "%Y-%m",
                        CountedDt = ymd(paste(CountedY, CountedM, "15", sep = "-")), "%Y-%m",
                        CombinedDt = as.Date(ifelse(SuspectDt < CountedDt, SuspectDt, CountedDt),
                                             origin = "1970-01-01"),
                        CombinedY = year(CombinedDt),
                        RecentArrival2 = as.numeric(CombinedY) - YearOfArrival <= 2,
                        RecentArrival1 = as.numeric(CombinedY) - YearOfArrival <= 1,
                        MIRUComb = paste0(MIRU, MIRU2))
             %>% replace_with_na(list(MIRUComb = "NANA"))
             %>% full_join(massContacts, by = "StudyID")
             %>% select(StudyID, SuspectDt, CountedDt, CombinedDt, CombinedY, County = COUNTY,
                        Sex = SEX, Age, Spoligotype, MIRU, MIRU2, MIRUComb, GENType, PCRType,
                        Lineage = `Genotyping Lineage`, CountryOfBirth = `Country of birth`,
                        USBorn, YearOfArrival, RecentArrival1, RecentArrival2,
                        ISUSINH, ISUSRIF, ISUSPZA, ISUSEMB, ISUSSM, ISUSETH, Smear, Culture,
                        AnyImmunoSup = `Any ImmunoSupression`, Contact = `Linked Case1`, ContactGroup,
                        HaveContInv = `Contact Investigation Done`, EpiLinkFound)
             %>% mutate(HaveContact = !is.na(ContactGroup))
             #Excluding cases with M. Bovis
             %>% filter(is.na(Lineage) | Lineage != "Bovis")
)

table(massInd2$CombinedY, massInd2$HaveContInv)

haveContact2 <- (massInd2
                 %>% filter(!is.na(ContactGroup) | !is.na(Contact))
                 %>% select(StudyID, ContactGroup, Contact, HaveContact, HaveContInv,
                            EpiLinkFound, RecentArrival1, RecentArrival2, YearOfArrival, CombinedY)
)

#### Creating Training Pairs ####

#Histogram of case counts
ggplot(data = massInd2) +
  geom_histogram(aes(x = CombinedDt))

#Histogram of cases with contact investigations
ggplot(data = massInd2 %>% filter(HaveContInv == "Yes" | HaveContact == TRUE)) +
  geom_bar(aes(x = factor(year(CombinedDt)), fill = HaveContact),
           position = "dodge")

#Randomly choosing 8 people per year who have contact investigations
#(so there is the same number of people each year)
contactSub <- (massInd2
               %>% filter(HaveContInv == "Yes", HaveContact == FALSE)
               %>% group_by(CombinedY)
               %>% sample_n(8)
)

massInd3 <- massInd2 %>% mutate(HaveContInvTrain = ifelse(StudyID %in% contactSub$StudyID |
                                                            HaveContact == TRUE, TRUE, FALSE))
table(massInd3$HaveContact, massInd3$HaveContInvTrain)

#Histogram of cases used in the training dataset
ggplot(data = massInd3 %>% filter(HaveContInvTrain == TRUE)) +
  geom_bar(aes(x = factor(year(CombinedDt)), fill = is.na(ContactGroup)),
           position = "dodge")


################### Creating a pairs dataframe #####################

#Finding all pairs of IDs (order matters)
pairs <- expand.grid(massInd3$StudyID, massInd3$StudyID)
pairs2 <- (pairs
           %>% rename(StudyID.1 = Var1,  StudyID.2 = Var2)
           %>% filter(StudyID.1 != StudyID.2)
           %>% mutate(EdgeID = paste(StudyID.1, StudyID.2, sep = "_"))
)

massPair <- (pairs2
             #Combining with the individual level data by StudyID.1
             %>% full_join(massInd3, by = c("StudyID.1" = "StudyID"))
             %>% full_join(massInd3, by = c("StudyID.2" = "StudyID"),
                           suffix = c(".1", ".2"))
             #Defining pair-level covariates
             %>% mutate(CombinedDiff = as.numeric(difftime(CombinedDt.2, CombinedDt.1, units = "days")),
                        CombinedDiffY = ifelse(CombinedDiff == 0, 1/365, CombinedDiff/ 365),
                        Lineage = Lineage.1 == Lineage.2,
                        Lineage = factor(Lineage, levels = c(FALSE, TRUE),
                                         labels = c("Different", "Same")),
                        Spoligotype = Spoligotype.1 == Spoligotype.2,
                        Spoligotype = factor(Spoligotype, levels = c(FALSE, TRUE),
                                             labels = c("Different", "Same")),
                        GENType = GENType.1 == GENType.2,
                        GENType = factor(GENType, levels = c(FALSE, TRUE),
                                         labels = c("Different", "Same")),
                        PCRType = PCRType.1 == PCRType.2,
                        PCRType = factor(PCRType, levels = c(FALSE, TRUE),
                                         labels = c("Different", "Same")),
                        Sex = ifelse(Sex.1 == "M" & Sex.2 == "M", 1,
                              ifelse(Sex.1 == "F" & Sex.2 == "F", 2,
                              ifelse(Sex.1 == "M" & Sex.2 == "F", 3, 4))),
                        Sex = factor(Sex, levels = c(2, 1, 3, 4),
                                     labels = c("f-f", "m-m", "m-f", "f-m")),
                        Smear = ifelse(Smear.1 == "POS", 2, 1),
                        Smear = factor(Smear, levels = c(1, 2),
                                       labels = c("InfectorSmear-", "InfectorSmear+")),
                        AnyImmunoSup = ifelse(AnyImmunoSup.1 == "Y", 2, 1),
                        AnyImmunoSup = factor(AnyImmunoSup, levels = c(1, 2),
                                       labels = c("InfectorImmuneSup-", "InfectorImmuneSup+")),
                        CountryOfBirth = ifelse(CountryOfBirth.1 == "UNITED STATES" &
                                                CountryOfBirth.2 == "UNITED STATES", 1,
                                         ifelse(CountryOfBirth.1 != "UNITED STATES" &
                                                CountryOfBirth.2 != "UNITED STATES" &
                                                CountryOfBirth.1 == CountryOfBirth.2, 2,
                                         ifelse(CountryOfBirth.1 != CountryOfBirth.2 &
                                                 (CountryOfBirth.1 == "UNITED STATES" |
                                                  CountryOfBirth.2 == "UNITED STATES"), 3, 4))),
                        CountryOfBirth = factor(CountryOfBirth, levels = c(3, 1, 4, 2),
                                             labels = c("Diff-USA","Same-USA", "Diff-Other",
                                                        "Same-Other")),
                        Age = Age.1 == Age.2,
                        Age = factor(Age, levels = c(FALSE, TRUE),
                                     labels = c("Different", "Same")),
                        TimeCat = ifelse(CombinedDiffY <= 1, 1,
                                  ifelse(CombinedDiffY > 1 & CombinedDiffY <= 2, 2,
                                  ifelse(CombinedDiffY > 2 & CombinedDiffY <= 3, 3,
                                  ifelse(CombinedDiffY > 3 & CombinedDiffY <= 4, 4, 5)))),
                        TimeCat = factor(TimeCat, levels = c(1, 2, 3, 4, 5),
                                         labels = c("<=1y", "1-2y", "2-3y", "3-4y", ">4y")),
                        ISUSINH = ifelse(ISUSINH.1 == "R" & ISUSINH.2 == "R", 1,
                                  ifelse(ISUSINH.1 == "S" & ISUSINH.2 == "S", 2,
                                  ifelse(ISUSINH.1 == "R" & ISUSINH.2 == "S", 3, 4))),
                        ISUSINH = factor(ISUSINH, levels = c(1, 2, 3, 4),
                                         labels = c("R-R", "S-S", "R-S", "S-R")),
                        ISUSRIF = ifelse(ISUSRIF.1 == "R" & ISUSRIF.2 == "R", 1,
                                  ifelse(ISUSRIF.1 == "S" & ISUSRIF.2 == "S", 2,
                                  ifelse(ISUSRIF.1 == "R" & ISUSRIF.2 == "S", 3, 4))),
                        ISUSRIF = factor(ISUSRIF, levels = c(1, 2, 3, 4),
                                         labels = c("R-R", "S-S", "R-S", "S-R")),
                        ISUSPZA = ifelse(ISUSPZA.1 == "R" & ISUSPZA.2 == "R", 1,
                                  ifelse(ISUSPZA.1 == "S" & ISUSPZA.2 == "S", 2,
                                  ifelse(ISUSPZA.1 == "R" & ISUSPZA.2 == "S", 3, 4))),
                        ISUSPZA = factor(ISUSPZA, levels = c(1, 2, 3, 4),
                                         labels = c("R-R", "S-S", "R-S", "S-R")),
                        ISUSEMB = ifelse(ISUSEMB.1 == "R" & ISUSEMB.2 == "R", 1,
                                  ifelse(ISUSEMB.1 == "S" & ISUSEMB.2 == "S", 2,
                                  ifelse(ISUSEMB.1 == "R" & ISUSEMB.2 == "S", 3, 4))),
                        ISUSEMB = factor(ISUSEMB, levels = c(1, 2, 3, 4),
                                         labels = c("R-R", "S-S", "R-S", "S-R")),
                        ISUSSM = ifelse(ISUSSM.1 == "R" & ISUSSM.2 == "R", 1,
                                 ifelse(ISUSSM.1 == "S" & ISUSSM.2 == "S", 2,
                                 ifelse(ISUSSM.1 == "R" & ISUSSM.2 == "S", 3, 4))),
                        ISUSSM = factor(ISUSSM, levels = c(1, 2, 3, 4),
                                         labels = c("R-R", "S-S", "R-S", "S-R")),
                        ISUSETH = ifelse(ISUSETH.1 == "R" & ISUSETH.2 == "R", 1,
                                  ifelse(ISUSETH.1 == "S" & ISUSETH.2 == "S", 2,
                                  ifelse(ISUSETH.1 == "R" & ISUSETH.2 == "S", 3, 4))),
                        ISUSETH = factor(ISUSETH, levels = c(1, 2, 3, 4),
                                         labels = c("R-R", "S-S", "R-S", "S-R")),
                        ContactGroup = ContactGroup.1 == ContactGroup.2,
                        HaveContInvTrain = HaveContInvTrain.1 == TRUE & HaveContInvTrain.2 == TRUE,
                        ContactTrain = ifelse(!is.na(ContactGroup) & ContactGroup == TRUE, TRUE,
                                       ifelse(HaveContInvTrain == TRUE, FALSE, NA)))
             %>% unite(Resistance, ISUSINH, ISUSRIF, ISUSPZA, ISUSEMB, ISUSSM, ISUSETH,
                       remove = FALSE)
             %>% mutate(SharedRes = str_count(Resistance, "R-R"),
                        SharedResG = ifelse(SharedRes >= 3, "3+", as.character(SharedRes)),
                        SharedResG = factor(SharedResG, levels = c("0", "1", "2", "3+")))
)

table(massPair$ContactTrain, useNA = "always")

#Adding county data
table(massInd2$County)
countyMatrix <- read.csv("../Datasets/DPH_Counties.csv", row.names = "RowNames",
                         stringsAsFactors = FALSE)
countyDf <- reshape2::melt(as.matrix(countyMatrix), varnames = c("County.1", "County.2"))
names(countyDf) <- c("County.1", "County.2", "County")
countyDf$County.1 <- as.character(countyDf$County.1)
countyDf$County.2 <- as.character(countyDf$County.2)

massPair2 <- (massPair
              %>% left_join(countyDf, by = c("County.1", "County.2"))
              %>% mutate(County = factor(County, levels = c(3, 2, 1),
                                         labels = c("Other", "Neighbor", "Same")))
  
)
table(massPair2$County, useNA = "ifany")


#Function to find the difference in MIRU
#Not including - (error) and % (mixed result) xyz (some other anomalous reading)
findDiff <- function(a, b){
  
  #Replacing characters meaning missing values with #
  a <- gsub("[%]|[-]|[x]|[y]|[z]", "#", a)
  b <- gsub("[%]|[-]|[x]|[y]|[z]", "#", b)
  
  #Splitting the strings into a list of two vectors
  strList <- strsplit(c(a, b), split = "")
  
  #Extracting the different locations in the first and second MIRU
  diffa <- strList[[1]][strList[[1]] != strList[[2]]]
  diffb <- strList[[2]][strList[[1]] != strList[[2]]]
  
  #Removing # signs representing unknown locations
  diffa <- diffa[diffa != "#"]
  diffb <- diffb[diffb != "#"]
  
  diffLength <- min(length(diffa), length(diffb))
  
  return(diffLength)
}

#Adding MIRU difference
massPair3 <- (massPair2
              %>% rowwise()
              %>% mutate(MIRUDiff = findDiff(MIRUComb.1, MIRUComb.2))
              %>% ungroup()
              %>% mutate(MIRUDiffG = ifelse(MIRUDiff >= 4, "4+", as.character(MIRUDiff)),
                         MIRUDiffG = factor(MIRUDiffG, levels = c("0", "1", "2", "3", "4+")))
)

#Saving datasets
saveRDS(massInd3, "../Datasets/MassInd.rds")
saveRDS(massPair3, "../Datasets/MassPair.rds") 


  