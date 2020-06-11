
library (plyr)
library (dplyr)
library (purrr)
library (tidyverse)
library (ggplot2)
library(matrixStats) 
library (plotly)


ControlMerit <- read.csv("CommMeritSumC.csv") 
ControlIAVcost <- read.csv("IAVcostSumC.csv") 

########################################################

productivity <- 6.6  ### loss in productiity from IAV 
vaccination <- 3.71  ### estimated vaccination costs IAV
GenImprov <- 4  #### genetic improvement annually - pigsite with PIC
HerdImmunity <- 90 ### level of herd immunity to stop vaccinating
InterestRate <- 0.05 

EditingCost <- 10
ZygoteDeathRate <- 0.25 #inverse of zygote survival
MarketValue <- 109.5 #rough USD value of carcass at slaughter. 10 year average
##### CommPop count must be compiled into table

####################################

MeritOutput <- read.csv("CommMeritSumDGtEP.csv")
PropOutput <- read.csv("CommPropSumDGtEP.csv")
EditsOutput <- read.csv("SPF_EditCountSumDGtEP.csv")
IAVcostOutput <- read.csv("IAVcostSumDGtEP.csv") 
CommPop <-read.csv("CommPopDGtEP.csv")
CommResistant <- read.csv("CommResistantDGtEP.csv")
##################################

EconAnalysis <- EconOutput(CommPop, CommResistant, MeritOutput, PropOutput, EditsOutput, IAVcostOutput, ControlMerit, ControlIAVcost, vaccination, productivity, GenImpov, HerdImmunity, InterestRate, EditingCost)
write.csv(EconAnalysis, "EconAnalysis_DGt_EP.csv", row.names = FALSE)

  
EconOutput <- function(CommPop, CommResistant, MeritOutput, PropOutput, EditsOutput, IAVcostOutput, ControlMerit, ConrolIAVcost, vaccination, productivity, GenImpov, HerdImmunity, InterestRate, EditingCost){

  MergeforEcon <- as.data.frame(cbind(MeritOutput$gen, CommPop$MeanAll, CommResistant$MeanAll, MeritOutput$MeanAll, PropOutput$MeanAll ,EditsOutput$MeanAll, IAVcostOutput$MeanAll, ControlMerit$MeanAll, ControlIAVcost$MeanAll))
  MergeforEcon <- setNames(MergeforEcon, c("gen", "CommPop", "CommResistant", "Merit", "Prop", "Edits", "IAV", "MeritCon", "IAVcon"))
  

  EconAnalysis <- setNames(as.data.frame(matrix(nrow = 10, ncol =9)), c("year", "EditCost", "LostPigletsValue", "MeritCost", "Health", "RealValue", "DiscountFactor", "PresentValue", "CumPv"))
  EconAnalysis$year <- seq(1,10)
  EconAnalysis$DiscountFactor <- (1/(1+InterestRate))^EconAnalysis$year
 
  ###### Create column of year ###
  MergeforEcon$year <- ifelse(MergeforEcon$gen == 1:12, 1, 
                           ifelse(MergeforEcon$gen == 13:24, 2, 
                                  ifelse(MergeforEcon$gen == 25:36, 3, 
                                         ifelse(MergeforEcon$gen == 37:48, 4, 
                                                ifelse(MergeforEcon$gen == 49:60, 5, 
                                                       ifelse(MergeforEcon$gen == 61:72, 6, 
                                                              ifelse(MergeforEcon$gen == 73:84, 7, 
                                                                     ifelse(MergeforEcon$gen == 85:96, 8, 
                                                                            ifelse(MergeforEcon$gen == 97:108, 9, 
                                                                                   ifelse(MergeforEcon$gen == 109:120, 10, 0))))))))))
  
  
  #### Count cost of editing per year, cummulative of each year group figure
  SumEdits <- MergeforEcon %>%  group_by(year) %>% summarise (EditCount = sum(Edits))
  EconAnalysis$EditCost <- SumEdits$EditCount * EditingCost
  EconAnalysis$LostPigletsValue <- SumEdits$EditCount * ZygoteDeathRate * MarketValue

  #### Genetic Merit Cost ####
  
  MergeforEconMerit <- MergeforEcon %>% 
    select (Merit, gen) %>%
    filter (gen %in% c(6, 18, 30, 42, 54, 66, 78, 90, 102, 114))
  
  ControlMerit <- ControlMerit %>% 
    select (MeanAll, gen) %>%
    filter (gen %in% c(6, 18, 30, 42, 54, 66, 78, 90, 102, 114))
  
  Z <- as.data.frame(((MergeforEconMerit$Merit - ControlMerit$MeanAll) / ControlMerit$MeanAll)) 
  
  CommPopAve <- MergeforEcon %>%  group_by(year) %>% summarise (CommPopAve = sum(CommPop)) #### need commpopave
  
  EconAnalysis$MeritCost <-  Z[,1] * (GenImprov * EconAnalysis$year) * CommPopAve$CommPopAve ### need to deterMIne the actual piglets born per year better 
  
  ################### Health Gain ###
  
  MergeforEcon$Health <- ifelse(MergeforEcon$Prop > HerdImmunity, MergeforEcon$CommPop * vaccination, 0)
  Health <- MergeforEcon %>%  group_by(year) %>% summarise (Health = sum(Health)) %>% select(Health = Health)
  EconAnalysis$Health <- Health$Health

 
  EconAnalysis$RealValue <- EconAnalysis$Health - EconAnalysis$EditCost + EconAnalysis$MeritCost - EconAnalysis$LostPigletsValue # Health is gain, Editing is cost, Merit is gain or cost
  
  EconAnalysis$PresentValue <- EconAnalysis$RealValue * EconAnalysis$DiscountFactor
  
  EconAnalysis$CumPv <- cumsum(EconAnalysis$PresentValue)
  
  
#  write.csv(EconAnalysis, paste0("EconAnalysis_",Output,".csv"))
  
  return (EconAnalysis) }


