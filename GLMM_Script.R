library(MCMCglmm)
library(coda)#gelman convergence diagnostic
library(imputeTS) # for na.replace function
library(tidyverse)
library(ggplot2)
library(ggregplot)
library(lme4)
library(postMCMCglmm)

rm(list = ls())

# Import and clean data ####
dt <- read.csv("Network_summary_master_file.csv", header=T)

dt <- dt[,c("filename","study_id","transmission_route","class","sociality", "edge_wt_type_reduced", "scale", "data_collection_duration_unit_day",
            "nodes","edges"
            ,"network.density","CV.degree",
            "deg.assort","avg.betw.centrality",
            "clustering", "cohesion", "num.modules","diameter")]
dt$data_collection_duration_unit_day <- as.numeric(dt$data_collection_duration_unit_day)
dt$data_collection_duration_unit_day <- na_replace(dt$data_collection_duration_unit_day,0)
dt<- subset(dt, dt$network.density>0) #only include networks with edges
dt<-subset(dt, dt$network.density<1)  #only data where density is less than 1 (self loops)

#make sure transmission modes are factors
dt$transmission_route.dv <- factor(dt$transmission_route, levels = c( "Nonphysical", "Physical", "Fluid", "Indirect")) ##making 
summary(factor(dt$transmission_route.dv))

###scaling and log transforming all the network metrics to encourage model fitting
dt$cnodes <- scale(dt$nodes, center=TRUE, scale=TRUE)
dt$cdata_collection_duration_unit_day <- scale(dt$data_collection_duration_unit_day, center=TRUE, scale=TRUE)
dt$clog.CV.degree <- scale(log(dt$CV.degree+1), center=TRUE, scale=TRUE)
dt$clog.network.density <- scale(log(dt$network.density), center=TRUE, scale=TRUE)
dt$clog.betw.centrality <- scale(log(dt$avg.betw.centrality+1), center=TRUE, scale=TRUE)
dt$clog.diameter <- scale(log(dt$diameter+2), center=TRUE, scale=TRUE)
dt$clog.clustering <- scale(log(dt$clustering+1), center=TRUE, scale=TRUE)
dt$clog.num.modules <- scale(log(dt$num.modules), center=TRUE, scale=TRUE)
dt$clog.deg.assort <- scale(log(dt$deg.assort+2), center=TRUE, scale=TRUE)
dt$clog.cohesion <- scale(log(dt$cohesion+2), center=TRUE, scale=TRUE)

####convert all sociality levels to binary
dt$sociality.dv <- ifelse(dt$sociality=="Relativelysolitary", "solitary",
                          ifelse(dt$sociality=="Fissionfusion", "social",
                                 ifelse(dt$sociality=="Social","social", "socialDom")))
dt$sociality.dv <- factor(dt$sociality.dv, levels = c( "solitary",  "social", "socialDom"))
summary(factor(dt$sociality.dv))

##Convert study scales to binary factors
dt$scale.dv <- ifelse(dt$scale== "captive", "cap",
                      ifelse(dt$scale== "social", "soc", "spat"))
dt$scale <- factor(dt$scale, levels = c("captive", "social", "spatial" ))


dt<-na.omit(dt) ##remove networks with Nan values for any metric

#############################################################################
# Multivariate GLMM model ####
FixedCovar1 <- c("transmission_route.dv",
                 "cdata_collection_duration_unit_day",
                 "edge_wt_type_reduced",
                 "scale", "sociality.dv",
                 "cnodes")    #Fixed Predictor variables

Resps1 <- c(
  "clog.CV.degree",
  "clog.betw.centrality",
  "clog.deg.assort",
  "clog.cohesion",
  "clog.num.modules",
  "clog.clustering",
  "clog.diameter", "clog.network.density"
) #Response Varaibles

Resps1<-cbind(Resps1)

Resps1 %>% c -> Resps1

MultivFormula_Base <- as.formula(paste0("cbind(", paste(Resps1, collapse = ","), ")", "~ - 1 + trait + trait:(transmission_route.dv +", "(",
                                        paste(FixedCovar1[-1], collapse = " + "), "))"))
MultivPrior1<-list(R=list(V=diag(length(Resps1)), nu=9.002),
                   G=list(G1=list(V=diag(length(Resps1)), nu=length(Resps1),
                                  alpha.mu=rep(0,length(Resps1)),
                                  alpha.V=diag(length(Resps1))*100)))

ModelList <- list()

N = 1000

x <- 1

for(x in x:N){ # for loop to perform 1000 GLMMs pulling max 15 networks from each study
  
  print(x)
  
  shd = data.frame()
  
  for(i in study_id_list){
    
    NID = as.numeric(nrow(subset(dt, study_id == i)))
    
    df <- if(NID >15){
      
      sample_n(subset(dt, study_id == i), 15)
      
    }
    
    else {
      
      subset(dt, study_id == i)
      
    }
    
    shd <- rbind(shd,df)
    
  }
  
  
  shd %>% dplyr::select(Resps1) %>% mutate_all(as.numeric)  ->
    
    sd1
  
  
  sd1 %>%
    
    bind_cols(dplyr::select(shd, -c(colnames(sd1)))) %>%
    
    mutate(transmission_route.dv = factor(transmission_route.dv, levels = c("Nonphysical",
                                                                            "Physical",
                                                                            "Fluid",
                                                                            "Indirect"))) -> sd1
  
  mod <- MCMCglmm(
    
    fixed = MultivFormula_Base,
    random =~ idh(trait):study_id,
    data = sd1,
    prior = MultivPrior1,
    rcov =~ us(trait):units,
    family = rep("gaussian", length(Resps1)),
    verbose = F,
    nitt = 10500, thin = 10, burnin = 500, DIC = TRUE)
  
  ModelList[[x]] <- mod
  
}

saveRDS(ModelList, "Model_List_1000_runs.rds")

# Visualise Results ####

ModelList <-readRDS("Model_List_1000_runs.rds")

#Function that will take one sample from each of the 1000 models run
RandomSlice <- function(a, Max = NULL){
  if(is.null(Max)){
    Max <- nrow(a)
  }
  
  if((Max>nrow(a))){
    Max <- nrow(a)
  }
  Sample <- sample(1:nrow(a), Max)
  slice(a, Sample)
}

#pull one sample from each of the 1000 models and store in a dataframe named SliceEstimates
ModelList %>% map_dfr(~{
  .x$Sol %>% as.data.frame() %>% RandomSlice(1)
}, .id = "Model") -> SliceDF

SliceDF %>% dplyr::select(-Model) %>% apply(2, function(a){
  a %>% as.mcmc -> MCMCa
  data.frame(
    Mode = posterior.mode(MCMCa),
    Lower = HPDinterval(MCMCa)[1],
    Upper = HPDinterval(MCMCa)[2]
  )
}) %>% bind_rows(.id = "Variable") -> SliceEstimates

#Create columns for transmission modes comparison
a1 <- paste0("trait",Resps1[1],":",FixedCovar1[1], c("Physical","Fluid", "Indirect"))
b1 <- paste0("trait",Resps1[2],":",FixedCovar1[1], c("Physical","Fluid", "Indirect"))
c1<- paste0("trait",Resps1[3],":",FixedCovar1[1], c("Physical","Fluid", "Indirect"))
d1<- paste0("trait",Resps1[4],":",FixedCovar1[1], c("Physical","Fluid", "Indirect"))
e1<- paste0("trait",Resps1[5],":",FixedCovar1[1], c("Physical","Fluid", "Indirect"))
f1<- paste0("trait",Resps1[6],":",FixedCovar1[1], c("Physical","Fluid", "Indirect"))
g1 <- paste0("trait",Resps1[7],":",FixedCovar1[1], c("Physical","Fluid", "Indirect"))
h1 <- paste0("trait",Resps1[8],":",FixedCovar1[1], c("Physical","Fluid", "Indirect"))

selected_columns1 <- c(a1, b1, c1, d1, e1, f1,g1,h1)
drop_columns1 <- paste0("trait",c(Resps1[1:8]), ":",FixedCovar1[1],"Nonphysical")

all_columns1<-c(selected_columns1, drop_columns1)

# Create function to calculate P values for differences in network metrics between transmission modes
PCalc <- function(Vector){
  
  table(Vector>0)
  min(table(Vector>0))/(length(Vector)/2)
  
}

#Export pvalues into csv file

SliceDF %>% as.data.frame() %>%
  dplyr::select(selected_columns1) ->
  
  SubSolSlice

SubSolSlice[,paste0("trait",c(Resps1[1:9]),":", FixedCovar1[1],"Nonphysical")] <- 0

SubSolSlice %>% dplyr::select(all_columns1) -> SubSolSlice

sapply(1:(ncol(SubSolSlice)), function(b){
  
  SubSolSlice %>% pull(b) -> Cols1
  
  (SubSolSlice - Cols1) %>% apply(2, PCalc)
  
}) -> SubListSlice
colnames(SubListSlice) <- rownames(SubListSlice) %>% str_remove("transmission_route.dv")
rownames(SubListSlice) <- rownames(SubListSlice) %>% str_remove("transmission_route.dv")

SubListSlice[SubListSlice==2] <- 0 # two-tailed test

pvalues <- as.data.frame(SubListSlice)

write_csv(pvalues, "pvalues.csv")

