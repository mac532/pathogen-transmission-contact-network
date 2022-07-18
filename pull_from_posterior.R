# libraries ####
library(lsmeans)
library(gridExtra)
library(grid)
# library(plyr)
library(lattice)
library(car) # for vif
library(MASS)
library(ape)
library(caper)
library(phytools)
library(MCMCglmm)
library(coda)#gelman convergence diagnostic
library(imputeTS) # for na.replace function
library(tidyverse)
library(INLA)
library(ggplot2)
library(ggregplot)
library(magrittr)
library(kader)
library(GGally)
library(ggforce)

## Import response and predictor variables
Resps1 <- c(
  "clog.CV.degree",
  "clog.betw.centrality",
  "clog.deg.assort",
  "clog.cohesion",
  "clog.num.modules",
  "clog.clustering",
  "clog.diameter", "clog.network.density"
)

Resps1<-cbind(Resps1)

FixedCovar1 <- c("transmission_route.dv",
                 "cdata_collection_duration_unit_day",
                 "edge_wt_type_reduced",
                 "scale", "sociality.dv",
                 "cnodes"
                 #,"cedges", "cavg.edge.strength"
)


##Import Model list from 1000 GLMMs

ModelList <-readRDS("Model_List_1000_runs.rds")

#Randmom Slice function to create SLiceDF with one posterior estimate per 1000 models
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

## Generate Slice DF
ModelList %>% map_dfr(~{
  
  .x$Sol %>% as.data.frame() %>% RandomSlice(1)
  
}, .id = "Model") -> SliceDF

Model <- ModelList[[1]]

#########################################################################

#Get density estimates from posteriors of 1000 model run ####

#posterior estimates
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
#model matrix
X <-Model$X

#Just want last 232 values of the matrix for the density response variable
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
densX <- tail(XDF, 232)

N <-1000 #1000 MCMC rows
i <- 1
#get predicted intercept densities (for nonphysical transmission)

nonphys_dens <- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(densX)) #first row of Sol x model matrix
  d <- mean(y) #mean of those values is one predicted estimate
  nonphys_dens[i] <- d #add it to the vector for droplet densities
  
}

nonphys_dens

###Now get the predicted densities for the other three transmission modes

a <- paste0("trait",Resps1[8],":",FixedCovar1[1], c("Physical","Fluid", "Indirect"))

Sol %>% as.data.frame() %>%
  dplyr::select(a) ->
  
  SubSol

phys_effects <- as.numeric(SubSol$`traitclog.network.density:transmission_route.dvPhysical`)
fluid_effects <- as.numeric(SubSol$`traitclog.network.density:transmission_route.dvFluid`)
ind_effects <- as.numeric(SubSol$`traitclog.network.density:transmission_route.dvIndirect`)

#get physical transmission densities...
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
densX <- tail(XDF, 232)

N <-1000
i <- 1
phys_densities<- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(densX))
  d <- mean(y) + phys_effects[i]
  phys_densities[i] <- d
  
}

phys_densities

#get fluid exchange transmission densities...
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
densX <- tail(XDF, 232)

N <-1000
i <- 1
fluid_densities<- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(densX))
  d <- mean(y) + fluid_effects[i]
  fluid_densities[i] <- d
  
}

fluid_densities

#get indirect transmission densities
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
densX <- tail(XDF, 232)

N <-1000
i <- 1
ind_densities<- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(densX))
  d_list <- y + ind_effects[i]
  d = mean(d_list)
  ind_densities[i] <- d
  
}

ind_densities

Densities <- c(nonphys_dens, phys_densities, fluid_densities, ind_densities)


# For deg het ####
#Just want deg het values from the model matrix
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
degX <- subset(XDF, traitclog.CV.degree==1)

#get predicted intercept deeg het (for nonphysical transmission)

N <-1000 #1000 MCMC rows
i <- 1
nonphys_deg <- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(degX)) #first row of Sol x model matrix
  d <- mean(y) #mean of those values is one predicted estimate
  nonphys_deg[i] <- d 
}

nonphys_deg

#get effects for other three transmission modes
a <- paste0("trait",Resps1[1],":",FixedCovar1[1], c("Physical","Fluid", "Indirect"))

Sol %>% as.data.frame() %>%
  dplyr::select(a) ->
  
  SubSol

phys_effects <- as.numeric(SubSol$`traitclog.CV.degree:transmission_route.dvPhysical`)
fluid_effects <- as.numeric(SubSol$`traitclog.CV.degree:transmission_route.dvFluid`)
ind_effects <- as.numeric(SubSol$`traitclog.CV.degree:transmission_route.dvIndirect`)


##For Physical Degree Het
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
degX <- subset(XDF, traitclog.CV.degree==1)


N <-1000
i <- 1
phys_deg<- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(degX))
  d <- mean(y) + phys_effects[i]
  phys_deg[i] <- d
  
}

phys_deg

#For fluid-exchange deg het
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
degX <- subset(XDF, traitclog.CV.degree==1)

N <-1000 #1000 MCMC rows
i <- 1
fluid_deg<- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(degX))
  d <- mean(y) + fluid_effects[i]
  fluid_deg[i] <- d
  
}

fluid_deg

#For indirect deg het
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
degX <- subset(XDF, traitclog.CV.degree==1)

N <-1000 #1000 MCMC rows
i <- 1
ind_deg<- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(degX))
  d_list <- y + ind_effects[i]
  d = mean(d_list)
  ind_deg[i] <- d
  
}

ind_deg

Deg_Het <- c(nonphys_deg, phys_deg, fluid_deg, ind_deg)



# For betw ####

Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
betwX <- subset(XDF, traitclog.betw.centrality==1)

#get predicted intercept betw (for non physical transmission)
N <-1000 #1000 MCMC rows
i <- 1
nonphys_betw <- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(betwX)) #first row of Sol x model matrix
  d <- mean(y) #mean of those values is one predicted estimate
  nonphys_betw[i] <- d 
  
}

nonphys_betw


#get effects for other three transmission modes
a <- paste0("trait",Resps1[2],":",FixedCovar1[1], c("Physical","Fluid", "Indirect"))

Sol %>% as.data.frame() %>%
  dplyr::select(a) ->
  
  SubSol

phys_effects <- as.numeric(SubSol$`traitclog.betw.centrality:transmission_route.dvPhysical`)
fluid_effects <- as.numeric(SubSol$`traitclog.betw.centrality:transmission_route.dvFluid`)
ind_effects <- as.numeric(SubSol$`traitclog.betw.centrality:transmission_route.dvIndirect`)


##for physical betw
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
betwX <- subset(XDF, traitclog.betw.centrality==1)

N <-1000
i <- 1
phys_betw<- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(betwX))
  d <- mean(y) + phys_effects[i]
  phys_betw[i] <- d
  
}

phys_betw

#for fluid betw
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
betwX <- subset(XDF, traitclog.betw.centrality==1)

N <-1000 #1000 MCMC rows
i <- 1
fluid_betw<- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(betwX))
  d <- mean(y) + fluid_effects[i]
  fluid_betw[i] <- d
  
}

fluid_betw

#For indirect betw
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
betwX <- subset(XDF, traitclog.betw.centrality==1)

N <-1000 #1000 MCMC rows
i <- 1
ind_betw<- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(betwX))
  d_list <- y + ind_effects[i]
  d = mean(d_list)
  ind_betw[i] <- d
  
}

ind_betw

Betw <- c(nonphys_betw, phys_betw, fluid_betw, ind_betw)

# For deg Assort ####

Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
assrtX <- subset(XDF, traitclog.deg.assort==1)

N <-1000 #1000 MCMC rows
i <- 1
#get predicted intercept deg assrt for nonphysical transmission

nonphys_assrt <- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(assrtX)) #first row of Sol x model matrix
  d <- mean(y) #mean of those values is one predicted estimate
  nonphys_assrt[i] <- d 
  
}

nonphys_assrt


#get effects for other three transmission modes
a <- paste0("trait",Resps1[3],":",FixedCovar1[1], c("Physical","Fluid", "Indirect"))
Sol %>% as.data.frame() %>%
  dplyr::select(a) ->
  
  SubSol

phys_effects <- as.numeric(SubSol$`traitclog.deg.assort:transmission_route.dvPhysical`)
fluid_effects <- as.numeric(SubSol$`traitclog.deg.assort:transmission_route.dvFluid`)
ind_effects <- as.numeric(SubSol$`traitclog.deg.assort:transmission_route.dvIndirect`)


##for physical deg assrt
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
assrtX <- subset(XDF, traitclog.deg.assort==1)

N <-1000 #1000 MCMC rows
i <- 1
#get predicted intercept densities (for droplet)
phys_assrt<- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(assrtX))
  d <- mean(y) + phys_effects[i]
  phys_assrt[i] <- d
  
}

phys_assrt

#For fluid exchange deg assrt
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
assrtX <- subset(XDF, traitclog.deg.assort==1)

N <-100
i <- 1
fluid_assrt<- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(assrtX))
  d <- mean(y) + fluid_effects[i]
  fluid_assrt[i] <- d
  
}

fluid_assrt

#For indirect deg assort
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
assrtX <- subset(XDF, traitclog.deg.assort==1)

N <-1000 
i <- 1
ind_assrt<- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(assrtX))
  d_list <- y + ind_effects[i]
  d = mean(d_list)
  ind_assrt[i] <- d
  
}

ind_assrt

Assrt <- c(nonphys_assrt, phys_assrt, fluid_assrt, ind_assrt)

# For clustering ####

Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
clustX <- subset(XDF, traitclog.clustering==1)

N <-1000 #1000 MCMC rows
i <- 1
#get predicted intercept clustering (for non physical transmission)

nonphys_clust <- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(clustX)) #first row of Sol x model matrix
  d <- mean(y) #mean of those values is one predicted estimate
  nonphys_clust[i] <- d 
  
}

nonphys_clust


#get effects of other three transmission modes
a <- paste0("trait",Resps1[6],":",FixedCovar1[1], c("Physical","Fluid", "Indirect"))
Sol %>% as.data.frame() %>%
  dplyr::select(a) ->
  
  SubSol
phys_effects <- as.numeric(SubSol$`traitclog.clustering:transmission_route.dvPhysical`)
fluid_effects <- as.numeric(SubSol$`traitclog.clustering:transmission_route.dvFluid`)
ind_effects <- as.numeric(SubSol$`traitclog.clustering:transmission_route.dvIndirect`)


##for physical clustering
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
clustX <- subset(XDF, traitclog.clustering==1)

N <-1000 #1000 MCMC rows
i <- 1
phys_clust<- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(clustX))
  d <- mean(y) + phys_effects[i]
  phys_clust[i] <- d
  
}

phys_clust

#get fluid exchange clustering
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
clustX <- subset(XDF, traitclog.clustering==1)

N <-1000 #1000 MCMC rows
i <- 1
fluid_clust<- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(clustX))
  d <- mean(y) + fluid_effects[i]
  fluid_clust[i] <- d
  
}

fluid_clust

#get indirect clustering
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
clustX <- subset(XDF, traitclog.clustering==1)

N <-1000 #1000 MCMC rows
i <- 1
ind_clust<- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(clustX))
  d_list <- y + ind_effects[i]
  d = mean(d_list)
  ind_clust[i] <- d
  
}

ind_clust

Clust <- c(nonphys_clust, phys_clust, fluid_clust, ind_clust)

# Diameter ####

Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
diamX <- subset(XDF, traitclog.diameter==1)

N <-1000 #1000 MCMC rows
i <- 1
#get predicted intercept diameters for nonphysical transmission

nonphys_diam <- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(diamX)) #first row of Sol x model matrix
  d <- mean(y) #mean of those values is one predicted estimate
  nonphys_diam[i] <- d 
  
}

nonphys_diam

a <- paste0("trait",Resps1[7],":",FixedCovar1[1], c("Physical","Fluid", "Indirect"))

#get effects of each transmission mode
Sol %>% as.data.frame() %>%
  dplyr::select(a) ->
  
  SubSol

phys_effects <- as.numeric(SubSol$`traitclog.diameter:transmission_route.dvPhysical`)
fluid_effects <- as.numeric(SubSol$`traitclog.diameter:transmission_route.dvFluid`)
ind_effects <- as.numeric(SubSol$`traitclog.diameter:transmission_route.dvIndirect`)


##For physical diameter
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
diamX <- subset(XDF, traitclog.diameter==1)

N <-1000 #1000 MCMC rows
i <- 1
phys_diam<- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(diamX))
  d <- mean(y) + phys_effects[i]
  phys_diam[i] <- d
  
}

phys_diam

#get fluid diameters
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)
diamX <- subset(XDF, traitclog.diameter==1)

N <-1000 #1000 MCMC rows
i <- 1
fluid_diam<- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(diamX))
  d <- mean(y) + fluid_effects[i]
  fluid_diam[i] <- d
  
}

fluid_diam

#get indirect diameters
Sol <- subset(SliceDF, select = -c(Model) )
Sol<- as.matrix(Sol)
Sol <- as.mcmc(Sol)
X <- Model$X
X_mat <- as.matrix(X)
XDF <- as.data.frame(X_mat)

N <-1000 #1000 MCMC rows
i <- 1
ind_diam<- c()
for(i in i:N){
  y <-as.numeric((Sol[i,]) %*% t(diamX))
  d_list <- y + ind_effects[i]
  d = mean(d_list)
  ind_diam[i] <- d
  
}

ind_diam

Diam <- c(nonphys_diam, phys_diam, fluid_diam, ind_diam)


# Write dataframe of predictions ####
nplist <- replicate(1000,"Nonphysical")
physlist <- replicate(1000,"Physical")
indlist <- replicate(1000,"Indirect")
fluidlist <- replicate(1000,"Fluid-Exchange")
Tmode <- c(nplist, physlist, fluidlist, indlist)

DF_pred <- as.data.frame(cbind(Densities, Deg_Het, Diam, Clust, Assrt, Betw, Tmode))
colnames(DF_pred) <- c("density", "CVdeg", "Diam", "Clust", "DegAssort", "Betw", "Tmode")
write.csv(DF_pred, "Model_Predictions.csv")




