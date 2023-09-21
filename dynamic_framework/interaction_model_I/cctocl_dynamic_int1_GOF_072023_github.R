############################################################################
# This script is for the re-analysis 
# MULTISEASON, MULTI SPECIES OCCUPANCY MODELS
# of large mammals in North Anatolia.

# GOODNESS OF FIT TEST FOR

# Capreolus capreolus on Canus lupus: multispecies dynamic occupancy models.

# Interaction model 1

# Data preparation for this script is in multispecies_data_prep_072023.R

# Date: September 2023
# Author: Dilsad Dagtekin
###########################################################################

###########################################################################
# 1. House keeping and loading libraries and data
###########################################################################

## 1.1. House keeping ----
##########################

rm(list = ls())

## 1.2. Loading libraries ----
##############################

load.libraries <- function(){
  library(jagsUI)
  library(MCMCvis)
}

load.libraries()

## 1.3. Loading data ----
#########################

setwd("/path_to_files/")
load("path_to_data/cctocl_dynamic_multisp_data_github.RData") # load prepared jags data

###########################################################################
# 2. Analysis of the model
###########################################################################

## 2.1. Specify model in BUGS language
######################################

cat(file = "cctocl_dynamic_int1_072023_gof.txt"," 
model {

# Model: cctocl_dynamic_int1_072023_gof
#
# initial occupancy (psi1) = habitat + popden + habitat:popden + re(area)
# colonization (gamma) = season + popden + season:popden + re(area)
# desertion (epsilon) = season + popden + season:popden + re(area)
# detection (p) = season + habitat + season:habitat + re(area)
# -----------------------------------------------------------------------------------------
#
# Parameters:
#
# psi1: initial occupancy probability
# gamma: colonization probability
# epsilon: desertion probability
# p: detection probability
# -----------------------------------------------------------------------------------------
#
# Others:
#
# narea: number of study areas, 10
# nsite: number of camera-trap stations, 171
# nprimary: total number of primary sampling periods, seasons, 22
# nobs: number of observations
# y: detection history
# z: true latent state variable, z=1 means occupied, z=0 not occupied.
# n.occ: number of occupied sites
# n.prop: proportion of occupied sites to all total sites
#
# area: study areas
# season: seasons, 2-level categorical, winter(1) and summer (0)
# popden: rural human population density, continuous, standardized
# elevation: meter, continuous, standardized
# habitat: habitat type, 4-level categorical:
## Broad Leaved Forest (BL) = 1
## Coniferous Forest (CF) = 2 
## Mixed Forest (MF) = 3
## Human Land-Use Areas (O) = 4
# 
# -----------------------------------------------------------------------------------------

## Model for Capreolus capreolus - Roe Deer (dominant species)
##############################################################

#### PRIORS ####

# psi1 ~ habitat + popden + habitat:popden + re.psi1(area)

for (h in 1:4) {
  alphapsi.cc[h] ~ dlogis(0, 1)                                    # intercept for each habitat type(4)
  
  # betap values should be according to habitat types (AHM II Chp15 pp.104-105)
  # we have 4 habitat types --> 4 betapsi values for 1 betapsi category
  # we have 1 betapsi category: popden effect
  
  betapsi.cc[h] ~ dnorm(0, 1)                                      # prior for popden coefficient w/ habitat types
}

# gamma ~ season + popden + elevation + re.gamma(area)

alphagamma.cc ~ dlogis(0, 1)                                       # intercept

for(i in 1:3){                                                  # prior for coefficients 
  betagamma.cc[i] ~ dnorm(0, 1)
}

# eps ~ season + popden + elevation + re.eps(area)

alphaeps.cc ~ dlogis(0, 1)                                         # intercept

for(i in 1:3){                                                  # prior for coefficients 
  betaeps.cc[i] ~ dnorm(0, 1)
}

# p ~ season*habitat

for (h in 1:4) {
  # alphap.cc[h] ~ dlogis(0, 1)                                      # intercept for each habitat type(4)
  
  alphap.cc[h] ~ dlogis(-2, 1)
    
  # betap values should be according to habitat types (AHM II Chp15 pp.104-105)
  # we have 4 habitat types --> 4 betap values for 1 betap category
  # we have 1 betap category: betap[p] season effect
  
  betap.cc[h] ~ dnorm(0, 1)                                        # betap season effect
}

# Random effect priors
  
for (t in 1:(narea)){
randompsi.cc[t] ~ dnorm(0, taupsi.cc)
randomgamma.cc[t] ~ dnorm(0, taugamma.cc)
randomeps.cc[t] ~ dnorm(0, taueps.cc)
randomp.cc[t] ~ dnorm(0, taup.cc)

    
}

sigmapsi.cc ~ dnorm(0, 1)T(0,)                    # Priors for standard deviations
sigmagamma.cc ~ dnorm(0, 1)T(0,)
sigmaeps.cc ~ dnorm(0, 1)T(0,)
sigmap.cc ~ dnorm(0, 1)T(0,)

taupsi.cc <- pow(sigmapsi.cc, -2)
taugamma.cc <- pow(sigmagamma.cc, -2)
taueps.cc <- pow(sigmaeps.cc, -2)
taup.cc <- pow(sigmap.cc, -2)

sigma2psi.cc <- pow(sigmapsi.cc, 2)                  # Temporal variances
sigma2gamma.cc <- pow(sigmagamma.cc, 2)
sigma2eps.cc <- pow(sigmaeps.cc, 2)
sigma2p.cc <- pow(sigmap.cc, 2)

#### MODELS ####

# Ecological submodel: intial occupancy as derived parameter

for (i in 1:nsite){
  z.cc[i,1] ~ dbern(psi1.cc[i])
  logit(psi1.cc[i]) <- alphapsi.cc[habitat[i]] 
  + betapsi.cc[habitat[i]]*popden[i]
  + randompsi.cc[area[i]]
}

# State transitions: colonization and extinction
# for col and ext

for (i in 1:nsite){
  for (s in 2:nprimary){
    logit(gamma.cc[i,s-1]) <- alphagamma.cc 
    + betagamma.cc[1]*season[i,s-1] + betagamma.cc[2]*popden[i] + betagamma.cc[3]*season[i,s-1]*popden[i]
    + randomgamma.cc[area[i]]
    
    logit(eps.cc[i,s-1]) <- alphaeps.cc 
    + betaeps.cc[1]*season[i,s-1] + betaeps.cc[2]*popden[i] + betaeps.cc[3]*season[i,s-1]*popden[i]
    + randomeps.cc[area[i]]
  }
}

# for z: true occupancy

for(i in 1:nsite){
  for (t in 2:nprimary){
    
    z.cc[i,t] ~ dbern((z.cc[i,t-1]*(1-eps.cc[i,t-1])) + ((1-z.cc[i,t-1])*gamma.cc[i,t-1]))
  }
}

# Observation model 

for (i in 1:nobs.cc){
  logit(p.cc[i]) <- alphap.cc[habitat.p.cc[i]]
  + betap.cc[habitat.p.cc[i]]*season.p.cc[site.cc[i]]
  + randomp.cc[area.cc[i]]
  
   #y.cc[i] ~ dbern(z.cc[site.cc[i],pocc.cc[i]]*p.cc[i])  # y: observed occupancy
  
   pstar.cc[i] <- 1-(1-p.cc[i])^ndays.cc[i] 
   y.cc[i] ~ dbern(z.cc[site.cc[i],pocc.cc[i]]*pstar.cc[i])  # y: observed occupancy
}

## Model for Canis lupus - Gray wolf (sub-ordinate species)
###########################################################

#### PRIORS ####

for (k in 1:2){                                               # Presence/absence of dominant species
  
  # psi1 ~ habitat + popden + habitat:popden + re.psi1(area)
  
  for (h in 1:4) {
    alphapsi.cl[k,h] ~ dlogis(0, 1)                                    # intercept for each habitat type(4)
    
    # betap values should be according to habitat types (AHM II Chp15 pp.104-105)
    # we have 4 habitat types --> 4 betapsi values for 1 betapsi category
    # we have 1 betapsi category: popden effect
    
    betapsi.cl[k,h] ~ dnorm(0, 1)                                      # prior for popden coefficient w/ habitat types
  }
  
  # gamma ~ season + popden + elevation + re.gamma(area)
  
  alphagamma.cl[k] ~ dlogis(0, 1)                                       # intercept
  
  for(i in 1:3){                                                  # prior for coefficients 
    betagamma.cl[i,k] ~ dnorm(0, 1)
  }
  
  # eps ~ season + popden + elevation + re.eps(area)
  
  alphaeps.cl[k] ~ dlogis(0, 1)                                         # intercept
  
  for(i in 1:3){                                                  # prior for coefficients 
    betaeps.cl[i,k] ~ dnorm(0, 1)
  }
  
  # p ~ season*habitat
  
  for (h in 1:4) {
    # alphap.cl[k,h] ~ dlogis(0, 1)                                      # intercept for each habitat type(4)
    
    alphap.cl[k,h] ~ dlogis(-2, 1)
    
    # betap values should be according to habitat types (AHM II Chp15 pp.104-105)
    # we have 4 habitat types --> 4 betap values for 1 betap category
    # we have 1 betap category: betap[p] season effect
    
    betap.cl[k,h] ~ dnorm(0, 1)                                        # betap season effect
  }
  
  # Random effect priors
  
  for (t in 1:(narea)){
    randompsi.cl[k,t] ~ dnorm(0, taupsi.cl[k])
    randomgamma.cl[k,t] ~ dnorm(0, taugamma.cl[k])
    randomeps.cl[k,t] ~ dnorm(0, taueps.cl[k])
    randomp.cl[k,t] ~ dnorm(0, taup.cl[k])
    
    
  }
  
  sigmapsi.cl[k] ~ dnorm(0, 1)T(0,)                      # Priors for standard deviations
  sigmagamma.cl[k] ~ dnorm(0, 1)T(0,)   
  sigmaeps.cl[k] ~ dnorm(0, 1)T(0,)   
  sigmap.cl[k] ~ dnorm(0, 1)T(0,)   
  
  taupsi.cl[k] <- pow(sigmapsi.cl[k], -2)
  taugamma.cl[k] <- pow(sigmagamma.cl[k], -2)
  taueps.cl[k] <- pow(sigmaeps.cl[k], -2)
  taup.cl[k] <- pow(sigmap.cl[k], -2)
  
  sigma2psi.cl[k] <- pow(sigmapsi.cl[k], 2)                  # Temporal variances
  sigma2gamma.cl[k] <- pow(sigmagamma.cl[k], 2)
  sigma2eps.cl[k] <- pow(sigmaeps.cl[k], 2)
  sigma2p.cl[k] <- pow(sigmap.cl[k], 2)
}

#### LIKELIHOOD ####

# Ecological submodel: intial occupancy as derived parameter

for (i in 1:nsite){
  z.cl[i,1] ~ dbern(psi1.cl[i])
  logit(psi1.cl[i]) <- alphapsi.cl[z.cc[site.cc[i],pocc.cc[i]]+1, habitat[i]] 
  + betapsi.cl[z.cc[site.cc[i],pocc.cc[i]]+1, habitat[i]]*popden[i]
  + randompsi.cl[z.cc[site.cc[i],pocc.cc[i]]+1, area[i]]
}

# State transitions: colonization and extinction
# for col and ext

for (i in 1:nsite){
  for (s in 2:nprimary){
    logit(gamma.cl[i,s-1]) <- alphagamma.cl[z.cc[site.cc[i],pocc.cc[i]]+1] 
    + betagamma.cl[1,z.cc[site.cc[i],pocc.cc[i]]+1]*season[i,s-1] + betagamma.cl[2,z.cc[site.cc[i],pocc.cc[i]]+1]*popden[i] 
    + betagamma.cl[3,z.cc[site.cc[i],pocc.cc[i]]+1]*season[i,s-1]*popden[i] + randomgamma.cl[z.cc[site.cc[i],pocc.cc[i]]+1,area[i]]
    
    logit(eps.cl[i,s-1]) <- alphaeps.cl[z.cc[site.cc[i],pocc.cc[i]]+1] 
    + betaeps.cl[1,z.cc[site.cc[i],pocc.cc[i]]+1]*season[i,s-1] + betaeps.cl[2,z.cc[site.cc[i],pocc.cc[i]]+1]*popden[i] 
    + betaeps.cl[3,z.cc[site.cc[i],pocc.cc[i]]+1]*season[i,s-1]*popden[i] + randomeps.cl[z.cc[site.cc[i],pocc.cc[i]]+1, area[i]]
  }
}

# for z: true occupancy

for(i in 1:nsite){
  for (t in 2:nprimary){
    
    z.cl[i,t] ~ dbern((z.cl[i,t-1]*(1-eps.cl[i,t-1])) + ((1-z.cl[i,t-1])*gamma.cl[i,t-1]))
  }
}

# Observation model 

for (i in 1:nobs.cl){
  logit(p.cl[i]) <- alphap.cl[z.cc[site.cc[i],pocc.cc[i]]+1,habitat.p.cl[i]]
  + betap.cl[z.cc[site.cc[i],pocc.cc[i]]+1,habitat.p.cl[i]]*season.p.cl[site.cl[i]]
  + randomp.cl[z.cc[site.cc[i],pocc.cc[i]]+1,area.cl[i]]
  
  # y.cl[i] ~ dbern(z.cl[site.cl[i],pocc.cl[i]]*p.cl[i])  # y: observed occupancy
  
  pstar.cl[i] <- 1-(1-p.cl[i])^ndays.cl[i] 
  y.cl[i] ~ dbern(z.cl[site.cl[i],pocc.cl[i]]*pstar.cl[i])  # y: observed occupancy
}


#### GOODNESS-OF-FIT GoF ####
# Based on posterior predictive distribution
# -------------------------------------------------------------------------

# derived occupancy state
for(i in 1:nsite){
 
 psi.cl[i, 1] <- psi1.cl[i] # Population occupancy as derived quantity
 for (t in 2:nprimary){
   
   psi.cl[i, t] <- (psi.cl[i,t-1]*(1-eps.cl[i,t-1])) + ((1-psi.cl[i,t-1])*gamma.cl[i,t-1])
 
 }
}

# Draw a replicate data set under the fitted model

for (i in 1:nobs.cl){
 
   yrep.cl[i] ~ dbern(z.cl[site.cl[i],pocc.cl[i]]*pstar.cl[i])
}

}")


# Initial values
inits<-function(){list(z.cl=matrix(1,171,22), z.cc=matrix(1,171,22))}

# Parameters monitored
params <- c("yrep.cl", "psi.cl", "gamma.cl","eps.cl","z.cl","pstar.cl")

# MCMC settings 
na <-1; ni <- 10; nt <- 1 ; nb <- 5 ; nc <- 3 
# na <- 50000 ; ni <- 500000; nt <- 20 ; nb <- 50000 ; nc <- 3

# parallel:::setDefaultClusterOptions(setup_strategy = "sequential") # after R update this is needed for parallel run
# Call JAGS from R and summarize posteriors

start <- Sys.time()

cctocl_dynamic_int1_072023_gof <- jags(jags.data, inits, params, "cctocl_dynamic_int1_072023_gof.txt", n.chains = nc, n.adapt = na, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = T)

end <- Sys.time()

end - start

print(cctocl_dynamic_int1_072023_gof)

str(cctocl_dynamic_int1_072023_gof$sims.list)


# GOF ---------------------------------------------------------------------

###########################################################################
# FOR ALL MCMC SAMPLES ----------------------------------------------------
###########################################################################

dat_gof <- clupus_dat

rownames(dat_gof) <- NULL
str(dat_gof)
dat_gof$site <- as.character(dat_gof$site)
dat_gof$socc <- as.character(dat_gof$socc)
dat_gof$pocc <- as.character(dat_gof$pocc)

dat_gof$soccNr <- ifelse(dat_gof$month == 12, 1,
                         ifelse(dat_gof$month == 1, 2,
                                ifelse(dat_gof$month == 2, 3,
                                       ifelse(dat_gof$month == 3, 4,
                                              ifelse(dat_gof$month == 5, 1,
                                                     ifelse(dat_gof$month == 6, 2,
                                                            ifelse(dat_gof$month == 7, 3,
                                                                   ifelse(dat_gof$month == 8, 4,
                                                                          ifelse(dat_gof$month == 9, 5, 
                                                                                 ifelse(dat_gof$month == 10, 6, NA))))))))))


# create all the output from the model

# taking out the replicated observation data from the model
yrep.cl <- list(NA)
dat_rep <- list(NA)


for (i in 1:cctocl_dynamic_int1_072023_gof$mcmc.info$n.samples){
  yrep.cl[[i]] <- cctocl_dynamic_int1_072023_gof$sims.list$yrep.cl[i,]
  
  dat_rep[[i]] <- dat_gof[,c(1:9, 25)]
  dat_rep[[i]]$y <- yrep.cl[[i]]
  
}
str(dat_rep[[1]])
str(dat_rep[[2]])

# taking out the estimated detection prob from the model
p_model <- list(NA)
dat_p <- list(NA)

for (i in 1:cctocl_dynamic_int1_072023_gof$mcmc.info$n.samples){
  p_model[[i]] <- cctocl_dynamic_int1_072023_gof$sims.list$pstar[i,]
  
  dat_p[[i]] <- dat_gof[,c(1:9, 25)]
  dat_p[[i]]$y <- yrep.cl[[i]]
  dat_p[[i]]$pstar <- round(p_model[[i]],2)
  
}

# we need to convert detection histories into 3d array again
# months start with 12/2007

# winter 12/1/2/3
# summer 5/6/7/8/9/10

sites_arr <- 1:171
months_arr <- 1:6
seasons_arr <- 1:11

## observed data: dat_gof$y, only one array
y3dobs <- array(NA, dim = c(171, 6, 22), dimnames = list(1:171, 1:6, 1:22))

for (i in 1:171){
  for (j in 1:6){
    for (t in 1:22){  
      
      y3dobs[i,j,t] <- ifelse(sites_arr[i] %in% dat_gof[which(dat_gof$y==1 & dat_gof$pocc==t & dat_gof$soccNr==j),]$site, 1,
                              ifelse(sites_arr[i] %in% dat_gof[which(dat_gof$y==0 & dat_gof$pocc==t & dat_gof$soccNr==j),]$site, 0, NA))
    }
  }
}


# replicated data: dat_rep, for all MCMC samples
y3drep <- vector(mode = 'list', length = cctocl_dynamic_int1_072023_gof$mcmc.info$n.samples)

for (k in 1:cctocl_dynamic_int1_072023_gof$mcmc.info$n.samples){
  y3drep[[k]] <- array(NaN, dim = c(171, 6, 22), dimnames = list(1:171, 1:6, 1:22))
}

for (k in 1:cctocl_dynamic_int1_072023_gof$mcmc.info$n.samples){
  for (i in 1:171){
    for (j in 1:6){
      for (t in 1:22){  
        
        y3drep[[k]][i,j,t] <- ifelse(sites_arr[i] %in% dat_rep[[k]][which(dat_rep[[k]]$y==1 & dat_rep[[k]]$pocc==t & dat_rep[[k]]$soccNr==j),]$site, 1,
                                     ifelse(sites_arr[i] %in% dat_rep[[k]][which(dat_rep[[k]]$y==0 & dat_rep[[k]]$pocc==t & dat_rep[[k]]$soccNr==j),]$site, 0, NA))
      }
    }
  }
}

# we also need to convert estimated pstar (detection prob) into 3d array

p_est <- vector(mode = 'list', length = cctocl_dynamic_int1_072023_gof$mcmc.info$n.samples)

for (k in 1:cctocl_dynamic_int1_072023_gof$mcmc.info$n.samples){
  p_est[[k]] <- array(NaN, dim = c(171, 6, 22), dimnames = list(1:171, 1:6, 1:22))
}

for (k in 1:cctocl_dynamic_int1_072023_gof$mcmc.info$n.samples){
  for (i in 1:171){
    for (j in 1:6){
      for (t in 1:22){  
        
        p_est[[k]][i,j,t] <- ifelse(length(dat_p[[k]][which(dat_p[[k]]$site==i & dat_p[[k]]$pocc ==t & dat_p[[k]]$soccNr==j),]$pstar)!= 0, 
                                    dat_p[[k]][which(dat_p[[k]]$site==i & dat_p[[k]]$pocc ==t & dat_p[[k]]$soccNr  ==j),]$pstar, NA)
      }
    }
  }
}

# Computations for the GoF of the open part of the model
# (based on number of state transitions)
# ----------------------------------------------------------

nsite = dim(y3drep[[1]])[1] # number of sites
nprimary = dim(y3drep[[1]])[3] # number of seasons (primar occasions)
e = 0.0001 # to avoid dividing by zero

psi.cl <- list(NA)
gamma.cl <- list(NA)
eps.cl <- list(NA)

for (k in 1:cctocl_dynamic_int1_072023_gof$mcmc.info$n.samples){
  psi.cl[[k]] <- matrix(NaN, nrow=nsite, ncol=nprimary)
  gamma.cl[[k]] <- matrix(NaN, nrow=nsite, ncol=nprimary-1)
  eps.cl[[k]] <- matrix(NaN, nrow=nsite, ncol=nprimary-1)
  
}


for (i in 1:cctocl_dynamic_int1_072023_gof$mcmc.info$n.samples){
  psi.cl[[i]] <- cctocl_dynamic_int1_072023_gof$sims.list$psi.cl[i,,]
  gamma.cl[[i]] <- cctocl_dynamic_int1_072023_gof$sims.list$gamma.cl[i,,]
  eps.cl[[i]] <- cctocl_dynamic_int1_072023_gof$sims.list$eps.cl[i,,]
  
}

# Create needed empty matrices

# observed data
zobs <- matrix(NaN, nrow = nsite, ncol=nprimary) # empty matrix to save latent state from observed data

des <- matrix(NA, nrow = nsite, ncol=nprimary-1) 
col <- matrix(NA, nrow = nsite, ncol=nprimary-1)
nondes <- matrix(NA, nrow = nsite, ncol=nprimary-1)
noncol <- matrix(NA, nrow = nsite, ncol=nprimary-1)

tm <- array(NA, dim = c(2,2,nprimary-1))

# replicated data
zobsrep <- list(NA)

desrep <- list(NA)
colrep <- list(NA)
nondesrep <- list(NA)
noncolrep <- list(NA)

tmrep <- list(NA)

noncol.exp <- list(NA)
col.exp <- list(NA)
des.exp <- list(NA)
nondes.exp <- list(NA)

Etm <- list(NA)

x2Open <- list(NA)
x2repOpen <- list(NA)

Chi2Open <- NA
Chi2repOpen <- NA
Chi2ratioOpen <- NA

for (k in 1:cctocl_dynamic_int1_072023_gof$mcmc.info$n.samples){
  
  zobsrep[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary)  # empty matrix to save latent state from replicated data
  
  desrep[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary-1)
  colrep[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary-1)
  nondesrep[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary-1)
  noncolrep[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary-1)
  
  tmrep[[k]] <- array(NaN, dim = c(2,2,nprimary-1))
  
  noncol.exp[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary-1)
  col.exp[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary-1)
  des.exp[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary-1)
  nondes.exp[[k]] <- matrix(NaN, nrow = nsite, ncol=nprimary-1)
  
  Etm[[k]] <- array(NaN, dim = c(2,2,nprimary-1))
  
  x2Open[[k]] <- array(NaN, dim = c(2,2,nprimary-1))
  x2repOpen[[k]] <- array(NaN, dim = c(2,2,nprimary-1))
  
}

# Calculate the Chi-square

# Compute observed z matrix for observed and replicated data
for (k in 1:cctocl_dynamic_int1_072023_gof$mcmc.info$n.samples){
  
  for (i in 1:nsite){
    for (t in 1:nprimary){
      zobs[i,t] <- max(y3dobs[i,,t], na.rm=T)  # For observed data
      zobsrep[[k]][i,t] <- max(y3drep[[k]][i,,t], na.rm = T) # For replicated data
    }
    
    # Identify desertions, non-desertions, colonization and non-colonizations
    for (t in 2:nprimary){
      # ... for observed data
      des[i,(t-1)] <- ifelse(zobs[i,t]==0 & zobs[i,t-1]==1,1,0)
      nondes[i,(t-1)] <- ifelse(zobs[i,t]==1 & zobs[i,t-1]==1,1,0)
      col[i,(t-1)] <- ifelse(zobs[i,t]==1 & zobs[i,t-1]==0,1,0)
      noncol[i,(t-1)] <- ifelse(zobs[i,t]==0 & zobs[i,t-1]==0,1,0)
      # ... for replicated data
      desrep[[k]][i,(t-1)] <- ifelse(zobsrep[[k]][i,t]==0 & zobsrep[[k]][i,t-1]==1,1,0)
      nondesrep[[k]][i,(t-1)] <- ifelse(zobsrep[[k]][i,t]==1 & zobsrep[[k]][i,t-1]==1,1,0)
      colrep[[k]][i,(t-1)] <- ifelse(zobsrep[[k]][i,t]==1 & zobsrep[[k]][i,t-1]==0,1,0)
      noncolrep[[k]][i,(t-1)] <- ifelse(zobsrep[[k]][i,t]==0 & zobsrep[[k]][i,t-1]==0,1,0)
    }
  }
  
  # Tally up number of transitions and put into a matrix for each year
  for(t in 1:(nprimary-1)){
    # ... for observed data
    tm[1,1,t] <- sum(noncol[,t]) # transition mat for obs. data
    tm[1,2,t] <- sum(col[,t])
    tm[2,1,t] <- sum(des[,t])
    tm[2,2,t] <- sum(nondes[,t])
    # ... for replicated data
    tmrep[[k]][1,1,t] <- sum(noncolrep[[k]][,t]) # transition mat for rep. data
    tmrep[[k]][1,2,t] <- sum(colrep[[k]][,t])
    tmrep[[k]][2,1,t] <- sum(desrep[[k]][,t])
    tmrep[[k]][2,2,t] <- sum(nondesrep[[k]][,t])
  }
  
  # Compute expected numbers of transitions under the model
  # Probability of each individual transition
  for(i in 1:nsite){
    for(t in 1:(nprimary-1)){
      noncol.exp[[k]][i,t] <- (1-psi.cl[[k]][i,t]) * (1-gamma.cl[[k]][i,t])
      col.exp[[k]][i,t] <- (1-psi.cl[[k]][i,t]) * gamma.cl[[k]][i,t]
      des.exp[[k]][i,t] <- psi.cl[[k]][i,t] * eps.cl[[k]][i,t]
      nondes.exp[[k]][i,t] <- psi.cl[[k]][i,t] * (1-eps.cl[[k]][i,t])
    }
  }
  
  
  
  # Sum up over sites to obtain the expected number of those transitions
  for(t in 1:(nprimary-1)){
    Etm[[k]][1,1,t] <- sum(noncol.exp[[k]][,t])
    Etm[[k]][1,2,t] <- sum(col.exp[[k]][,t])
    Etm[[k]][2,1,t] <- sum(des.exp[[k]][,t])
    Etm[[k]][2,2,t] <- sum(nondes.exp[[k]][,t])
  }
  
  
  # Compute Chi-square discrepancy
  for(t in 1:(nprimary-1)){
    # ... for observed data
    x2Open[[k]][1,1,t] <- ((tm[1,1,t] - Etm[[k]][1,1,t]) ^ 2) / (Etm[[k]][1,1,t]+e)
    x2Open[[k]][1,2,t] <- ((tm[1,2,t] - Etm[[k]][1,2,t]) ^ 2) / (Etm[[k]][1,2,t]+e)
    x2Open[[k]][2,1,t] <- ((tm[2,1,t] - Etm[[k]][2,1,t]) ^ 2) / (Etm[[k]][2,1,t]+e)
    x2Open[[k]][2,2,t] <- ((tm[2,2,t] - Etm[[k]][2,2,t]) ^ 2) / (Etm[[k]][2,2,t]+e)
    # ... for replicated data
    x2repOpen[[k]][1,1,t] <- ((tmrep[[k]][1,1,t]-Etm[[k]][1,1,t]) ^ 2)/(Etm[[k]][1,1,t]+e)
    x2repOpen[[k]][1,2,t] <- ((tmrep[[k]][1,2,t]-Etm[[k]][1,2,t]) ^ 2)/(Etm[[k]][1,2,t]+e)
    x2repOpen[[k]][2,1,t] <- ((tmrep[[k]][2,1,t]-Etm[[k]][2,1,t]) ^ 2)/(Etm[[k]][2,1,t]+e)
    x2repOpen[[k]][2,2,t] <- ((tmrep[[k]][2,2,t]-Etm[[k]][2,2,t]) ^ 2)/(Etm[[k]][2,2,t]+e)
  }
  
  # Add up overall test statistic and compute fit stat ratio (open part)
  Chi2Open[[k]] <- sum(x2Open[[k]][,,])       # Chisq. statistic for observed data
  Chi2repOpen[[k]] <- sum(x2repOpen[[k]][,,]) # Chisq. statistic for replicated data
  Chi2ratioOpen[[k]] <- Chi2Open[[k]] / Chi2repOpen[[k]]
}

# plot Open part
pl <- range(c(Chi2Open, Chi2repOpen))
plot(Chi2Open, Chi2repOpen,
     xlab = "Chi2 observed data", ylab = "Chi2 expected data",
     main = "Open part of model ", xlim = pl, ylim = pl, frame.plot = FALSE)
abline(0, 1, lwd = 2)
text(2685, 2755, paste('Bpv = ', round(mean(Chi2repOpen >
                                              Chi2Open), 2)), cex = 2)