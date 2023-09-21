############################################################################
# This script is for the analysis of large-mammals in North Anatolia.
# with dynamic multispecies occupancy models

## This is the script to apply multi-season multi-species occupancy models.
# According to  Waddle et al 2010, Ecological Applications
# Capreolus capreolus as dominant and Canis lupus as sub-ordinant

# interaction 2 models: season*elevation
# initial occupancy (psi1) = habitat + popden + habitat:popden + re(area)
# colonization (gamma) = season + elevation + season:elevation + re(area)
# desertion (epsilon) = season + elevation + season:elevation + re(area)
# detection (p) = season + habitat + season:habitat + re(area)

# Date: July 2023
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

cat(file = "cctocl_dynamic_int2_072023.txt"," 
model {

# Model: cctocl_dynamic_int2_072023
#
# initial occupancy (psi1) = habitat + popden + habitat:popden + re(area)
# colonization (gamma) = season + elevation + season:elevation + re(area)
# desertion (epsilon) = season + elevation + season:elevation + re(area)
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
    + betagamma.cc[1]*season[i,s-1] + betagamma.cc[2]*elevation[i] + betagamma.cc[3]*season[i,s-1]*elevation[i]
    + randomgamma.cc[area[i]]
    
    logit(eps.cc[i,s-1]) <- alphaeps.cc 
    + betaeps.cc[1]*season[i,s-1] + betaeps.cc[2]*elevation[i] + betaeps.cc[3]*season[i,s-1]*elevation[i]
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

#### DERIVED PARAMETERS ####
# Compute population and sample occupancy

n.occ.cc[1] <- sum(z.cc[1:nsite,1])  # Number of occupied sites in sample
n.prop.cc[1] <- n.occ.cc[1]/nsite

for (t in 2:nprimary){
  n.occ.cc[t] <- sum(z.cc[1:nsite,t])
  n.prop.cc[t] <- n.occ.cc[t]/nsite
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
    + betagamma.cl[1,z.cc[site.cc[i],pocc.cc[i]]+1]*season[i,s-1] + betagamma.cl[2,z.cc[site.cc[i],pocc.cc[i]]+1]*elevation[i] 
    + betagamma.cl[3,z.cc[site.cc[i],pocc.cc[i]]  +1]*season[i,s-1]*elevation[i] + randomgamma.cl[z.cc[site.cc[i],pocc.cc[i]]+1,area[i]]
    
    logit(eps.cl[i,s-1]) <- alphaeps.cl[z.cc[site.cc[i],pocc.cc[i]]+1] 
    + betaeps.cl[1,z.cc[site.cc[i],pocc.cc[i]]+1]*season[i,s-1] + betaeps.cl[2,z.cc[site.cc[i],pocc.cc[i]]+1]*elevation[i] 
    + betaeps.cl[3,z.cc[site.cc[i],pocc.cc[i]]+1]*season[i,s-1]*elevation[i] + randomeps.cl[z.cc[site.cc[i],pocc.cc[i]]+1, area[i]]
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

#### DERIVED PARAMETERS ####
# Compute population and sample occupancy

n.occ.cl[1] <- sum(z.cl[1:nsite,1])  # Number of occupied sites in sample
n.prop.cl[1] <- n.occ.cl[1]/nsite

for (t in 2:nprimary){
  n.occ.cl[t] <- sum(z.cl[1:nsite,t])
  n.prop.cl[t] <- n.occ.cl[t]/nsite
}

# Derived parameters for species interactions

# Number of sites in each of 4 states
for(t in 1:nprimary){                   
  for (i in 1:nsite){
    WD[i,t] <- z.cl[i,t]*z.cc[i,t]
    Wd[i,t] <- z.cl[i,t]*(1-z.cc[i,t])
    wD[i,t] <- (1-z.cl[i,t])*z.cc[i,t]
    wd[i,t] <- (1-z.cl[i,t])*(1-z.cc[i,t])
  }
  nWD[t] <- sum(WD[,t])                                           # Number of sites with both species
  nWd[t] <- sum(Wd[,t])                                           # Number of sites with only clupus
  nwD[t] <- sum(wD[,t])                                           # Number of sites with only ccapreolus
  nwd[t] <- sum(wd[,t])                                           # Number of sites with neither species
}
  
}")


# Initial values
inits<-function(){list(z.cl=matrix(1,171,22), z.cc=matrix(1,171,22))}

# Parameters monitored
params <- c("alphapsi.cl", "betapsi.cl", 
            "alphagamma.cl", "betagamma.cl",
            "alphaeps.cl","betaeps.cl", 
            "alphap.cl", "betap.cl",
            "randompsi.cl", "randomgamma.cl", "randomeps.cl", "randomp.cl",
            "n.occ.cl", "n.prop.cl",
            "alphapsi.cc", "betapsi.cc", 
            "alphagamma.cc", "betagamma.cc",
            "alphaeps.cc","betaeps.cc", 
            "alphap.cc", "betap.cc",
            "randompsi.cc", "randomgamma.cc", "randomeps.cc", "randomp.cc",
            "n.occ.cc", "n.prop.cc",
            "nWD", "nWd", "nwD","nwd")

# to try the model
# na <-1; ni <- 10; nt <- 1 ; nb <- 5 ; nc <- 3

# to run the model
na <- 50000 ; ni <- 400000; nt <- 20 ; nb <- 50000 ; nc <- 3

# Call JAGS from R and summarize posteriors

start <- Sys.time()

parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

cctocl_dynamic_int2_072023 <- jags(jags.data, inits, params, "cctocl_dynamic_int2_072023.txt", n.chains = nc, n.adapt = na, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = T)

end <- Sys.time()

end - start

print(cctocl_dynamic_int2_072023)
