############################################################################
# This script is for the analysis of large-mammals in North Anatolia.
# with static multispecies occupancy models

## This is the script to apply single-season multi-species occupancy models.
# According to  Waddle et al 2010, Ecological Applications
# Capreolus capreolus as dominant and Canis lupus as sub-ordinant

# static models:
# occupancy (psi) = habitat + popden + habitat:popden + re(area)
# detection (p) = habitat + re(area) 
# p model was season*habitat before, 
# but since I am analyzing winter and summer sites separately, we can just do habitat and re(area)

# additionally:
# we need to add camera trap effort into the detection probability
# for (i in 1:nobs){
#   logit(p[i]) <- alphap + betap
#   Pstar[i]<-1-(1-p[i]) ^ Ndays[i]
#   y[i] ~ dbern(z[site[i]]*Pstar[i]) # y: observed occupancy
# } 

# We analyzed two seasons seperately to account for season effect.
# This is the code for SUMMER season

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
load("path_to_data/cctocl_static_summer_data_github.RData") # load prepared jags data

###########################################################################
# 2. Analysis of the model
###########################################################################

## 2.1. Specify model in BUGS language
######################################

cat(file = "cctocl_static_summer_072023.txt"," 
model {

# Model: cctocl_static_summer_072023
#
# psi ~ habitat + popden + habitat:popden + re.psi1(area)
# p ~ habitat + re.p(area)
# -----------------------------------------------------------------------------------------
#
# Parameters:
#
# psi: occupancy probability
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

# psi ~ habitat + popden + habitat:popden + re.psi(area)

for (h in 1:4) {
  alphapsi.cc[h] ~ dlogis(0, 1)                                    # intercept for each habitat type(4)
  
  # betap values should be according to habitat types (AHM II Chp15 pp.104-105)
  # we have 4 habitat types --> 4 betapsi values for 1 betapsi category
  # we have 1 betapsi category: popden effect
  
  betapsi.cc[h] ~ dnorm(0, 1)                                      # prior for popden coefficient w/ habitat types
}


# p ~ habitat + re.p(area)

for (h in 1:4) {
  alphap.cc[h] ~ dlogis(-2, 1)                                      # intercept for each habitat type(4)
  
}

# Random effect priors
  
for (t in 1:(narea)){
randompsi.cc[t] ~ dnorm(0, taupsi.cc)
randomp.cc[t] ~ dnorm(0, taup.cc)

    
}

sigmapsi.cc ~ dnorm(0, 1)T(0,)                     # Priors for standard deviations
sigmap.cc ~ dnorm(0, 1)T(0,) 

taupsi.cc <- pow(sigmapsi.cc, -2)
taup.cc <- pow(sigmap.cc, -2)

sigma2psi.cc <- pow(sigmapsi.cc, 2)                  # Temporal variances
sigma2p.cc<- pow(sigmap.cc, 2)

#### MODELS ####

# Ecological submodel: intial occupancy as derived parameter

for (i in 1:nsite){

  logit(psi.cc[i]) <- alphapsi.cc[habitat.cc[i]] + betapsi.cc[habitat.cc[i]]*popden.cc[i] + randompsi.cc[area.cc[i]]
  
  z.cc[i] ~ dbern(psi.cc[i])                                # True occupancy state
}

# Observation model 

for (i in 1:nobs.cc){
  logit(p.cc[i]) <- alphap.cc[habitat.p.cc[i]] + randomp.cc[area.p.cc[i]]
  
    # y.cc[i] ~ dbern(z.cc[site.cc[i]]*p.cc[i])  # y: observed occupancy
    
    pstar.cc[i] <- 1-(1-p.cc[i])^ndays.cc[i] 
    y.cc[i] ~ dbern(z.cc[site.cc[i]]*pstar.cc[i])  # y: observed occupancy
}

## Model for Canis lupus - Gray wolf (sub-ordinate species)
###########################################################

#### PRIORS ####

for (k in 1:2){                                               # Presence/absence of dominant species
  
  # psi ~ habitat + popden + habitat:popden + re.psi(area)
  
  for (h in 1:4) {
    alphapsi.cl[k,h] ~ dlogis(0, 1)                                    # intercept for each habitat type(4)
    
    # betap values should be according to habitat types (AHM II Chp15 pp.104-105)
    # we have 4 habitat types --> 4 betapsi values for 1 betapsi category
    # we have 1 betapsi category: popden effect
    
    betapsi.cl[k,h] ~ dnorm(0, 1)                                      # prior for popden coefficient w/ habitat types
  }
  
  
  # p ~ habitat + re.p(area)
  
  for (h in 1:4) {
    alphap.cl[k,h] ~ dlogis(-2, 1)                                      # intercept for each habitat type(4)
    
  }
  
  # Random effect priors
  
  for (t in 1:(narea)){
    randompsi.cl[k,t] ~ dnorm(0, taupsi.cl[k])
    randomp.cl[k,t] ~ dnorm(0, taup.cl[k])
    
    
  }
  
  sigmapsi.cl[k] ~ dnorm(0, 1)T(0,)                     # Priors for standard deviations
  sigmap.cl[k] ~ dnorm(0, 1)T(0,) 
  
  taupsi.cl[k] <- pow(sigmapsi.cl[k], -2)
  taup.cl[k] <- pow(sigmap.cl[k], -2)
  
  sigma2psi.cl[k] <- pow(sigmapsi.cl[k], 2)                  # Temporal variances
  sigma2p.cl[k] <- pow(sigmap.cl[k], 2)
  
}
#### MODELS ####

# Ecological submodel: intial occupancy as derived parameter

for (i in 1:nsite){
  
  logit(psi.cl[i]) <- alphapsi.cl[z.cc[site.cc[i]]+1,habitat.cl[i]] + betapsi.cl[z.cc[site.cc[i]]+1,habitat.cl[i]]*popden.cl[i] + randompsi.cl[z.cc[site.cc[i]]+1,area.cl[i]]
  
  z.cl[i] ~ dbern(psi.cl[i])                                # True occupancy state
}

# Observation model 

for (i in 1:nobs.cl){
  logit(p.cl[i]) <- alphap.cl[z.cc[site.cc[i]]+1,habitat.p.cl[i]] + randomp.cl[z.cc[site.cc[i]]+1,area.p.cl[i]]
  
  # y.cl[i] ~ dbern(z.cl[site.cl[i]]*p.cl[i])  # y: observed occupancy
  
  pstar.cl[i] <- 1-(1-p.cl[i])^ndays.cl[i] 
  y.cl[i] ~ dbern(z.cl[site.cl[i]]*pstar.cl[i])  # y: observed occupancy
}

## Differences in Clupus parameters between Ccapreolus presence and absence##

diff_betapsi_1 <- betapsi.cl[2,1] - betapsi.cl[1,1] # hab BL
diff_betapsi_2 <- betapsi.cl[2,2] - betapsi.cl[1,2] # hab CF
diff_betapsi_3 <- betapsi.cl[2,3] - betapsi.cl[1,3] # hab MF
diff_betapsi_4 <- betapsi.cl[2,4] - betapsi.cl[1,4] # hab O (HLU)

diff_alphapsi_1 <- alphapsi.cl[2,1] - alphapsi.cl[1,1] # hab BL
diff_alphapsi_2 <- alphapsi.cl[2,2] - alphapsi.cl[1,2] # hab CF
diff_alphapsi_3 <- alphapsi.cl[2,3] - alphapsi.cl[1,3] # hab MF
diff_alphapsi_4 <- alphapsi.cl[2,4] - alphapsi.cl[1,4] # hab O (HLU)

diff_alphap_1 <- alphap.cl[2,1] - alphap.cl[1,1] # hab BL
diff_alphap_2 <- alphap.cl[2,2] - alphap.cl[1,2] # hab CF
diff_alphap_3 <- alphap.cl[2,3] - alphap.cl[1,3] # hab MF
diff_alphap_4 <- alphap.cl[2,4] - alphap.cl[1,4] # hab O (HLU)

## Total number of sites for each state ##

for (i in 1:nsite){
  DW_S[i] <- z.cc[i] * z.cl[i]
  Dw_S[i] <- z.cc[i] * (1-z.cl[i])
  dW_S[i] <- (1-z.cc[i]) * z.cl[i]
  dw_S[i] <- (1-z.cc[i]) * (1-z.cl[i])
}

nDW_S <- sum(DW_S[])        # Number of sites with both species present

nDw_S <- sum(Dw_S[])        # Number of sites with only roe deer present

ndW_S <- sum(dW_S[])        # Number of sites with only gray wolf present

ndw_S <- sum(dw_S[])        # Number of sites with both species absent

## Species Interaction Factor - SIF ##

for(i in 1:nsite){

  # derived occupancy prob for wolf when roe deer absent
  
  dpsi.cl[i] <- ilogit(alphapsi.cl[1,habitat.cl[i]] + betapsi.cl[1,habitat.cl[i]] + randompsi.cl[1,area.cl[i]])

  # derived occupancy prob for wolf when roe deer present
  
  dpsi.cl.D[i] <- ilogit(alphapsi.cl[2,habitat.cl[i]] + betapsi.cl[2,habitat.cl[i]] + randompsi.cl[2,area.cl[i]])
  
  # derived occupancy prob for roe deer
  
  dpsi.cc[i] <- ilogit(alphapsi.cc[habitat.cc[i]] + betapsi.cc[habitat.cc[i]] + randompsi.cc[area.cc[i]])
  
  SIF[i] <- (dpsi.cc[i]*dpsi.cl.D[i])/(dpsi.cc[i]*((dpsi.cc[i]*dpsi.cl.D[i])+((1-dpsi.cc[i])*dpsi.cl[i])))
}

}")


## Model specifics ----
######################################

# Initial values
inits<-function(){list(z.cc=rep(1,171),z.cl=rep(1,171))}

# Parameters monitored
params <- c("alphapsi.cc", "betapsi.cc",
            "alphap.cc",
            "randompsi.cc","randomp.cc",
            "alphapsi.cl", "betapsi.cl",
            "alphap.cl", 
            "randompsi.cl","randomp.cl",
            "diff_betapsi_1","diff_betapsi_2","diff_betapsi_3","diff_betapsi_4",
            "diff_alphapsi_1", "diff_alphapsi_2", "diff_alphapsi_3", "diff_alphapsi_4", 
            "diff_alphap_1", "diff_alphap_2", "diff_alphap_3", "diff_alphap_4", 
            "nDW_S", "nDw_S", "ndW_S", "ndw_S",
            "SIF")

# MCMC settings 

# to try the model
# na <-1; ni <- 10; nt <- 1 ; nb <- 5 ; nc <- 3

# to run the model
na <- 50000 ; ni <- 400000; nt <- 20 ; nb <- 50000 ; nc <- 3

# Call JAGS from R and summarize posteriors

start <- Sys.time()

parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
cctocl_static_summer_072023 <- jags(jags.data, inits, params, "cctocl_static_summer_072023.txt", n.chains = nc, n.adapt = na, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = T)

end <- Sys.time()

end - start

print(cctocl_static_summer_072023)
