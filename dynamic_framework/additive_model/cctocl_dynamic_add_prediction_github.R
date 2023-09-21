############################################################################
# This script is for the analysis of large-mammals in North Anatolia.
# with static multispecies occupancy models

## This is the script to apply multi-season multi-species occupancy models.
# According to  Waddle et al 2010, Ecological Applications
# Capreolus capreolus as dominant and Canis lupus as sub-ordinant

# Prediction script for model outputs of dynamic multispecies modesl

# Additive models

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
  library(ggplot2)
  library(ggthemes)
  library(MCMCvis)
  library(rphylopic)
  library(bayesplot)
  library(cowplot)
  
}

load.libraries()

## 1.3. Loading data ----
#########################

setwd("path_to_data_files/")

# load data
load("path_to_data/cctocl_dynamic_multisp_data_github.RData")

# load model output
load("path_to_model_output/cctocl_dynamic_add_model_out.RData")

###########################################################################
# 2. Predictions
###########################################################################

# colonization (gamma) = season + popden + elevation + re(area)
# desertion (epsilon) = season + popden + elevation + re(area)

# Winter = 1
# Summer = 0

str(cctocl_dynamic_add_072023$sims.list)


popden_d <- data.frame(popden = seq(min(cov_multiseason_072023$popden), max(cov_multiseason_072023$popden), length.out = 10))
popden_d$popden_std <- scale(popden_d$popden)[,1]

elevation_d <- data.frame(elevation = seq(min(cov_multiseason_072023$elevation), max(cov_multiseason_072023$elevation), length.out = 10))
elevation_d$elevation_std <- scale(elevation_d$elevation)[,1]

## 2.1. Colonization ~ season + popden + elevation + re(area) ----
##################################################################

## 2.1.1. popden effect (elevation = 0 (mean)) ----
##########################################################

# season X prey p/a X popden.length.out X mcmc list
# we don't need to plot different areas, it was just for accounting in the model

col_cctocl_dynamic_add_popden <- array(NA, dim = c(2, 2, length(popden_d$popden_std),cctocl_dynamic_add_072023$mcmc.info$n.samples))

for (i in 1:length(popden_d$popden_std)){ # for popden
  for (a in 1:length(unique(cl_ms_dat$area))){ # for area
    
    # summer
    
    # prey absent
    
    col_cctocl_dynamic_add_popden[1,1,i,] <- plogis(cctocl_dynamic_add_072023$sims.list$alphagamma.cl[,1]
                                                    +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,1,1]*0 # season = summer
                                                    +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,2,1]*popden_d$popden_std[i] # popden
                                                    +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,3,1]*0 # elev = 0, mean
                                                    +cctocl_dynamic_add_072023$sims.list$randomgamma.cl[,1,a]) 
    
    # prey present
    
    col_cctocl_dynamic_add_popden[1,2,i,] <- plogis(cctocl_dynamic_add_072023$sims.list$alphagamma.cl[,2]
                                                    +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,1,2]*0 # season = summer
                                                    +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,2,2]*popden_d$popden_std[i] # popden
                                                    +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,3,2]*0 # elev = 0, mean
                                                    +cctocl_dynamic_add_072023$sims.list$randomgamma.cl[,2,a]) 
    
    # winter
    
    # prey absent
    
    col_cctocl_dynamic_add_popden[2,1,i,] <- plogis(cctocl_dynamic_add_072023$sims.list$alphagamma.cl[,1]
                                                    +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,1,1]*1 # season = winter
                                                    +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,2,1]*popden_d$popden_std[i] # popden
                                                    +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,3,1]*0 # elev = 0, mean
                                                    +cctocl_dynamic_add_072023$sims.list$randomgamma.cl[,1,a]) 
    
    # prey present
    
    col_cctocl_dynamic_add_popden[2,2,i,] <- plogis(cctocl_dynamic_add_072023$sims.list$alphagamma.cl[,2]
                                                    +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,1,2]*1 # season = winter
                                                    +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,2,2]*popden_d$popden_std[i] # popden
                                                    +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,3,2]*0 # elev = 0, mean
                                                    +cctocl_dynamic_add_072023$sims.list$randomgamma.cl[,2,a]) 
    
  }
}

str(col_cctocl_dynamic_add_popden)  

# then we take the mean of the mcmc list
pm.col_cctocl_dynamic_add_popden <- apply(col_cctocl_dynamic_add_popden, c(1,2,3), mean)
str(pm.col_cctocl_dynamic_add_popden)

# then calculate the credible intervals
CRI.col_cctocl_dynamic_add_popden <- apply(col_cctocl_dynamic_add_popden, c(1,2,3), function(x) quantile(x, c(0.025, 0.975)))
str(CRI.col_cctocl_dynamic_add_popden)


col_prob_cctocl_dynamic_add_popden <-expand.grid(season = c("Summer to Winter", "Winter to Summer"), roedeer = c("Roe Deer Absent","Roe Deer Present"), popden=popden_d$popden)
col_prob_cctocl_dynamic_add_popden$popden_std <- rep(popden_d$popden_std, each = 4)
col_prob_cctocl_dynamic_add_popden$pred <- c(pm.col_cctocl_dynamic_add_popden)
col_prob_cctocl_dynamic_add_popden$lower <- c(CRI.col_cctocl_dynamic_add_popden[1,,,])
col_prob_cctocl_dynamic_add_popden$upper <- c(CRI.col_cctocl_dynamic_add_popden[2,,,])

head(col_prob_cctocl_dynamic_add_popden)


## 2.1.2. elevation effect (popden = 0 (mean)) ----
##########################################################

# season X prey p/a X elevation.length.out X mcmc list
# we don't need to plot different areas, it was just for accounting in the model

col_cctocl_dynamic_add_elevation <- array(NA, dim = c(2, 2, length(elevation_d$elevation_std),cctocl_dynamic_add_072023$mcmc.info$n.samples))

for (i in 1:length(elevation_d$elevation_std)){ # for popden
  for (a in 1:length(unique(cl_ms_dat$area))){ # for area
    
    # summer
    
    # prey absent
    
    col_cctocl_dynamic_add_elevation[1,1,i,] <- plogis(cctocl_dynamic_add_072023$sims.list$alphagamma.cl[,1]
                                                       +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,1,1]*0 # season = summer
                                                       +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,2,1]*0 # popden = 0, mean
                                                       +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,3,1]*elevation_d$elevation_std[i] # elev
                                                       +cctocl_dynamic_add_072023$sims.list$randomgamma.cl[,1,a]) 
    
    # prey present
    
    col_cctocl_dynamic_add_elevation[1,2,i,] <- plogis(cctocl_dynamic_add_072023$sims.list$alphagamma.cl[,2]
                                                       +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,1,2]*0 # season = summer
                                                       +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,2,2]*0 # popden = 0, mean
                                                       +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,3,2]*elevation_d$elevation_std[i] # elev
                                                       +cctocl_dynamic_add_072023$sims.list$randomgamma.cl[,2,a]) 
    
    # winter
    
    # prey absent
    
    col_cctocl_dynamic_add_elevation[2,1,i,] <- plogis(cctocl_dynamic_add_072023$sims.list$alphagamma.cl[,1]
                                                       +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,1,1]*1 # season = winter
                                                       +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,2,1]*0 # popden = 0, mean
                                                       +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,3,1]*elevation_d$elevation_std[i] # elev
                                                       +cctocl_dynamic_add_072023$sims.list$randomgamma.cl[,1,a]) 
    
    # prey present
    
    col_cctocl_dynamic_add_elevation[2,2,i,] <- plogis(cctocl_dynamic_add_072023$sims.list$alphagamma.cl[,2]
                                                       +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,1,2]*1 # season = winter
                                                       +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,2,2]*0 # popden = 0, mean
                                                       +cctocl_dynamic_add_072023$sims.list$betagamma.cl[,3,2]*elevation_d$elevation_std[i] # elev
                                                       +cctocl_dynamic_add_072023$sims.list$randomgamma.cl[,2,a]) 
    
  }
}

str(col_cctocl_dynamic_add_elevation)  

# then we take the mean of the mcmc list
pm.col_cctocl_dynamic_add_elevation <- apply(col_cctocl_dynamic_add_elevation, c(1,2,3), mean)
str(pm.col_cctocl_dynamic_add_elevation)

# then calculate the credible intervals
CRI.col_cctocl_dynamic_add_elevation <- apply(col_cctocl_dynamic_add_elevation, c(1,2,3), function(x) quantile(x, c(0.025, 0.975)))
str(CRI.col_cctocl_dynamic_add_elevation)


col_prob_cctocl_dynamic_add_elevation <-expand.grid(season = c("Summer to Winter", "Winter to Summer"), roedeer = c("Roe Deer Absent","Roe Deer Present"), elevation=elevation_d$elevation)
col_prob_cctocl_dynamic_add_elevation$elevation_std <- rep(elevation_d$elevation_std, each = 4)
col_prob_cctocl_dynamic_add_elevation$pred <- c(pm.col_cctocl_dynamic_add_elevation)
col_prob_cctocl_dynamic_add_elevation$lower <- c(CRI.col_cctocl_dynamic_add_elevation[1,,,])
col_prob_cctocl_dynamic_add_elevation$upper <- c(CRI.col_cctocl_dynamic_add_elevation[2,,,])

head(col_prob_cctocl_dynamic_add_elevation)

## 2.1.3. Plot colonization ----
##########################################################

col_prob_cctocl_dynamic_add_popden_plot <- ggplot(col_prob_cctocl_dynamic_add_popden)+
  geom_ribbon(aes(x= popden, ymin = lower, ymax = upper, fill= season),alpha=0.5) +
  #scale_fill_brewer(type = "qual",palette=6,direction = 1)+
  geom_line(aes(popden, pred,group=season, colour=season), lwd=2, linetype=1)+
  facet_grid(season~roedeer) + 
  #scale_colour_brewer(type = "qual",palette=6,direction = 1)+
  ylim(c(0,1)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())+
  ylab(expression(paste("Colonization Probabilitiy", " ", (gamma)))) +
  xlab(expression("Rural human population density "(per~'4'~km^2))) +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        panel.spacing = unit(1, "lines"),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size=10))

col_prob_cctocl_dynamic_add_popden_plot


col_prob_cctocl_dynamic_add_elevation_plot <- ggplot(col_prob_cctocl_dynamic_add_elevation)+
  geom_ribbon(aes(x= elevation, ymin = lower, ymax = upper, fill= season),alpha=0.5) +
  #scale_fill_brewer(type = "qual",palette=6,direction = 1)+
  geom_line(aes(elevation, pred,group=season, colour=season), lwd=2, linetype=1)+
  facet_grid(season~roedeer) + 
  #scale_colour_brewer(type = "qual",palette=6,direction = 1)+
  ylim(c(0,1)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())+
  ylab(expression(paste("Colonization Probabilitiy", " ", (gamma)))) +
  xlab("Elevation (m)") +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        panel.spacing = unit(1, "lines"),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size=10))

col_prob_cctocl_dynamic_add_elevation_plot


## 2.2. Desertion ~ season + popden + elevation + re(area) ----
##################################################################

## 2.2.1. popden effect (elevation = 0 (mean)) ----
##########################################################

# season X prey p/a X popden.length.out X mcmc list
# we don't need to plot different areas, it was just for accounting in the model

ext_cctocl_dynamic_add_popden <- array(NA, dim = c(2, 2, length(popden_d$popden_std),cctocl_dynamic_add_072023$mcmc.info$n.samples))

for (i in 1:length(popden_d$popden_std)){ # for popden
  for (a in 1:length(unique(cl_ms_dat$area))){ # for area
    
    # summer
    
    # prey absent
    
    ext_cctocl_dynamic_add_popden[1,1,i,] <- plogis(cctocl_dynamic_add_072023$sims.list$alphaeps.cl[,1]
                                                    +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,1,1]*0 # season = summer
                                                    +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,2,1]*popden_d$popden_std[i] # popden
                                                    +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,3,1]*0 # elev = 0, mean
                                                    +cctocl_dynamic_add_072023$sims.list$randomeps.cl[,1,a]) 
    
    # prey present
    
    ext_cctocl_dynamic_add_popden[1,2,i,] <- plogis(cctocl_dynamic_add_072023$sims.list$alphaeps.cl[,2]
                                                    +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,1,2]*0 # season = summer
                                                    +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,2,2]*popden_d$popden_std[i] # popden
                                                    +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,3,2]*0 # elev = 0, mean
                                                    +cctocl_dynamic_add_072023$sims.list$randomeps.cl[,2,a]) 
    
    # winter
    
    # prey absent
    
    ext_cctocl_dynamic_add_popden[2,1,i,] <- plogis(cctocl_dynamic_add_072023$sims.list$alphaeps.cl[,1]
                                                    +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,1,1]*1 # season = winter
                                                    +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,2,1]*popden_d$popden_std[i] # popden
                                                    +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,3,1]*0 # elev = 0, mean
                                                    +cctocl_dynamic_add_072023$sims.list$randomeps.cl[,1,a]) 
    
    # prey present
    
    ext_cctocl_dynamic_add_popden[2,2,i,] <- plogis(cctocl_dynamic_add_072023$sims.list$alphaeps.cl[,2]
                                                    +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,1,2]*1 # season = winter
                                                    +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,2,2]*popden_d$popden_std[i] # popden
                                                    +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,3,2]*0 # elev = 0, mean
                                                    +cctocl_dynamic_add_072023$sims.list$randomeps.cl[,2,a]) 
    
  }
}

str(ext_cctocl_dynamic_add_popden)  

# then we take the mean of the mcmc list
pm.ext_cctocl_dynamic_add_popden <- apply(ext_cctocl_dynamic_add_popden, c(1,2,3), mean)
str(pm.ext_cctocl_dynamic_add_popden)

# then calculate the credible intervals
CRI.ext_cctocl_dynamic_add_popden <- apply(ext_cctocl_dynamic_add_popden, c(1,2,3), function(x) quantile(x, c(0.025, 0.975)))
str(CRI.ext_cctocl_dynamic_add_popden)


ext_prob_cctocl_dynamic_add_popden <-expand.grid(season = c("Summer to Winter", "Winter to Summer"), roedeer = c("Roe Deer Absent","Roe Deer Present"), popden=popden_d$popden)
ext_prob_cctocl_dynamic_add_popden$popden_std <- rep(popden_d$popden_std, each = 4)
ext_prob_cctocl_dynamic_add_popden$pred <- c(pm.ext_cctocl_dynamic_add_popden)
ext_prob_cctocl_dynamic_add_popden$lower <- c(CRI.ext_cctocl_dynamic_add_popden[1,,,])
ext_prob_cctocl_dynamic_add_popden$upper <- c(CRI.ext_cctocl_dynamic_add_popden[2,,,])

head(ext_prob_cctocl_dynamic_add_popden)


## 2.2.2. elevation effect (popden = 0 (mean)) ----
##########################################################

# season X prey p/a X elevation.length.out X mcmc list
# we don't need to plot different areas, it was just for accounting in the model

ext_cctocl_dynamic_add_elevation <- array(NA, dim = c(2, 2, length(elevation_d$elevation_std),cctocl_dynamic_add_072023$mcmc.info$n.samples))

for (i in 1:length(elevation_d$elevation_std)){ # for popden
  for (a in 1:length(unique(cl_ms_dat$area))){ # for area
    
    # summer
    
    # prey absent
    
    ext_cctocl_dynamic_add_elevation[1,1,i,] <- plogis(cctocl_dynamic_add_072023$sims.list$alphaeps.cl[,1]
                                                       +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,1,1]*0 # season = summer
                                                       +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,2,1]*0 # popden = 0, mean
                                                       +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,3,1]*elevation_d$elevation_std[i] # elev
                                                       +cctocl_dynamic_add_072023$sims.list$randomeps.cl[,1,a]) 
    
    # prey present
    
    ext_cctocl_dynamic_add_elevation[1,2,i,] <- plogis(cctocl_dynamic_add_072023$sims.list$alphaeps.cl[,2]
                                                       +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,1,2]*0 # season = summer
                                                       +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,2,2]*0 # popden = 0, mean
                                                       +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,3,2]*elevation_d$elevation_std[i] # elev
                                                       +cctocl_dynamic_add_072023$sims.list$randomeps.cl[,2,a]) 
    
    # winter
    
    # prey absent
    
    ext_cctocl_dynamic_add_elevation[2,1,i,] <- plogis(cctocl_dynamic_add_072023$sims.list$alphaeps.cl[,1]
                                                       +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,1,1]*1 # season = winter
                                                       +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,2,1]*0 # popden = 0, mean
                                                       +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,3,1]*elevation_d$elevation_std[i] # elev
                                                       +cctocl_dynamic_add_072023$sims.list$randomeps.cl[,1,a]) 
    
    # prey present
    
    ext_cctocl_dynamic_add_elevation[2,2,i,] <- plogis(cctocl_dynamic_add_072023$sims.list$alphaeps.cl[,2]
                                                       +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,1,2]*1 # season = winter
                                                       +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,2,2]*0 # popden = 0, mean
                                                       +cctocl_dynamic_add_072023$sims.list$betaeps.cl[,3,2]*elevation_d$elevation_std[i] # elev
                                                       +cctocl_dynamic_add_072023$sims.list$randomeps.cl[,2,a]) 
    
  }
}

str(ext_cctocl_dynamic_add_elevation)  

# then we take the mean of the mcmc list
pm.ext_cctocl_dynamic_add_elevation <- apply(ext_cctocl_dynamic_add_elevation, c(1,2,3), mean)
str(pm.ext_cctocl_dynamic_add_elevation)

# then calculate the credible intervals
CRI.ext_cctocl_dynamic_add_elevation <- apply(ext_cctocl_dynamic_add_elevation, c(1,2,3), function(x) quantile(x, c(0.025, 0.975)))
str(CRI.ext_cctocl_dynamic_add_elevation)


ext_prob_cctocl_dynamic_add_elevation <-expand.grid(season = c("Summer to Winter", "Winter to Summer"), roedeer = c("Roe Deer Absent","Roe Deer Present"), elevation=elevation_d$elevation)
ext_prob_cctocl_dynamic_add_elevation$elevation_std <- rep(elevation_d$elevation_std, each = 4)
ext_prob_cctocl_dynamic_add_elevation$pred <- c(pm.ext_cctocl_dynamic_add_elevation)
ext_prob_cctocl_dynamic_add_elevation$lower <- c(CRI.ext_cctocl_dynamic_add_elevation[1,,,])
ext_prob_cctocl_dynamic_add_elevation$upper <- c(CRI.ext_cctocl_dynamic_add_elevation[2,,,])

head(ext_prob_cctocl_dynamic_add_elevation)

## 2.2.3. Plot desertion ----
##########################################################

ext_prob_cctocl_dynamic_add_popden_plot <- ggplot(ext_prob_cctocl_dynamic_add_popden)+
  geom_ribbon(aes(x= popden, ymin = lower, ymax = upper, fill= season),alpha=0.5) +
  #scale_fill_brewer(type = "qual",palette=6,direction = 1)+
  geom_line(aes(popden, pred,group=season, colour=season), lwd=2, linetype=1)+
  facet_grid(season~roedeer) + 
  #scale_colour_brewer(type = "qual",palette=6,direction = 1)+
  ylim(c(0,1)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())+
  ylab(expression(paste("Desertion Probabilitiy", " ", (epsilon)))) +
  xlab(expression("Rural human population density "(per~'4'~km^2))) +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        panel.spacing = unit(1, "lines"),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size=10))

ext_prob_cctocl_dynamic_add_popden_plot


ext_prob_cctocl_dynamic_add_elevation_plot <- ggplot(ext_prob_cctocl_dynamic_add_elevation)+
  geom_ribbon(aes(x= elevation, ymin = lower, ymax = upper, fill= season),alpha=0.5) +
  #scale_fill_brewer(type = "qual",palette=6,direction = 1)+
  geom_line(aes(elevation, pred,group=season, colour=season), lwd=2, linetype=1)+
  facet_grid(season~roedeer) + 
  #scale_colour_brewer(type = "qual",palette=6,direction = 1)+
  ylim(c(0,1)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())+
  ylab(expression(paste("Desertion Probabilitiy", " ", (epsilon)))) +
  xlab("Elevation (m)") +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        panel.spacing = unit(1, "lines"),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size=10))

ext_prob_cctocl_dynamic_add_elevation_plot

## 2.3. Site occupancy states ----
##########################################################
plot(cctocl_dynamic_add_072023$mean$nWD, col="red", type='b', ylim=c(0,100))
lines(cctocl_dynamic_add_072023$mean$nWd, col="blue", type='b', ylim=c(0,100))
lines(cctocl_dynamic_add_072023$mean$nwD, col="green", type='b', ylim=c(0,100))
lines(cctocl_dynamic_add_072023$mean$nwd, col="black", type='b', ylim=c(0,100))

range_fin_occu <- range(c(cctocl_dynamic_add_072023$mean$nWD,
                          cctocl_dynamic_add_072023$mean$nWd,
                          cctocl_dynamic_add_072023$mean$nwD,
                          cctocl_dynamic_add_072023$mean$nwd))
range_fin_occu

fin_occu_total_cctocl_dynamic_add <- data.frame(season_nr=rep(1:22,4), season=rep(c("W","S")))

fin_occu_total_cctocl_dynamic_add$pred[1:22] <- cctocl_dynamic_add_072023$mean$nWD
fin_occu_total_cctocl_dynamic_add$pred[23:44] <- cctocl_dynamic_add_072023$mean$nWd
fin_occu_total_cctocl_dynamic_add$pred[45:66] <- cctocl_dynamic_add_072023$mean$nwD
fin_occu_total_cctocl_dynamic_add$pred[67:88] <- cctocl_dynamic_add_072023$mean$nwd

fin_occu_total_cctocl_dynamic_add$state <- rep(c("RoeDeer(1) Wolf(1)", "RoeDeer(0) Wolf(1)", "RoeDeer(1) Wolf(0)","RoeDeer(0) Wolf(0)"), 
                                               each=22, length.out=88)
fin_occu_total_cctocl_dynamic_add$state_wrapped <- stringr::str_wrap(fin_occu_total_cctocl_dynamic_add$state, width = 1)

fin_occu_total_cctocl_dynamic_add$state2 <- rep(c("Prey(1) Wolf(1)", "Prey(0) Wolf(1)", "Prey(1) Wolf(0)","Prey(0) Wolf(0)"), 
                                                each=22, length.out=88)
fin_occu_total_cctocl_dynamic_add$state2_wrapped <- stringr::str_wrap(fin_occu_total_cctocl_dynamic_add$state2, width = 1)

fin_occu_total_cctocl_dynamic_add$state3 <- rep(c("Co-occupied", "Predator only", "Prey only","Unoccupied"), 
                                                each=22, length.out=88)
fin_occu_total_cctocl_dynamic_add$state3_wrapped <- stringr::str_wrap(fin_occu_total_cctocl_dynamic_add$state3, width = 1)

fin_occu_total_cctocl_dynamic_add$lower[1:22] <- cctocl_dynamic_add_072023$q2.5$nWD
fin_occu_total_cctocl_dynamic_add$lower[23:44] <- cctocl_dynamic_add_072023$q2.5$nWd
fin_occu_total_cctocl_dynamic_add$lower[45:66] <- cctocl_dynamic_add_072023$q2.5$nwD
fin_occu_total_cctocl_dynamic_add$lower[67:88] <- cctocl_dynamic_add_072023$q2.5$nwd

fin_occu_total_cctocl_dynamic_add$upper[1:22] <- cctocl_dynamic_add_072023$q97.5$nWD
fin_occu_total_cctocl_dynamic_add$upper[23:44] <- cctocl_dynamic_add_072023$q97.5$nWd
fin_occu_total_cctocl_dynamic_add$upper[45:66] <- cctocl_dynamic_add_072023$q97.5$nwD
fin_occu_total_cctocl_dynamic_add$upper[67:88] <- cctocl_dynamic_add_072023$q97.5$nwd

fin_occu_total_cctocl_dynamic_add_plot <- ggplot(fin_occu_total_cctocl_dynamic_add) + 
  geom_line(aes(season_nr, pred, color=state3, group=state3), lwd=1)+ 
  geom_ribbon(aes(season_nr, ymin=lower, ymax=upper, fill=state3), alpha=0.2) +
  scale_fill_manual(values = c("#7570b3","#0571b0","#e66101","gray40")) + 
  scale_color_manual(values = c("#7570b3","#0571b0","#e66101","gray40")) + 
  #facet_wrap(.~state, ncol=2) +
  scale_y_continuous(breaks=seq(0,130, by=10)) + 
  scale_x_continuous(breaks = 1:23,labels=c(rep(c("W","S"), 11),"W"))+
  theme_bw() + 
  theme(panel.grid.minor = element_blank())+
  ylab(expression("Number of sites")) +
  xlab(expression("Seasons")) +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        axis.title.y = element_text(size=20, face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20,  margin = margin(t = 20, r = 0, b = 0, l = 0)),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position = "none",
        legend.direction = "horizontal")+
  guides(fill=guide_legend(nrow=2,byrow=TRUE, title.position = "top"))


fin_occu_total_cctocl_dynamic_add_plot

