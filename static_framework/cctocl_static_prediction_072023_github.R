############################################################################
# This script is for the analysis of large-mammals in North Anatolia.
# with static multispecies occupancy models

## This is the script to apply single-season multi-species occupancy models.
# According to  Waddle et al 2010, Ecological Applications
# Capreolus capreolus as dominant and Canis lupus as sub-ordinant

# Prediction script for model outputs of static multispecies modesl
# SUMMER and WINTER together

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
  library(eulerr)
  
}

load.libraries()

## 1.3. Loading data ----
#########################

setwd("path_to_data_files/")

# load data
load("path_to_data/cctocl_static_summer_data_github.RData")
load("path_to_data/cctocl_static_winter_data_github.RData")

# load model output
load("path_to_model_output/cctocl_static_summer_model_out.RData")
load("path_to_model_output/cctocl_static_winter_model_out.RData")

###########################################################################
# 2. Predictions
###########################################################################

# Broad Leaved Forest (BL) = 1
# Coniferous Forest (CF) = 2 
# Mixed Forest (MF) = 3
# Human Land-Use Areas (O) = 4

# Winter = 1
# Summer = 0

str(cctocl_static_summer_072023$sims.list)
str(cctocl_static_winter_072023$sims.list)

## 2.1. Detection ~ habitat + re.p(area) ----
##########################################################
# Broad Leaved Forest (BL) = 1
# Coniferous Forest (CF) = 2 
# Mixed Forest (MF) = 3
# Human Land-Use Areas (O) = 4

## 2.1.1. Summer ----
##########################################################

# prey p/a X habitat X mcmc list
# we don't need to plot different areas, it was just for accounting in the model

det_cctocl_static_summer <- array(NA, dim = c(2, 4,cctocl_static_summer_072023$mcmc.info$n.samples))

for (i in 1:length(unique(cc_ss_dat_summer$area))){ # for area
  
  # prey absent
  
  det_cctocl_static_summer[1,1,] <- plogis(cctocl_static_summer_072023$sims.list$alphap.cl[,1,1] # BL
                                           +cctocl_static_summer_072023$sims.list$randomp.cl[,1,i]) 
  
  det_cctocl_static_summer[1,2,] <- plogis(cctocl_static_summer_072023$sims.list$alphap.cl[,1,2] # CF
                                           +cctocl_static_summer_072023$sims.list$randomp.cl[,1,i]) 
  
  det_cctocl_static_summer[1,3,] <- plogis(cctocl_static_summer_072023$sims.list$alphap.cl[,1,3] # MF
                                           +cctocl_static_summer_072023$sims.list$randomp.cl[,1,i]) 
  
  det_cctocl_static_summer[1,4,] <- plogis(cctocl_static_summer_072023$sims.list$alphap.cl[,1,4] # O
                                           +cctocl_static_summer_072023$sims.list$randomp.cl[,1,i]) 
  
  # prey present
  
  det_cctocl_static_summer[2,1,] <- plogis(cctocl_static_summer_072023$sims.list$alphap.cl[,2,1] # BL
                                           +cctocl_static_summer_072023$sims.list$randomp.cl[,2,i]) 
  
  det_cctocl_static_summer[2,2,] <- plogis(cctocl_static_summer_072023$sims.list$alphap.cl[,2,2] # CF
                                           +cctocl_static_summer_072023$sims.list$randomp.cl[,2,i]) 
  
  det_cctocl_static_summer[2,3,] <- plogis(cctocl_static_summer_072023$sims.list$alphap.cl[,2,3] # MF
                                           +cctocl_static_summer_072023$sims.list$randomp.cl[,2,i]) 
  
  det_cctocl_static_summer[2,4,] <- plogis(cctocl_static_summer_072023$sims.list$alphap.cl[,2,4] # O
                                           +cctocl_static_summer_072023$sims.list$randomp.cl[,2,i]) 
}

str(det_cctocl_static_summer)  

# then we take the mean of the mcmc list
pm.det_cctocl_static_summer <- apply(det_cctocl_static_summer, c(1,2), mean)
str(pm.det_cctocl_static_summer)

# then calculate the credible intervals
CRI.det_cctocl_static_summer <- apply(det_cctocl_static_summer, c(1,2), function(x) quantile(x, c(0.025, 0.975)))
str(CRI.det_cctocl_static_summer)

det_prob_cctocl_static_summer <-expand.grid(roedeer = c("Roe Deer Absent","Roe Deer Present"), habitat=c("Broad Leaved Forest", "Coniferous Forest", "Mixed Forest", "Human Land-Use Areas"))
det_prob_cctocl_static_summer$season <- "Summer"
det_prob_cctocl_static_summer$pred <- c(pm.det_cctocl_static_summer)
det_prob_cctocl_static_summer$lower <- c(CRI.det_cctocl_static_summer[1,,])
det_prob_cctocl_static_summer$upper <- c(CRI.det_cctocl_static_summer[2,,])
det_prob_cctocl_static_summer$species <- "Roe Deer"

head(det_prob_cctocl_static_summer)


## 2.1.2. Winter ----
##########################################################

# prey p/a X habitat X mcmc list
# we don't need to plot different areas, it was just for accounting in the model

det_cctocl_static_winter <- array(NA, dim = c(2, 4,cctocl_static_winter_072023$mcmc.info$n.samples))

for (i in 1:length(unique(cc_ss_dat_winter$area))){ # for area
  
  # prey absent
  
  det_cctocl_static_winter[1,1,] <- plogis(cctocl_static_winter_072023$sims.list$alphap.cl[,1,1] # BL
                                           +cctocl_static_winter_072023$sims.list$randomp.cl[,1,i]) 
  
  det_cctocl_static_winter[1,2,] <- plogis(cctocl_static_winter_072023$sims.list$alphap.cl[,1,2] # CF
                                           +cctocl_static_winter_072023$sims.list$randomp.cl[,1,i]) 
  
  det_cctocl_static_winter[1,3,] <- plogis(cctocl_static_winter_072023$sims.list$alphap.cl[,1,3] # MF
                                           +cctocl_static_winter_072023$sims.list$randomp.cl[,1,i]) 
  
  det_cctocl_static_winter[1,4,] <- plogis(cctocl_static_winter_072023$sims.list$alphap.cl[,1,4] # O
                                           +cctocl_static_winter_072023$sims.list$randomp.cl[,1,i]) 
  
  # prey present
  
  det_cctocl_static_winter[2,1,] <- plogis(cctocl_static_winter_072023$sims.list$alphap.cl[,2,1] # BL
                                           +cctocl_static_winter_072023$sims.list$randomp.cl[,2,i]) 
  
  det_cctocl_static_winter[2,2,] <- plogis(cctocl_static_winter_072023$sims.list$alphap.cl[,2,2] # CF
                                           +cctocl_static_winter_072023$sims.list$randomp.cl[,2,i]) 
  
  det_cctocl_static_winter[2,3,] <- plogis(cctocl_static_winter_072023$sims.list$alphap.cl[,2,3] # MF
                                           +cctocl_static_winter_072023$sims.list$randomp.cl[,2,i]) 
  
  det_cctocl_static_winter[2,4,] <- plogis(cctocl_static_winter_072023$sims.list$alphap.cl[,2,4] # O
                                           +cctocl_static_winter_072023$sims.list$randomp.cl[,2,i]) 
}

str(det_cctocl_static_winter)  

# then we take the mean of the mcmc list
pm.det_cctocl_static_winter <- apply(det_cctocl_static_winter, c(1,2), mean)
str(pm.det_cctocl_static_winter)

# then calculate the credible intervals
CRI.det_cctocl_static_winter <- apply(det_cctocl_static_winter, c(1,2), function(x) quantile(x, c(0.025, 0.975)))
str(CRI.det_cctocl_static_winter)

det_prob_cctocl_static_winter <-expand.grid(roedeer = c("Roe Deer Absent","Roe Deer Present"), habitat=c("Broad Leaved Forest", "Coniferous Forest", "Mixed Forest", "Human Land-Use Areas"))
det_prob_cctocl_static_winter$season <- "Winter"
det_prob_cctocl_static_winter$pred <- c(pm.det_cctocl_static_winter)
det_prob_cctocl_static_winter$lower <- c(CRI.det_cctocl_static_winter[1,,])
det_prob_cctocl_static_winter$upper <- c(CRI.det_cctocl_static_winter[2,,])
det_prob_cctocl_static_winter$species <- "Roe Deer"

head(det_prob_cctocl_static_winter)

## 2.1.3. Merge and plot detection ----
##########################################################

det_prob_cctocl_static <- rbind(det_prob_cctocl_static_summer, det_prob_cctocl_static_winter)

det_prob_cctocl_static$habitat_wrapped <- stringr::str_wrap(det_prob_cctocl_static$habitat, width = 3)

det_prob_cctocl_static_plot <- ggplot(det_prob_cctocl_static)+
  geom_errorbar(aes(habitat_wrapped ,ymin = lower, ymax = upper, col=season, group=habitat_wrapped), size = 1, width = 0.5, position = position_dodge(width = 0.8))  +
  scale_fill_continuous(type = "viridis") +
  geom_point(aes(habitat_wrapped, pred, col=season, group=habitat_wrapped), size=5, position = position_dodge(width = 0.8))  +
  facet_grid(season ~ roedeer) + 
  ylim(c(0,1)) + 
  theme_bw() + 
  ylab(expression(paste("Detection Probabilitiy", " ", (p)))) +
  xlab('Habitat Type') +
  theme(axis.text.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=16, face = "bold"),
        panel.spacing = unit(2, "lines"),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))

det_prob_cctocl_static_plot

## 2.2. Occupancy ~ habitat + popden + habitat:popden ----
##########################################################

popden_d <- data.frame(popden = seq(min(cl_ss_dat_summer$popden), max(cl_ss_dat_summer$popden), length.out = 10))
popden_d$popden_std <- scale(popden_d$popden)[,1]

# Broad Leaved Forest (BL) = 1
# Coniferous Forest (CF) = 2 
# Mixed Forest (MF) = 3
# Human Land-Use Areas (O) = 4

## 2.2.1. Summer ----
##########################################################

# prey p/a X habitat X length.out.popden X mcmc list

occ_cctocl_static_summer <- array(NA, dim = c(2,4,length(popden_d$popden_std),cctocl_static_summer_072023$mcmc.info$n.samples))

for (i in 1:length(popden_d$popden_std)){ # for popden
  for(a in 1:length(unique(cl_ss_dat_summer$area))){ # for area
    
    # prey absent
    
    occ_cctocl_static_summer[1,1,i,] <- plogis(cctocl_static_summer_072023$sims.list$alphapsi.cl[,1,1] # BL
                                               + cctocl_static_summer_072023$sims.list$betapsi.cl[,1,1]*popden_d$popden_std[i]
                                               + cctocl_static_summer_072023$sims.list$randompsi.cl[,1,a]) 
    
    occ_cctocl_static_summer[1,2,i,] <- plogis(cctocl_static_summer_072023$sims.list$alphapsi.cl[,1,2] # CF
                                               + cctocl_static_summer_072023$sims.list$betapsi.cl[,1,2]*popden_d$popden_std[i]
                                               + cctocl_static_summer_072023$sims.list$randompsi.cl[,1,a]) 
    
    occ_cctocl_static_summer[1,3,i,] <- plogis(cctocl_static_summer_072023$sims.list$alphapsi.cl[,1,3] # MF
                                               + cctocl_static_summer_072023$sims.list$betapsi.cl[,1,3]*popden_d$popden_std[i]
                                               + cctocl_static_summer_072023$sims.list$randompsi.cl[,1,a]) 
    
    occ_cctocl_static_summer[1,4,i,] <- plogis(cctocl_static_summer_072023$sims.list$alphapsi.cl[,1,4] # O
                                               + cctocl_static_summer_072023$sims.list$betapsi.cl[,1,4]*popden_d$popden_std[i]
                                               + cctocl_static_summer_072023$sims.list$randompsi.cl[,1,a]) 
    
    # prey present
    
    occ_cctocl_static_summer[2,1,i,] <- plogis(cctocl_static_summer_072023$sims.list$alphapsi.cl[,2,1] # BL
                                               + cctocl_static_summer_072023$sims.list$betapsi.cl[,2,1]*popden_d$popden_std[i]
                                               + cctocl_static_summer_072023$sims.list$randompsi.cl[,2,a]) 
    
    occ_cctocl_static_summer[2,2,i,] <- plogis(cctocl_static_summer_072023$sims.list$alphapsi.cl[,2,2] # CF
                                               + cctocl_static_summer_072023$sims.list$betapsi.cl[,2,2]*popden_d$popden_std[i]
                                               + cctocl_static_summer_072023$sims.list$randompsi.cl[,2,a]) 
    
    occ_cctocl_static_summer[2,3,i,] <- plogis(cctocl_static_summer_072023$sims.list$alphapsi.cl[,2,3] # MF
                                               + cctocl_static_summer_072023$sims.list$betapsi.cl[,2,3]*popden_d$popden_std[i]
                                               + cctocl_static_summer_072023$sims.list$randompsi.cl[,2,a]) 
    
    occ_cctocl_static_summer[2,4,i,] <- plogis(cctocl_static_summer_072023$sims.list$alphapsi.cl[,2,4] # O
                                               + cctocl_static_summer_072023$sims.list$betapsi.cl[,2,4]*popden_d$popden_std[i]
                                               + cctocl_static_summer_072023$sims.list$randompsi.cl[,2,a]) 
  }
}

str(occ_cctocl_static_summer)  

# then we take the mean of the mcmc list
pm.occ_cctocl_static_summer <- apply(occ_cctocl_static_summer, c(1,2,3), mean)
str(pm.occ_cctocl_static_summer)

# then calculate the credible intervals
CRI.occ_cctocl_static_summer <- apply(occ_cctocl_static_summer, c(1,2,3), function(x) quantile(x, c(0.025, 0.975)))
str(CRI.occ_cctocl_static_summer)

occ_prob_cctocl_static_summer <- expand.grid(roedeer = c("Roe Deer Absent","Roe Deer Present"), habitat=c("Broad Leaved Forest", "Coniferous Forest", "Mixed Forest", "Human Land-Use Areas"), popden_std=popden_d$popden_std)
occ_prob_cctocl_static_summer$popden <- rep(popden_d$popden, each = 8)
occ_prob_cctocl_static_summer$pred <- c(pm.occ_cctocl_static_summer)
occ_prob_cctocl_static_summer$lower<- c(CRI.occ_cctocl_static_summer[1,,,])
occ_prob_cctocl_static_summer$upper<- c(CRI.occ_cctocl_static_summer[2,,,])
occ_prob_cctocl_static_summer$season <- "Summer"
occ_prob_cctocl_static_summer$species <- "Roe Deer"

head(occ_prob_cctocl_static_summer)

## 2.2.2. Winter ----
##########################################################

# prey p/a X habitat X length.out.popden X mcmc list

occ_cctocl_static_winter <- array(NA, dim = c(2,4,length(popden_d$popden_std),cctocl_static_winter_072023$mcmc.info$n.samples))

for (i in 1:length(popden_d$popden_std)){ # for popden
  for(a in 1:length(unique(cl_ss_dat_winter$area))){ # for area
    
    # prey absent
    
    occ_cctocl_static_winter[1,1,i,] <- plogis(cctocl_static_winter_072023$sims.list$alphapsi.cl[,1,1] # BL
                                               + cctocl_static_winter_072023$sims.list$betapsi.cl[,1,1]*popden_d$popden_std[i]
                                               + cctocl_static_winter_072023$sims.list$randompsi.cl[,1,a]) 
    
    occ_cctocl_static_winter[1,2,i,] <- plogis(cctocl_static_winter_072023$sims.list$alphapsi.cl[,1,2] # CF
                                               + cctocl_static_winter_072023$sims.list$betapsi.cl[,1,2]*popden_d$popden_std[i]
                                               + cctocl_static_winter_072023$sims.list$randompsi.cl[,1,a]) 
    
    occ_cctocl_static_winter[1,3,i,] <- plogis(cctocl_static_winter_072023$sims.list$alphapsi.cl[,1,3] # MF
                                               + cctocl_static_winter_072023$sims.list$betapsi.cl[,1,3]*popden_d$popden_std[i]
                                               + cctocl_static_winter_072023$sims.list$randompsi.cl[,1,a]) 
    
    occ_cctocl_static_winter[1,4,i,] <- plogis(cctocl_static_winter_072023$sims.list$alphapsi.cl[,1,4] # O
                                               + cctocl_static_winter_072023$sims.list$betapsi.cl[,1,4]*popden_d$popden_std[i]
                                               + cctocl_static_winter_072023$sims.list$randompsi.cl[,1,a]) 
    
    # prey present
    
    occ_cctocl_static_winter[2,1,i,] <- plogis(cctocl_static_winter_072023$sims.list$alphapsi.cl[,2,1] # BL
                                               + cctocl_static_winter_072023$sims.list$betapsi.cl[,2,1]*popden_d$popden_std[i]
                                               + cctocl_static_winter_072023$sims.list$randompsi.cl[,2,a]) 
    
    occ_cctocl_static_winter[2,2,i,] <- plogis(cctocl_static_winter_072023$sims.list$alphapsi.cl[,2,2] # CF
                                               + cctocl_static_winter_072023$sims.list$betapsi.cl[,2,2]*popden_d$popden_std[i]
                                               + cctocl_static_winter_072023$sims.list$randompsi.cl[,2,a]) 
    
    occ_cctocl_static_winter[2,3,i,] <- plogis(cctocl_static_winter_072023$sims.list$alphapsi.cl[,2,3] # MF
                                               + cctocl_static_winter_072023$sims.list$betapsi.cl[,2,3]*popden_d$popden_std[i]
                                               + cctocl_static_winter_072023$sims.list$randompsi.cl[,2,a]) 
    
    occ_cctocl_static_winter[2,4,i,] <- plogis(cctocl_static_winter_072023$sims.list$alphapsi.cl[,2,4] # O
                                               + cctocl_static_winter_072023$sims.list$betapsi.cl[,2,4]*popden_d$popden_std[i]
                                               + cctocl_static_winter_072023$sims.list$randompsi.cl[,2,a]) 
  }
}

str(occ_cctocl_static_winter)  

# then we take the mean of the mcmc list
pm.occ_cctocl_static_winter <- apply(occ_cctocl_static_winter, c(1,2,3), mean)
str(pm.occ_cctocl_static_winter)

# then calculate the credible intervals
CRI.occ_cctocl_static_winter <- apply(occ_cctocl_static_winter, c(1,2,3), function(x) quantile(x, c(0.025, 0.975)))
str(CRI.occ_cctocl_static_winter)

occ_prob_cctocl_static_winter <- expand.grid(roedeer = c("Roe Deer Absent","Roe Deer Present"), habitat=c("Broad Leaved Forest", "Coniferous Forest", "Mixed Forest", "Human Land-Use Areas"), popden_std=popden_d$popden_std)
occ_prob_cctocl_static_winter$popden <- rep(popden_d$popden, each = 8)
occ_prob_cctocl_static_winter$pred <- c(pm.occ_cctocl_static_winter)
occ_prob_cctocl_static_winter$lower<- c(CRI.occ_cctocl_static_winter[1,,,])
occ_prob_cctocl_static_winter$upper<- c(CRI.occ_cctocl_static_winter[2,,,])
occ_prob_cctocl_static_winter$season <- "Winter"
occ_prob_cctocl_static_winter$species <- "Roe Deer"

head(occ_prob_cctocl_static_winter)

## 2.2.3. Merge and plot occupancy ----
##########################################################

occ_prob_cctocl_static <- rbind(occ_prob_cctocl_static_summer, occ_prob_cctocl_static_winter)
occ_prob_cctocl_static$habitat_wrapped <- stringr::str_wrap(occ_prob_cctocl_static$habitat, width = 3)

occ_prob_cctocl_static_plot <- ggplot(occ_prob_cctocl_static)+
  geom_ribbon(aes(x= popden, ymin = lower, ymax = upper, fill=season, group=season),alpha=0.3)+
  geom_line(aes(popden, pred, col=season, group=season), lwd=2) +
  facet_grid(roedeer ~ habitat_wrapped) + 
  ylim(c(0,1)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())+
  ylab(expression(paste("Occupancy Probabilitiy", " ", (psi)))) +
  xlab(expression("Rural human population density "(per~'4'~km^2))) +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=12, face = "bold"),
        panel.spacing = unit(1, "lines"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size=10))

occ_prob_cctocl_static_plot


## 2.3. Site occupancy states ----
##########################################################

fin_occu_total_cctocl_static <- expand.grid(state=c("RoeDeer(1) Wolf(1)", "RoeDeer(1) Wolf(0)", "RoeDeer(0) Wolf(1)", "RoeDeer(0) Wolf(0)"), season = c("Summer","Winter"))
fin_occu_total_cctocl_static$state_wrapped <- stringr::str_wrap(fin_occu_total_cctocl_static$state, width = 1)

fin_occu_total_cctocl_static$pred[1] <- cctocl_static_summer_072023$mean$nDW_S
fin_occu_total_cctocl_static$pred[2] <- cctocl_static_summer_072023$mean$nDw_S
fin_occu_total_cctocl_static$pred[3] <- cctocl_static_summer_072023$mean$ndW_S
fin_occu_total_cctocl_static$pred[4] <- cctocl_static_summer_072023$mean$ndw_S
fin_occu_total_cctocl_static$pred[5] <- cctocl_static_winter_072023$mean$nDW_W
fin_occu_total_cctocl_static$pred[6] <- cctocl_static_winter_072023$mean$nDw_W
fin_occu_total_cctocl_static$pred[7] <- cctocl_static_winter_072023$mean$ndW_W
fin_occu_total_cctocl_static$pred[8] <- cctocl_static_winter_072023$mean$ndw_W

fin_occu_total_cctocl_static$lower[1] <- cctocl_static_summer_072023$q2.5$nDW_S
fin_occu_total_cctocl_static$lower[2] <- cctocl_static_summer_072023$q2.5$nDw_S
fin_occu_total_cctocl_static$lower[3] <- cctocl_static_summer_072023$q2.5$ndW_S
fin_occu_total_cctocl_static$lower[4] <- cctocl_static_summer_072023$q2.5$ndw_S
fin_occu_total_cctocl_static$lower[5] <- cctocl_static_winter_072023$q2.5$nDW_W
fin_occu_total_cctocl_static$lower[6] <- cctocl_static_winter_072023$q2.5$nDw_W
fin_occu_total_cctocl_static$lower[7] <- cctocl_static_winter_072023$q2.5$ndW_W
fin_occu_total_cctocl_static$lower[8] <- cctocl_static_winter_072023$q2.5$ndw_W

fin_occu_total_cctocl_static$upper[1] <- cctocl_static_summer_072023$q97.5$nDW_S
fin_occu_total_cctocl_static$upper[2] <- cctocl_static_summer_072023$q97.5$nDw_S
fin_occu_total_cctocl_static$upper[3] <- cctocl_static_summer_072023$q97.5$ndW_S
fin_occu_total_cctocl_static$upper[4] <- cctocl_static_summer_072023$q97.5$ndw_S
fin_occu_total_cctocl_static$upper[5] <- cctocl_static_winter_072023$q97.5$nDW_W
fin_occu_total_cctocl_static$upper[6] <- cctocl_static_winter_072023$q97.5$nDw_W
fin_occu_total_cctocl_static$upper[7] <- cctocl_static_winter_072023$q97.5$ndW_W
fin_occu_total_cctocl_static$upper[8] <- cctocl_static_winter_072023$q97.5$ndw_W

head(fin_occu_total_cctocl_static)


fin_occu_total_cctocl_static_plot <- ggplot(fin_occu_total_cctocl_static) + 
  geom_bar(aes(state_wrapped, pred, fill=season),stat = "identity", position="dodge", colour="black")+ 
  geom_errorbar(aes(state_wrapped, ymin=lower, ymax=upper, group=season), position=position_dodge(.9),width=0.3) +
  scale_y_continuous(breaks=seq(0,210, by=10)) +
  theme_classic() + 
  theme(panel.grid.minor = element_blank())+
  ylab("Number of occupied sites") +
  xlab('States') +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.text = element_text(size=12, face = "bold"),
        panel.spacing = unit(2, "lines"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.width = unit(1,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size=14))

fin_occu_total_cctocl_static_plot

## 2.4. Species Interaction Factor  ----
##########################################################

str(cctocl_static_summer_072023$sims.list$SIF)
str(cctocl_static_summer_072023$mean$SIF)

SIF_cctocl_static <- data.frame(site = rep(1:171,2))

SIF_cctocl_static$pred <- c(cctocl_static_summer_072023$mean$SIF,
                            cctocl_static_winter_072023$mean$SIF)

SIF_cctocl_static$lower <- c(cctocl_static_summer_072023$q2.5$SIF,
                             cctocl_static_winter_072023$q2.5$SIF)

SIF_cctocl_static$upper <- c(cctocl_static_summer_072023$q97.5$SIF,
                             cctocl_static_winter_072023$q97.5$SIF)

SIF_cctocl_static$season[1:171] <- "Summer" 
SIF_cctocl_static$season[172:342] <- "Winter"

head(SIF_cctocl_static)

SIF_cctocl_static_plot <- ggplot(SIF_cctocl_static, aes(season, pred, fill=season)) + 
  geom_violin(aes(group=season,), width=1,position = position_dodge(.9))+ 
  geom_boxplot(aes(group=season),colour="black", width=.1,size=.5) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=seq(0.5,2,0.5),limits = c(0.5,2))+
  #scale_x_continuous(breaks=seq(1,171,by=1)) +
  ylab("Species Interaction Factor") +
  #xlab('Season') +
  theme(axis.text = element_text(size=16),
        axis.title.y = element_text(size=18),
        axis.title.x = element_blank(),
        legend.text = element_text(size=22),
        legend.title = element_blank(),
        legend.position = "none",
        legend.direction = "horizontal")

SIF_cctocl_static_plot


## 2.5. Co-occurrence diagram  ----
##########################################################

# THIS NEEDS TO BE RUN AFTER CALCULATING THE FIN_OCCU PARAMETERS (CHECK SECTION 2.3)

###########################################################################
# colors in euler diagrams should indicate these:
# PURPLE = CO-OCCURRENCE
# BLUE = ONLY PREDATOR
# ORANGE = ONLY PREDATOR
# WHITE = NONE
###########################################################################

cctocl_static_finoccu <- fin_occu_total_cctocl_static
cctocl_static_finoccu$pred <- round(cctocl_static_finoccu$pred)
cctocl_static_finoccu$percentage[1:4] <- round((cctocl_static_finoccu$pred[1:4]*100)/sum(cctocl_static_finoccu$pred[1:4])) # to make sure its same with euler, summer
cctocl_static_finoccu$percentage[5:8] <- round((cctocl_static_finoccu$pred[5:8]*100)/sum(cctocl_static_finoccu$pred[5:8])) # to make sure its same with euler, winter

# summer


cctocl_summer <- euler(c("Total" = cctocl_static_finoccu$pred[4],"RoeDeer" = 0, "GWolf" = 0,
                         "Total&RoeDeer" = cctocl_static_finoccu$pred[2],
                         "Total&GWolf" = cctocl_static_finoccu$pred[3],
                         "RoeDeer&GWolf" = 0,
                         "Total&RoeDeer&GWolf" = cctocl_static_finoccu$pred[1]))

cctocl_summer_euler <- plot(cctocl_summer,
                            fills = list(fill = c("white","#e66101","#0571b0",NA,NA,NA,"#7570b3")),
                            edges = list(lwd=3),
                            quantities = list(type="percent", cex=4, font=4, hjust=0.5),
                            #labels = list(labels=c("Unoccupied", "Roe Deer","Gray Wolf"), cex=1.8, font=2),
                            labels = list(labels=c("","",""), cex=1.8, font=2),
                            #main=list(label="Summer", cex=3),
                            padding=list(unit(10, "line")))
#legend = list(labels=c("Unoccupied", "Roe Deer","Gray Wolf"), side="top",vgap=1 ,cex=2))

cctocl_summer_euler

# winter

cctocl_winter <- euler(c("Total" = cctocl_static_finoccu$pred[8],"RoeDeer" = 0, "GWolf" = 0,
                         "Total&RoeDeer" = cctocl_static_finoccu$pred[6],
                         "Total&GWolf" = cctocl_static_finoccu$pred[7],
                         "RoeDeer&GWolf" = 0,
                         "Total&RoeDeer&GWolf" = cctocl_static_finoccu$pred[5]))

cctocl_winter_euler <- plot(cctocl_winter,
                            fills = list(fill = c("white","#e66101","#0571b0",NA,NA,NA,"#7570b3")),
                            edges = list(lwd=3),
                            quantities = list(type="percent", cex=4, font=4, hjust=0.5),
                            #labels = list(labels=c("Unoccupied", "Roe Deer","Gray Wolf"), cex=1.8, font=2),
                            labels = list(labels=c("","",""), cex=1.8, font=2),
                            #main=list(label="Summer", cex=3),
                            padding=list(unit(10, "line")))
#legend = list(labels=c("Unoccupied", "Roe Deer","Gray Wolf"), side="top",vgap=1 ,cex=2))

cctocl_winter_euler
