############################################################################
# This script is for the analysis of large-mammals in North Anatolia.
# with dynamic multispecies occupancy models

# single-season multi-species occupancy models.
# According to  Waddle et al 2010, Ecological Applications
# Capreolus capreolus as dominant and Canis lupus as sub-ordinant

# Script for model diagnostics

# Summer

# Date: July 2023
# Author: Dilsad Dagtekin
###########################################################################

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
load("path_to_data/cctocl_static_summer_data_github.RData")

# load model output
load("path_to_model_output/cctocl_static_summer_model_out.RData")

###########################################################################
# 2. Model diagnostics
###########################################################################

## Convergence and distribution check ----
###############################################

ni = cctocl_static_summer_072023$mcmc.info$n.iter

# checking Rhat values, if there ara any over 1.1
hist(cctocl_static_summer_072023$summary[,8])
length(which(cctocl_static_summer_072023$summary[,8]>1.1))

rhats <- data.frame(params = names(cctocl_static_summer_072023$summary[,8]), rhat = as.numeric(cctocl_static_summer_072023$summary[,8]))

# traceplots
MCMCtrace(cctocl_static_summer_072023$samples, params = c("alphapsi.cc", "betapsi.cc", "alphap.cc"), pdf=F, ind = T, Rhat = T, n.eff = T, iter=ni)
MCMCtrace(cctocl_static_summer_072023$samples, params = c("randompsi.cc", "randomp.cc"), pdf=F, ind = T, Rhat = T, n.eff = T,iter=ni)

MCMCtrace(cctocl_static_summer_072023$samples, params = c("alphapsi.cl", "betapsi.cl", "alphap.cl"), pdf=F, ind = T, Rhat = T, n.eff = T, iter=ni)
MCMCtrace(cctocl_static_summer_072023$samples, params = c("randompsi.cl", "randomp.cl"), pdf=F, ind = T, Rhat = T, n.eff = T,iter=ni)

par(mfrow=c(1,1))
# caterpillar plots
MCMCplot(cctocl_static_summer_072023$samples, params = c("alphapsi.cc", "betapsi.cc", "alphap.cc"), ref_ovl = T)
MCMCplot(cctocl_static_summer_072023$samples, params = c("randompsi.cc", "randomp.cc"), ref_ovl = T)

MCMCplot(cctocl_static_summer_072023$samples, params = c("alphapsi.cl", "betapsi.cl", "alphap.cl"), ref_ovl = T)
MCMCplot(cctocl_static_summer_072023$samples, params = c("randompsi.cl", "randomp.cl"), ref_ovl = T)

# posterior distributions and traceplots (better visualization)

color_scheme_set("purple")

plot_title <- ggtitle("Roe deer with gray wolf | Posterior distributions",
                      "medians and 95% intervals")

panel_background <- panel_bg(fill = 'white')

mcmc_areas(cctocl_static_summer_072023$samples,
           pars = c("alphapsi.cc[1]", "alphapsi.cc[2]", "alphapsi.cc[3]", "alphapsi.cc[4]",
                    "betapsi.cc[1]", "betapsi.cc[2]", "betapsi.cc[3]", "betapsi.cc[4]",
                    "alphap.cc[1]", "alphap.cc[2]", "alphap.cc[3]", "alphap.cc[4]"),
           area_method = "equal height",
           prob = 0.95,
           rhat = c(rhats[1:4,2], rhats[9:12,2], rhats[5:8,2])) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(alpha, "_", psi[1])),
                              expression(paste(alpha, "_", psi[2])),
                              expression(paste(alpha, "_", psi[3])),
                              expression(paste(alpha, "_", psi[4])),
                              expression(paste(beta, "_", psi[1])),
                              expression(paste(beta, "_", psi[2])),
                              expression(paste(beta, "_", psi[3])),
                              expression(paste(beta, "_", psi[4])),
                              expression(paste(alpha, "_", p[1])),
                              expression(paste(alpha, "_", p[2])),
                              expression(paste(alpha, "_", p[3])),
                              expression(paste(alpha, "_", p[4])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))

mcmc_areas(cctocl_static_summer_072023$samples,
           pars = c("alphapsi.cl[1,1]", "alphapsi.cl[1,2]", "alphapsi.cl[1,3]", "alphapsi.cl[1,4]",
                    "betapsi.cl[1,1]", "betapsi.cl[1,2]", "betapsi.cl[1,3]", "betapsi.cl[1,4]",
                    "alphap.cl[1,1]", "alphap.cl[1,2]", "alphap.cl[1,3]", "alphap.cl[1,4]"),
           area_method = "equal height",
           prob = 0.95,
           rhat = c(rhats[1:4,2], rhats[9:12,2], rhats[5:8,2])) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(alpha, "_", psi[1])),
                              expression(paste(alpha, "_", psi[2])),
                              expression(paste(alpha, "_", psi[3])),
                              expression(paste(alpha, "_", psi[4])),
                              expression(paste(beta, "_", psi[1])),
                              expression(paste(beta, "_", psi[2])),
                              expression(paste(beta, "_", psi[3])),
                              expression(paste(beta, "_", psi[4])),
                              expression(paste(alpha, "_", p[1])),
                              expression(paste(alpha, "_", p[2])),
                              expression(paste(alpha, "_", p[3])),
                              expression(paste(alpha, "_", p[4])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))

mcmc_areas(cctocl_static_summer_072023$samples,
           pars = c("alphapsi.cl[2,1]", "alphapsi.cl[2,2]", "alphapsi.cl[2,3]", "alphapsi.cl[2,4]",
                    "betapsi.cl[2,1]", "betapsi.cl[2,2]", "betapsi.cl[2,3]", "betapsi.cl[2,4]",
                    "alphap.cl[2,1]", "alphap.cl[2,2]", "alphap.cl[2,3]", "alphap.cl[2,4]"),
           area_method = "equal height",
           prob = 0.95,
           rhat = c(rhats[1:4,2], rhats[9:12,2], rhats[5:8,2])) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(alpha, "_", psi[1])),
                              expression(paste(alpha, "_", psi[2])),
                              expression(paste(alpha, "_", psi[3])),
                              expression(paste(alpha, "_", psi[4])),
                              expression(paste(beta, "_", psi[1])),
                              expression(paste(beta, "_", psi[2])),
                              expression(paste(beta, "_", psi[3])),
                              expression(paste(beta, "_", psi[4])),
                              expression(paste(alpha, "_", p[1])),
                              expression(paste(alpha, "_", p[2])),
                              expression(paste(alpha, "_", p[3])),
                              expression(paste(alpha, "_", p[4])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))

mcmc_areas(cctocl_static_summer_072023$samples,
           pars = c("randompsi.cc[1]", "randompsi.cc[2]", "randompsi.cc[3]", "randompsi.cc[4]", "randompsi.cc[5]", 
                    "randompsi.cc[6]", "randompsi.cc[7]", "randompsi.cc[8]", "randompsi.cc[9]", "randompsi.cc[10]", 
                    "randomp.cc[1]", "randomp.cc[2]", "randomp.cc[3]", "randomp.cc[4]", "randomp.cc[5]", 
                    "randomp.cc[6]", "randomp.cc[7]", "randomp.cc[8]", "randomp.cc[9]", "randomp.cc[10]"),
           area_method = "equal height",
           prob = 0.95) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(eps, "_", psi[1])),
                              expression(paste(eps, "_", psi[2])),
                              expression(paste(eps, "_", psi[3])),
                              expression(paste(eps, "_", psi[4])),
                              expression(paste(eps, "_", psi[5])),
                              expression(paste(eps, "_", psi[6])),
                              expression(paste(eps, "_", psi[7])),
                              expression(paste(eps, "_", psi[8])),
                              expression(paste(eps, "_", psi[9])),
                              expression(paste(eps, "_", psi[10])),
                              expression(paste(eps, "_", p[1])),
                              expression(paste(eps, "_", p[2])),
                              expression(paste(eps, "_", p[3])),
                              expression(paste(eps, "_", p[4])),
                              expression(paste(eps, "_", p[5])),
                              expression(paste(eps, "_", p[6])),
                              expression(paste(eps, "_", p[7])),
                              expression(paste(eps, "_", p[8])),
                              expression(paste(eps, "_", p[9])),
                              expression(paste(eps, "_", p[10])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))

mcmc_areas(cctocl_static_summer_072023$samples,
           pars = c("randompsi.cl[1,1]", "randompsi.cl[1,2]", "randompsi.cl[1,3]", "randompsi.cl[1,4]", "randompsi.cl[1,5]", 
                    "randompsi.cl[1,6]", "randompsi.cl[1,7]", "randompsi.cl[1,8]", "randompsi.cl[1,9]", "randompsi.cl[1,10]", 
                    "randomp.cl[1,1]", "randomp.cl[1,2]", "randomp.cl[1,3]", "randomp.cl[1,4]", "randomp.cl[1,5]", 
                    "randomp.cl[1,6]", "randomp.cl[1,7]", "randomp.cl[1,8]", "randomp.cl[1,9]", "randomp.cl[1,10]"),
           area_method = "equal height",
           prob = 0.95) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(eps, "_", psi[1])),
                              expression(paste(eps, "_", psi[2])),
                              expression(paste(eps, "_", psi[3])),
                              expression(paste(eps, "_", psi[4])),
                              expression(paste(eps, "_", psi[5])),
                              expression(paste(eps, "_", psi[6])),
                              expression(paste(eps, "_", psi[7])),
                              expression(paste(eps, "_", psi[8])),
                              expression(paste(eps, "_", psi[9])),
                              expression(paste(eps, "_", psi[10])),
                              expression(paste(eps, "_", p[1])),
                              expression(paste(eps, "_", p[2])),
                              expression(paste(eps, "_", p[3])),
                              expression(paste(eps, "_", p[4])),
                              expression(paste(eps, "_", p[5])),
                              expression(paste(eps, "_", p[6])),
                              expression(paste(eps, "_", p[7])),
                              expression(paste(eps, "_", p[8])),
                              expression(paste(eps, "_", p[9])),
                              expression(paste(eps, "_", p[10])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))


mcmc_areas(cctocl_static_summer_072023$samples,
           pars = c("randompsi.cl[2,1]", "randompsi.cl[2,2]", "randompsi.cl[2,3]", "randompsi.cl[2,4]", "randompsi.cl[2,5]", 
                    "randompsi.cl[2,6]", "randompsi.cl[2,7]", "randompsi.cl[2,8]", "randompsi.cl[2,9]", "randompsi.cl[2,10]", 
                    "randomp.cl[2,1]", "randomp.cl[2,2]", "randomp.cl[2,3]", "randomp.cl[2,4]", "randomp.cl[2,5]", 
                    "randomp.cl[2,6]", "randomp.cl[2,7]", "randomp.cl[2,8]", "randomp.cl[2,9]", "randomp.cl[2,10]"),
           area_method = "equal height",
           prob = 0.95) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(eps, "_", psi[1])),
                              expression(paste(eps, "_", psi[2])),
                              expression(paste(eps, "_", psi[3])),
                              expression(paste(eps, "_", psi[4])),
                              expression(paste(eps, "_", psi[5])),
                              expression(paste(eps, "_", psi[6])),
                              expression(paste(eps, "_", psi[7])),
                              expression(paste(eps, "_", psi[8])),
                              expression(paste(eps, "_", psi[9])),
                              expression(paste(eps, "_", psi[10])),
                              expression(paste(eps, "_", p[1])),
                              expression(paste(eps, "_", p[2])),
                              expression(paste(eps, "_", p[3])),
                              expression(paste(eps, "_", p[4])),
                              expression(paste(eps, "_", p[5])),
                              expression(paste(eps, "_", p[6])),
                              expression(paste(eps, "_", p[7])),
                              expression(paste(eps, "_", p[8])),
                              expression(paste(eps, "_", p[9])),
                              expression(paste(eps, "_", p[10])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))

