############################################################################
# This script is for the analysis of large-mammals in North Anatolia.
# with dynamic multispecies occupancy models

# multi-season multi-species occupancy models.
# According to  Waddle et al 2010, Ecological Applications
# Capreolus capreolus as dominant and Canis lupus as sub-ordinant

# Script for model diagnostics

# Interaction model 1

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
load("path_to_data/cctocl_dynamic_multisp_data_github.RData")

# load model output
load("path_to_model_output/cctocl_dynamic_int1_model_out.RData")

###########################################################################
# 2. Model diagnostics
###########################################################################

## Convergence and distribution check ----
###############################################

ni = cctocl_dynamic_int1_072023$mcmc.info$n.iter

# checking Rhat values, if there ara any over 1.1
hist(cctocl_dynamic_int1_072023$summary[,8])
length(which(cctocl_dynamic_int1_072023$summary[,8]>1.1))

rhats <- data.frame(params = names(cctocl_dynamic_int1_072023$summary[,8]), rhat = as.numeric(cctocl_dynamic_int1_072023$summary[,8]))

# traceplots
MCMCtrace(cctocl_dynamic_int1_072023$samples, params = c("alphapsi.cc", "betapsi.cc",
                                                                 "alphagamma.cc", "betagamma.cc",
                                                                 "alphaeps.cc","betaeps.cc", 
                                                                 "alphap.cc", "betap.cc"), pdf=F, ind = T, Rhat = T, n.eff = T, iter=ni)

MCMCtrace(cctocl_dynamic_int1_072023$samples, params = "randompsi.cc", pdf=F, ind = T, Rhat = T, n.eff = T,iter=ni)
MCMCtrace(cctocl_dynamic_int1_072023$samples, params = "randomgamma.cc", pdf=F, ind = T, Rhat = T, n.eff = T,iter=ni)
MCMCtrace(cctocl_dynamic_int1_072023$samples, params = "randomeps.cc", pdf=F, ind = T, Rhat = T, n.eff = T,iter=ni)
MCMCtrace(cctocl_dynamic_int1_072023$samples, params = "randomp.cc", pdf=F, ind = T, Rhat = T, n.eff = T,iter=ni)

MCMCtrace(cctocl_dynamic_int1_072023$samples, params = c("alphapsi.cl", "betapsi.cl",
                                                                 "alphagamma.cl", "betagamma.cl",
                                                                 "alphaeps.cl","betaeps.cl", 
                                                                 "alphap.cl", "betap.cl"), pdf=F, ind = T, Rhat = T, n.eff = T, iter=ni)

MCMCtrace(cctocl_dynamic_int1_072023$samples, params = "randompsi.cl", pdf=F, ind = T, Rhat = T, n.eff = T,iter=ni)
MCMCtrace(cctocl_dynamic_int1_072023$samples, params = "randomgamma.cl", pdf=F, ind = T, Rhat = T, n.eff = T,iter=ni)
MCMCtrace(cctocl_dynamic_int1_072023$samples, params = "randomeps.cl", pdf=F, ind = T, Rhat = T, n.eff = T,iter=ni)
MCMCtrace(cctocl_dynamic_int1_072023$samples, params = "randomp.cl", pdf=F, ind = T, Rhat = T, n.eff = T,iter=ni)

# caterpillar plots

par(mfrow=c(1,1))
MCMCplot(cctocl_dynamic_int1_072023$samples, params = c("alphapsi.cc", "betapsi.cc",
                                                                "alphagamma.cc", "betagamma.cc",
                                                                "alphaeps.cc","betaeps.cc", 
                                                                "alphap.cc", "betap.cc"), ref_ovl = T)

par(mfrow=c(2,2))
MCMCplot(cctocl_dynamic_int1_072023$samples, params = "randompsi.cc", ref_ovl = T)
MCMCplot(cctocl_dynamic_int1_072023$samples, params = "randomgamma.cc", ref_ovl = T)
MCMCplot(cctocl_dynamic_int1_072023$samples, params = "randomeps.cc", ref_ovl = T)
MCMCplot(cctocl_dynamic_int1_072023$samples, params = "randomp.cc", ref_ovl = T)

par(mfrow=c(1,1))
MCMCplot(cctocl_dynamic_int1_072023$samples, params = c("alphapsi.cl", "betapsi.cl",
                                                                "alphagamma.cl", "betagamma.cl",
                                                                "alphaeps.cl","betaeps.cl", 
                                                                "alphap.cl", "betap.cl"), ref_ovl = T)

par(mfrow=c(2,2))
MCMCplot(cctocl_dynamic_int1_072023$samples, params = "randompsi.cl", ref_ovl = T)
MCMCplot(cctocl_dynamic_int1_072023$samples, params = "randomgamma.cl", ref_ovl = T)
MCMCplot(cctocl_dynamic_int1_072023$samples, params = "randomeps.cl", ref_ovl = T)
MCMCplot(cctocl_dynamic_int1_072023$samples, params = "randomp.cl", ref_ovl = T)

par(mfrow=c(1,1))

# posterior distributions and traceplots (better visualization)

color_scheme_set("purple")

plot_title <- ggtitle("Roe deer with gray wolf | int1 | Posterior distributions",
                      "medians and 95% intervals")

panel_background <- panel_bg(fill = 'white')


# parameters prey ---------------------------------------------------------

mcmc_areas(cctocl_dynamic_int1_072023$samples,
           pars = c("alphapsi.cc[1]", "alphapsi.cc[2]", "alphapsi.cc[3]", "alphapsi.cc[4]",
                    "betapsi.cc[1]", "betapsi.cc[2]", "betapsi.cc[3]", "betapsi.cc[4]",
                    "alphagamma.cc", "betagamma.cc[1]", "betagamma.cc[2]", "betagamma.cc[3]",
                    "alphaeps.cc", "betaeps.cc[1]", "betaeps.cc[2]", "betaeps.cc[3]",
                    "alphap.cc[1]", "alphap.cc[2]", "alphap.cc[3]", "alphap.cc[4]",
                    "betap.cc[1]", "betap.cc[2]", "betap.cc[3]", "betap.cc[4]"),
           
           area_method = "equal height",
           prob = 0.95) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(alpha, "_", psi[1])),
                              expression(paste(alpha, "_", psi[2])),
                              expression(paste(alpha, "_", psi[3])),
                              expression(paste(alpha, "_", psi[4])),
                              expression(paste(beta, "_", psi[1])),
                              expression(paste(beta, "_", psi[2])),
                              expression(paste(beta, "_", psi[3])),
                              expression(paste(beta, "_", psi[4])),
                              expression(paste(alpha, "_", gamma)),
                              expression(paste(beta, "_", gamma[1])),
                              expression(paste(beta, "_", gamma[2])),
                              expression(paste(beta, "_", gamma[3])),
                              expression(paste(beta, "_", gamma[4])),
                              expression(paste(alpha, "_", epsilon)),
                              expression(paste(beta, "_", epsilon[1])),
                              expression(paste(beta, "_", epsilon[2])),
                              expression(paste(beta, "_", epsilon[3])),
                              expression(paste(beta, "_", epsilon[4])),
                              expression(paste(alpha, "_", p[1])),
                              expression(paste(alpha, "_", p[2])),
                              expression(paste(alpha, "_", p[3])),
                              expression(paste(alpha, "_", p[4])),
                              expression(paste(beta, "_", p[1])),
                              expression(paste(beta, "_", p[2])),
                              expression(paste(beta, "_", p[3])),
                              expression(paste(beta, "_", p[4])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))


# parameters predator -----------------------------------------------------


mcmc_areas(cctocl_dynamic_int1_072023$samples,
           pars = c("alphapsi.cl[1,1]", "alphapsi.cl[1,2]", "alphapsi.cl[1,3]", "alphapsi.cl[1,4]",
                    "betapsi.cl[1,1]", "betapsi.cl[1,2]", "betapsi.cl[1,3]", "betapsi.cl[1,4]",
                    "alphagamma.cl[1]", "betagamma.cl[1,1]", "betagamma.cl[2,1]", "betagamma.cl[3,1]",
                    "alphaeps.cl[1]", "betaeps.cl[1,1]", "betaeps.cl[2,1]", "betaeps.cl[3,1]",
                    "alphap.cl[1,1]", "alphap.cl[1,2]", "alphap.cl[1,3]", "alphap.cl[1,4]",
                    "betap.cl[1,1]", "betap.cl[1,2]", "betap.cl[1,3]", "betap.cl[1,4]"),
           
           area_method = "equal height",
           prob = 0.95) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(alpha, "_", psi[1])),
                              expression(paste(alpha, "_", psi[2])),
                              expression(paste(alpha, "_", psi[3])),
                              expression(paste(alpha, "_", psi[4])),
                              expression(paste(beta, "_", psi[1])),
                              expression(paste(beta, "_", psi[2])),
                              expression(paste(beta, "_", psi[3])),
                              expression(paste(beta, "_", psi[4])),
                              expression(paste(alpha, "_", gamma)),
                              expression(paste(beta, "_", gamma[1])),
                              expression(paste(beta, "_", gamma[2])),
                              expression(paste(beta, "_", gamma[3])),
                              expression(paste(beta, "_", gamma[4])),
                              expression(paste(alpha, "_", epsilon)),
                              expression(paste(beta, "_", epsilon[1])),
                              expression(paste(beta, "_", epsilon[2])),
                              expression(paste(beta, "_", epsilon[3])),
                              expression(paste(beta, "_", epsilon[4])),
                              expression(paste(alpha, "_", p[1])),
                              expression(paste(alpha, "_", p[2])),
                              expression(paste(alpha, "_", p[3])),
                              expression(paste(alpha, "_", p[4])),
                              expression(paste(beta, "_", p[1])),
                              expression(paste(beta, "_", p[2])),
                              expression(paste(beta, "_", p[3])),
                              expression(paste(beta, "_", p[4])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))

mcmc_areas(cctocl_dynamic_int1_072023$samples,
           pars = c("alphapsi.cl[2,1]", "alphapsi.cl[2,2]", "alphapsi.cl[2,3]", "alphapsi.cl[2,4]",
                    "betapsi.cl[2,1]", "betapsi.cl[2,2]", "betapsi.cl[2,3]", "betapsi.cl[2,4]",
                    "alphagamma.cl[2]", "betagamma.cl[1,2]", "betagamma.cl[2,2]", "betagamma.cl[3,2]",
                    "alphaeps.cl[2]", "betaeps.cl[1,2]", "betaeps.cl[2,2]", "betaeps.cl[3,2]",
                    "alphap.cl[2,1]", "alphap.cl[2,2]", "alphap.cl[2,3]", "alphap.cl[2,4]",
                    "betap.cl[2,1]", "betap.cl[2,2]", "betap.cl[2,3]", "betap.cl[2,4]"),
           
           area_method = "equal height",
           prob = 0.95) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(alpha, "_", psi[1])),
                              expression(paste(alpha, "_", psi[2])),
                              expression(paste(alpha, "_", psi[3])),
                              expression(paste(alpha, "_", psi[4])),
                              expression(paste(beta, "_", psi[1])),
                              expression(paste(beta, "_", psi[2])),
                              expression(paste(beta, "_", psi[3])),
                              expression(paste(beta, "_", psi[4])),
                              expression(paste(alpha, "_", gamma)),
                              expression(paste(beta, "_", gamma[1])),
                              expression(paste(beta, "_", gamma[2])),
                              expression(paste(beta, "_", gamma[3])),
                              expression(paste(beta, "_", gamma[4])),
                              expression(paste(alpha, "_", epsilon)),
                              expression(paste(beta, "_", epsilon[1])),
                              expression(paste(beta, "_", epsilon[2])),
                              expression(paste(beta, "_", epsilon[3])),
                              expression(paste(beta, "_", epsilon[4])),
                              expression(paste(alpha, "_", p[1])),
                              expression(paste(alpha, "_", p[2])),
                              expression(paste(alpha, "_", p[3])),
                              expression(paste(alpha, "_", p[4])),
                              expression(paste(beta, "_", p[1])),
                              expression(paste(beta, "_", p[2])),
                              expression(paste(beta, "_", p[3])),
                              expression(paste(beta, "_", p[4])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))



# re prey -----------------------------------------------------------------


mcmc_areas(cctocl_dynamic_int1_072023$samples,
           pars = c("randompsi.cc[1]", "randompsi.cc[2]", "randompsi.cc[3]", "randompsi.cc[4]", "randompsi.cc[5]", 
                    "randompsi.cc[6]", "randompsi.cc[7]", "randompsi.cc[8]", "randompsi.cc[9]", "randompsi.cc[10]"),
           area_method = "equal height",
           prob = 0.95) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(re(area), "_", psi[1])),
                              expression(paste(re(area), "_", psi[2])),
                              expression(paste(re(area), "_", psi[3])),
                              expression(paste(re(area), "_", psi[4])),
                              expression(paste(re(area), "_", psi[5])),
                              expression(paste(re(area), "_", psi[6])),
                              expression(paste(re(area), "_", psi[7])),
                              expression(paste(re(area), "_", psi[8])),
                              expression(paste(re(area), "_", psi[9])),
                              expression(paste(re(area), "_", psi[10])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))

mcmc_areas(cctocl_dynamic_int1_072023$samples,
           pars = c("randomgamma.cc[1]", "randomgamma.cc[2]", "randomgamma.cc[3]", "randomgamma.cc[4]", "randomgamma.cc[5]", 
                    "randomgamma.cc[6]", "randomgamma.cc[7]", "randomgamma.cc[8]", "randomgamma.cc[9]", "randomgamma.cc[10]"),
           area_method = "equal height",
           prob = 0.95) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(re(area), "_", gamma[1])),
                              expression(paste(re(area), "_", gamma[2])),
                              expression(paste(re(area), "_", gamma[3])),
                              expression(paste(re(area), "_", gamma[4])),
                              expression(paste(re(area), "_", gamma[5])),
                              expression(paste(re(area), "_", gamma[6])),
                              expression(paste(re(area), "_", gamma[7])),
                              expression(paste(re(area), "_", gamma[8])),
                              expression(paste(re(area), "_", gamma[9])),
                              expression(paste(re(area), "_", gamma[10])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))

mcmc_areas(cctocl_dynamic_int1_072023$samples,
           pars = c("randomeps.cc[1]", "randomeps.cc[2]", "randomeps.cc[3]", "randomeps.cc[4]", "randomeps.cc[5]", 
                    "randomeps.cc[6]", "randomeps.cc[7]", "randomeps.cc[8]", "randomeps.cc[9]", "randomeps.cc[10]"),
           area_method = "equal height",
           prob = 0.95) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(re(area), "_", epsilon[1])),
                              expression(paste(re(area), "_", epsilon[2])),
                              expression(paste(re(area), "_", epsilon[3])),
                              expression(paste(re(area), "_", epsilon[4])),
                              expression(paste(re(area), "_", epsilon[5])),
                              expression(paste(re(area), "_", epsilon[6])),
                              expression(paste(re(area), "_", epsilon[7])),
                              expression(paste(re(area), "_", epsilon[8])),
                              expression(paste(re(area), "_", epsilon[9])),
                              expression(paste(re(area), "_", epsilon[10])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))

mcmc_areas(cctocl_dynamic_int1_072023$samples,
           pars = c("randomp.cc[1]", "randomp.cc[2]", "randomp.cc[3]", "randomp.cc[4]", "randomp.cc[5]", 
                    "randomp.cc[6]", "randomp.cc[7]", "randomp.cc[8]", "randomp.cc[9]", "randomp.cc[10]"),
           area_method = "equal height",
           prob = 0.95) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(re(area), "_", p[1])),
                              expression(paste(re(area), "_", p[2])),
                              expression(paste(re(area), "_", p[3])),
                              expression(paste(re(area), "_", p[4])),
                              expression(paste(re(area), "_", p[5])),
                              expression(paste(re(area), "_", p[6])),
                              expression(paste(re(area), "_", p[7])),
                              expression(paste(re(area), "_", p[8])),
                              expression(paste(re(area), "_", p[9])),
                              expression(paste(re(area), "_", p[10])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))



# re predator -------------------------------------------------------------


mcmc_areas(cctocl_dynamic_int1_072023$samples,
           pars = c("randompsi.cl[1,1]", "randompsi.cl[1,2]", "randompsi.cl[1,3]", "randompsi.cl[1,4]", "randompsi.cl[1,5]", 
                    "randompsi.cl[1,6]", "randompsi.cl[1,7]", "randompsi.cl[1,8]", "randompsi.cl[1,9]", "randompsi.cl[1,10]"),
           area_method = "equal height",
           prob = 0.95) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(re(area), "_", psi[1])),
                              expression(paste(re(area), "_", psi[2])),
                              expression(paste(re(area), "_", psi[3])),
                              expression(paste(re(area), "_", psi[4])),
                              expression(paste(re(area), "_", psi[5])),
                              expression(paste(re(area), "_", psi[6])),
                              expression(paste(re(area), "_", psi[7])),
                              expression(paste(re(area), "_", psi[8])),
                              expression(paste(re(area), "_", psi[9])),
                              expression(paste(re(area), "_", psi[10])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))

mcmc_areas(cctocl_dynamic_int1_072023$samples,
           pars = c("randomgamma.cl[1,1]", "randomgamma.cl[1,2]", "randomgamma.cl[1,3]", "randomgamma.cl[1,4]", "randomgamma.cl[1,5]", 
                    "randomgamma.cl[1,6]", "randomgamma.cl[1,7]", "randomgamma.cl[1,8]", "randomgamma.cl[1,9]", "randomgamma.cl[1,10]"),
           area_method = "equal height",
           prob = 0.95) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(re(area), "_", gamma[1])),
                              expression(paste(re(area), "_", gamma[2])),
                              expression(paste(re(area), "_", gamma[3])),
                              expression(paste(re(area), "_", gamma[4])),
                              expression(paste(re(area), "_", gamma[5])),
                              expression(paste(re(area), "_", gamma[6])),
                              expression(paste(re(area), "_", gamma[7])),
                              expression(paste(re(area), "_", gamma[8])),
                              expression(paste(re(area), "_", gamma[9])),
                              expression(paste(re(area), "_", gamma[10])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))

mcmc_areas(cctocl_dynamic_int1_072023$samples,
           pars = c("randomeps.cl[1,1]", "randomeps.cl[1,2]", "randomeps.cl[1,3]", "randomeps.cl[1,4]", "randomeps.cl[1,5]", 
                    "randomeps.cl[1,6]", "randomeps.cl[1,7]", "randomeps.cl[1,8]", "randomeps.cl[1,9]", "randomeps.cl[1,10]"),
           area_method = "equal height",
           prob = 0.95) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(re(area), "_", epsilon[1])),
                              expression(paste(re(area), "_", epsilon[2])),
                              expression(paste(re(area), "_", epsilon[3])),
                              expression(paste(re(area), "_", epsilon[4])),
                              expression(paste(re(area), "_", epsilon[5])),
                              expression(paste(re(area), "_", epsilon[6])),
                              expression(paste(re(area), "_", epsilon[7])),
                              expression(paste(re(area), "_", epsilon[8])),
                              expression(paste(re(area), "_", epsilon[9])),
                              expression(paste(re(area), "_", epsilon[10])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))

mcmc_areas(cctocl_dynamic_int1_072023$samples,
           pars = c("randomp.cl[1,1]", "randomp.cl[1,2]", "randomp.cl[1,3]", "randomp.cl[1,4]", "randomp.cl[1,5]", 
                    "randomp.cl[1,6]", "randomp.cl[1,7]", "randomp.cl[1,8]", "randomp.cl[1,9]", "randomp.cl[1,10]"),
           area_method = "equal height",
           prob = 0.95) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(re(area), "_", p[1])),
                              expression(paste(re(area), "_", p[2])),
                              expression(paste(re(area), "_", p[3])),
                              expression(paste(re(area), "_", p[4])),
                              expression(paste(re(area), "_", p[5])),
                              expression(paste(re(area), "_", p[6])),
                              expression(paste(re(area), "_", p[7])),
                              expression(paste(re(area), "_", p[8])),
                              expression(paste(re(area), "_", p[9])),
                              expression(paste(re(area), "_", p[10])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))


mcmc_areas(cctocl_dynamic_int1_072023$samples,
           pars = c("randompsi.cl[2,1]", "randompsi.cl[2,2]", "randompsi.cl[2,3]", "randompsi.cl[2,4]", "randompsi.cl[2,5]", 
                    "randompsi.cl[2,6]", "randompsi.cl[2,7]", "randompsi.cl[2,8]", "randompsi.cl[2,9]", "randompsi.cl[2,10]"),
           area_method = "equal height",
           prob = 0.95) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(re(area), "_", psi[1])),
                              expression(paste(re(area), "_", psi[2])),
                              expression(paste(re(area), "_", psi[3])),
                              expression(paste(re(area), "_", psi[4])),
                              expression(paste(re(area), "_", psi[5])),
                              expression(paste(re(area), "_", psi[6])),
                              expression(paste(re(area), "_", psi[7])),
                              expression(paste(re(area), "_", psi[8])),
                              expression(paste(re(area), "_", psi[9])),
                              expression(paste(re(area), "_", psi[10])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))

mcmc_areas(cctocl_dynamic_int1_072023$samples,
           pars = c("randomgamma.cl[2,1]", "randomgamma.cl[2,2]", "randomgamma.cl[2,3]", "randomgamma.cl[2,4]", "randomgamma.cl[2,5]", 
                    "randomgamma.cl[2,6]", "randomgamma.cl[2,7]", "randomgamma.cl[2,8]", "randomgamma.cl[2,9]", "randomgamma.cl[2,10]"),
           area_method = "equal height",
           prob = 0.95) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(re(area), "_", gamma[1])),
                              expression(paste(re(area), "_", gamma[2])),
                              expression(paste(re(area), "_", gamma[3])),
                              expression(paste(re(area), "_", gamma[4])),
                              expression(paste(re(area), "_", gamma[5])),
                              expression(paste(re(area), "_", gamma[6])),
                              expression(paste(re(area), "_", gamma[7])),
                              expression(paste(re(area), "_", gamma[8])),
                              expression(paste(re(area), "_", gamma[9])),
                              expression(paste(re(area), "_", gamma[10])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))

mcmc_areas(cctocl_dynamic_int1_072023$samples,
           pars = c("randomeps.cl[2,1]", "randomeps.cl[2,2]", "randomeps.cl[2,3]", "randomeps.cl[2,4]", "randomeps.cl[2,5]", 
                    "randomeps.cl[2,6]", "randomeps.cl[2,7]", "randomeps.cl[2,8]", "randomeps.cl[2,9]", "randomeps.cl[2,10]"),
           area_method = "equal height",
           prob = 0.95) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(re(area), "_", epsilon[1])),
                              expression(paste(re(area), "_", epsilon[2])),
                              expression(paste(re(area), "_", epsilon[3])),
                              expression(paste(re(area), "_", epsilon[4])),
                              expression(paste(re(area), "_", epsilon[5])),
                              expression(paste(re(area), "_", epsilon[6])),
                              expression(paste(re(area), "_", epsilon[7])),
                              expression(paste(re(area), "_", epsilon[8])),
                              expression(paste(re(area), "_", epsilon[9])),
                              expression(paste(re(area), "_", epsilon[10])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))

mcmc_areas(cctocl_dynamic_int1_072023$samples,
           pars = c("randomp.cl[2,1]", "randomp.cl[2,2]", "randomp.cl[2,3]", "randomp.cl[2,4]", "randomp.cl[2,5]", 
                    "randomp.cl[2,6]", "randomp.cl[2,7]", "randomp.cl[2,8]", "randomp.cl[2,9]", "randomp.cl[2,10]"),
           area_method = "equal height",
           prob = 0.95) + 
  plot_title + panel_background +
  scale_y_discrete(labels = c(expression(paste(re(area), "_", p[1])),
                              expression(paste(re(area), "_", p[2])),
                              expression(paste(re(area), "_", p[3])),
                              expression(paste(re(area), "_", p[4])),
                              expression(paste(re(area), "_", p[5])),
                              expression(paste(re(area), "_", p[6])),
                              expression(paste(re(area), "_", p[7])),
                              expression(paste(re(area), "_", p[8])),
                              expression(paste(re(area), "_", p[9])),
                              expression(paste(re(area), "_", p[10])))) + 
  theme_classic() +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))
