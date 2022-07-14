#=============================================================
#Example use of samsim for time varying simulation evaluation
#using the coho data as it is compliant with the most recent
#samSim updates
#Catarina Wor
# March 2022 
#=============================================================

#TODO

#add estimation routines - in progress
 

#install samsim 
#remotes::install_github("Pacific-salmon-assess/samSim", force=TRUE)

library(ggplot2)
library(devtools)
library(gridExtra)
library(dplyr)
library(here)
source("dlm-wrapper.R")
source("sim_est_func.R")



#here::here()
## Load relevant input data
# Simulation run parameters describing different scenarios
simPar <- read.csv("../data/samsimIFcoho/cohoSimPars_test.csv")
#simPar <- read.csv("../data/samsimIFcoho/cohoSimPars_ER.csv")
# CU-specific parameters
cuPar <- read.csv("../data/samsimIFcoho/cohoCUPars_test.csv")

# Stock-recruit and catch data that are used to populate the simulation priming
# period
srDat  <- read.csv("../data/samsimIFcoho/cohoRecDatTrim.csv")
 
# Posterior values of  CU-specific stock recruitment parameters for Ricker and
# Larkin models; when available, passed and used to calculate alpha, beta and
# sigma parameters rather than passing point values 
ricPars <- read.csv("../data/samsimIFcoho/cohoRickerpars.csv")

corrmat <- read.csv("../data/samsimIFcoho/cohoCorrMat.csv", header=F)


out <- sim_est(simPar=simPar[1,], cuPar=cuPar, srDat=srDat, ricPars=ricPars,
 corrmat=corrmat, outDir="example", iteration=2, simrun=FALSE)

dfall <- calcpbias(simData=out$simData,estData= out$estData)
unique(dfall$estmodel)

dfc<-dfall[dfall$convergence==0,]
unique(dfc$estmodel)

myplot<- list()

scen <- unique(dfc$scenario)
dfc[dfc$parameter=="beta",]

dfc[dfc$parameter=="beta"&dfc$estmodel=="TMB Kalman Filter a vary",]
a=1
for(a in seq_along(scen)){
  dfp<-dfc[dfc$scenario==scen[a],]
  
  #mypl<-
  ggplot(dfp,aes(x=estmodel, y=pbias)) +
    geom_boxplot()+
    coord_cartesian(ylim = c(-1,1))+
    geom_hline(yintercept=0) +
    theme_bw(14)+
    ggtitle(scen[a])+
    facet_wrap(~parameter)

 ggsave(
      filename = paste0("../plots/",scen[a],"_pbias.pdf"), 
      plot = mypl,#marrangeGrob( mypl, nrow=1, ncol=1), 
      width = 12, height = 5
    )


}

ggsave(
      filename = "../plots/scenarios_pbias.pdf", 
      plot = marrangeGrob( myplot, nrow=1, ncol=1), 
      width = 12, height = 5
    )

ggplot(dfc,aes(x=estmodel, y=pbias)) +
geom_boxplot()+
coord_cartesian(ylim = c(-1,1))+
geom_hline(yintercept=0) +
theme_bw(14)+
#stat_summary(fun.data = give.n, geom = "text", hjust = 0.5,
#    vjust = 0.9)+
facet_wrap(~parameter)

#add more estimation routines
#



#=============================================
#cmdstanr

library(cmdstanr);library(loo)

dat<-out$simData[[1]]  %>% 
        dplyr::filter( iteration==n, CU==1) %>% 
        dplyr::select(byr, obsSpawners,obsRecruits ) %>%
        dplyr::rename(spwn=obsSpawners, rec=obsRecruits)


data=list(R_S = log(dat$rec/dat$spwn),
          N=nrow(dat),
          TT=as.numeric(factor(seq_len(nrow(dat)))),
          S=c(dat$spwn))

fit_mle$mle()
system(paste("cp stan/ricker_linear_varying_a.stan", paste0(cmdstan_path(),"/stanmodels")))
file1 <- file.path(cmdstan_path(),'stanmodels',"ricker_linear_varying_a.stan")

mod <- cmdstan_model(file1)
#mcmc
fit<- mod$sample(
    data = data,
    seed = 123, 
    chains = 6, 
    parallel_chains = 6,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 500,
    adapt_delta = 0.99,
    max_treedepth = 20 # print update every 500 iters
  )

params2<- fit$draws(format='df',variables=c('log_a','b','log_b'))
summary(params2)

#nll
fit_mle <- mod$optimize(data = data, seed = 123) 
fit_mle$summary()

