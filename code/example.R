#=============================================================
#Example use of samsim for time varying simulation evaluation
#using the coho data as it is compliant with the most recent
#samSim updates
#Catarina Wor
# March 2022 
#=============================================================

#TODO
#add estimation routines
#

#install samsim 
#remotes::install_github("Pacific-salmon-assess/samSim", ref="timevar", force=TRUE)

#install samest
#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst')

library(samEst)
library(samSim)
library(ggplot2)
library(devtools)
library(gridExtra)
library(dplyr)
library(here)
source("sgen_functions.R")
source("utils.R")
#TODO estimate only for 40 yrs of data.

#here::here()
## Load relevant input data
# Simulation run parameters describing different scenarios
#simPar <- read.csv("../data/samsimIFcoho/cohoSimPars_test.csv")
simPar <- read.csv("../data/samsimHarCk/harcnkSimPars.csv")
# CU-specific parameters
#cuPar <- read.csv("../data/samsimIFcoho/cohoCUPars_test.csv")
cuPar <- read.csv("../data/samsimHarCk/harcnkCUPars.csv")

# Stock-recruit and catch data that are used to populate the simulation priming
# period
#srDat  <- read.csv("../data/samsimIFcoho/cohoRecDatTrim.csv")
 
# Posterior values of  CU-specific stock recruitment parameters for Ricker and
# Larkin models; when available, passed and used to calculate alpha, beta and
# sigma parameters rather than passing point values 
#ricPars <- read.csv("../data/samsimIFcoho/cohoRickerpars.csv")

#corrmat <- read.csv("../data/samsimIFcoho/cohoCorrMat.csv", header=F)

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simPar$species),sep = "_"))

## First check to ensure that a single scenario can be run (only a small number
# of trials necessary)
plotscn <- TRUE
p <- list()
simData <- list()

genericRecoverySim(simPar=simPar[1,], cuPar=cuPar, catchDat=NULL, srDat=NULL,
            variableCU=FALSE, ricPars=NULL , larkPars=NULL,cuCustomCorrMat= NULL,
            outDir="test", nTrials=100, makeSubDirs=TRUE, random=FALSE)

simData <- readRDS(paste0("../test/SamSimOutputs/simData/", simPar$nameOM[1],"/",simPar$scenario[1],"/",
                         paste(simPar$nameOM[1],"_", simPar$nameMP[1], "_", "CUsrDat.RData",sep="")))$srDatout
  

simData<-simData[simData$CU==1,]

lfodf<-matrix(NA,nrow=length(unique(simData$iteration)),ncol=14,
  dimnames = list(unique(simData$iteration),
    c("simple", "autocorr", "rwa_lastparam","rwa_last3paramavg","rwa_last5paramavg",
      "rwb_lastparam","rwb_last3paramavg","rwb_last5paramavg",
      "hmm_regime_pick","hmm_regime_average","hmma_regime_pick",
      "hmma_regime_average","hmmb_regime_pick","hmmb_regime_average")))

simest<-list()
rmse<-list()

#compiled Bayesian models
simple_mod <- sr_mod(type='static', ac=FALSE, par='n', loglik=FALSE, modelcode=TRUE)
simpleac_mod <- sr_mod(type='static', ac=TRUE, par='n', loglik=FALSE, modelcode=TRUE)
rwa_mod <- sr_mod(type='rw',ac=FALSE,par="a",loglik=FALSE, modelcode=TRUE)
rwb_mod <- sr_mod(type='rw',ac=FALSE,par="b",loglik=FALSE, modelcode=TRUE)
rwab_mod <- sr_mod(type='rw',ac=FALSE,par="both",loglik=FALSE, modelcode=TRUE)
hmma_mod<-sr_mod(type='hmm',ac=FALSE,par="a",loglik=FALSE, modelcode=TRUE)
hmmb_mod<-sr_mod(type='hmm',ac=FALSE,par="b",loglik=FALSE, modelcode=TRUE)
hmmab_mod<-sr_mod(type='hmm',ac=FALSE,par="both",loglik=FALSE, modelcode=TRUE)



for(u in unique(simData$iteration)){
  #u=1
  dat<-simData[simData$iteration==u,]
  dat<-dat[dat$year>(max(dat$year)-46),]

  dat <- dat[!is.na(dat$obsRecruits),]
  df <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))
 

  p <- ricker_TMB(data=df)
  pac <- ricker_TMB(data=df, AC=TRUE)
  ptva <- ricker_rw_TMB(data=df,tv.par='a')
  ptvb <- ricker_rw_TMB(data=df, tv.par='b')
  ptvab <- ricker_rw_TMB(data=df, tv.par='both')
  phmma <- ricker_hmm_TMB(data=df, tv.par='a')
  phmmb <- ricker_hmm_TMB(data=df, tv.par='b')
  phmm  <- ricker_hmm_TMB(data=df, tv.par='both')

  b <- ricker_stan(data=df,iter = 2000,sm_ext=simple_mod)
  #ricker autocorr
  bac <- ricker_stan(data=df,iter = 2000, AC=TRUE, sm_ext=simpleac_mod)
  #ricker tva
  btva <- ricker_rw_stan(data=df, par="a",iter = 2000, sm_ext=rwa_mod)
  #ricker tvb
  btvb <- ricker_rw_stan(data=df, par="b",iter = 2000, sm_ext=rwb_mod)
  #ricker tvab
  btvab <- ricker_rw_stan(data=df, par="both",iter = 2000, sm_ext=rwab_mod)
  
  #ricker tvhmma
  bhmma <- ricker_hmm_stan(data=df, par="a",iter = 2000, sm_ext=hmma_mod)
  #ricker tvhmmb
  bhmmb <- ricker_hmm_stan(data=df, par="b",iter = 2000, sm_ext=hmmb_mod)
  #ricker tvhmmab
  bhmmab <- ricker_hmm_stan(data=df, par="both",iter = 2000, sm_ext=hmmab_mod)
  bhmmab[1:7]
  names(btvab)
  #a
  dfa<- data.frame(parameter="alpha",
    iteration=u,
    method=rep(c(rep("MLE",10),rep("MCMC",10)),each=nrow(df)),
    model=rep(c("simple","autocorr","rwa","rwb","rwab","hmma_regime","hmma_average",
      "hmmb_regime","hmmab_regime","hmmab_average","simple","autocorr",
      "rwa","rwb","rwab","hmma_regime","hmma_average",
      "hmmb_regime","hmmab_regime","hmmab_average"),each=nrow(df)),
    by=rep(dat$year,20),
    sim=rep(dat$alpha,20),
    est=c(rep(p$alpha,nrow(df)),
      rep(pac$alpha,nrow(df)),
      ptva$alpha,
      rep(ptvb$alpha,nrow(df)),
      ptvab$alpha,
      phmma$alpha[phmma$regime],
      c(phmma$alpha%*%phmma$probregime),
      rep(phmmb$alpha,nrow(df)),
      phmm$alpha[phmm$regime],
      phmm$alpha%*%phmm$probregime,
      rep(b$alpha,nrow(df)),
      rep(bac$alpha,nrow(df)),
      btva$alpha[-1],
      rep(btvb$alpha,nrow(df)),
      btvab$alpha[-1],
      bhmma$alpha_regime,
      bhmma$alpha_wgt,
      rep(bhmmb$alpha,nrow(df)),
      bhmmab$alpha_regime,
      bhmmab$alpha_wgt
      ),
     convergence=c(rep(c(p$model$convergence + p$conv_problem,
      pac$model$convergence + pac$conv_problem,
      ptva$model$convergence + ptva$conv_problem,
      ptvb$model$convergence + ptvb$conv_problem,
      ptvab$model$convergence + ptvab$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      as.numeric(abs(b$mcmcsummary["log_a","Rhat"]-1)>.1),
      as.numeric(abs(bac$mcmcsummary["log_a","Rhat"]-1)>.1)
      ),each=nrow(df)),
     as.numeric(abs(btva$mcmcsummary[grep("log_a\\[",rownames(btva$mcmcsummary)),"Rhat"]-1)>.1),
     rep(as.numeric(abs(btvb$mcmcsummary[grep("log_a",rownames(btvb$mcmcsummary)),"Rhat"]-1)>.1),nrow(df)),
     as.numeric(abs(btvab$mcmcsummary[grep("log_a\\[",rownames(btvab$mcmcsummary)),"Rhat"]-1)>.1),
     as.numeric(abs(bhmma$mcmcsummary[grep("log_a_t\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1),
     as.numeric(abs(bhmma$mcmcsummary[grep("log_a_wt\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1),
     rep(as.numeric(abs(bhmmb$mcmcsummary[grep("log_a",rownames(bhmmb$mcmcsummary)),"Rhat"]-1)>.1),nrow(df)),
     as.numeric(abs(bhmmab$mcmcsummary[grep("log_a_t\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1),
     as.numeric(abs(bhmmab$mcmcsummary[grep("log_a_wt\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1)
     ))
   dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
   
   rmsea<-aggregate((dfa$est-dfa$sim)^2, list(model=dfa$model,method=dfa$method), function(x)sqrt(mean(x)))
   rmsea$iteration<-u
   rmsea$parameter<-"alpha"
   rmsea$convergence <- c(p$model$convergence + p$conv_problem,
      pac$model$convergence + pac$conv_problem,
      ptva$model$convergence + ptva$conv_problem,
      ptvb$model$convergence + ptvb$conv_problem,
      ptvab$model$convergence + ptvab$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      as.numeric(abs(b$mcmcsummary["log_a","Rhat"]-1)>.1),
      as.numeric(abs(bac$mcmcsummary["log_a","Rhat"]-1)>.1),
      sum(as.numeric(abs(btva$mcmcsummary[grep("log_a\\[",rownames(btva$mcmcsummary)),"Rhat"]-1)>.1)),
      as.numeric(abs(btvb$mcmcsummary[grep("log_a",rownames(btvb$mcmcsummary)),"Rhat"]-1)>.1),
      sum(as.numeric(abs(btvab$mcmcsummary[grep("log_a\\[",rownames(btvab$mcmcsummary)),"Rhat"]-1)>.1)),
      sum(as.numeric(abs(bhmma$mcmcsummary[grep("log_a_t\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1)),
      sum(as.numeric(abs(bhmma$mcmcsummary[grep("log_a_wt\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1)),
      as.numeric(abs(bhmmb$mcmcsummary[grep("log_a",rownames(bhmmb$mcmcsummary)),"Rhat"]-1)>.1),
      sum(as.numeric(abs(bhmmab$mcmcsummary[grep("log_a_t\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1)),
      sum(as.numeric(abs(bhmmab$mcmcsummary[grep("log_a_wt\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1))  
      )
  


  #Smax
   dfsmax<- data.frame(parameter="Smax",
    iteration=u,
    method=rep(c(rep("MLE",10),rep("MCMC",10)),each=nrow(df)),
    model=rep(c("simple","autocorr","rwa","rwb","rwab","hmma_regime",
      "hmmb_regime","hmmb_average","hmmab_regime","hmmab_average",
      "simple","autocorr","rwa","rwb","rwab","hmma_regime",
      "hmmb_regime","hmmb_average","hmmab_regime","hmmab_average"),each=nrow(df)),
    by=rep(dat$year,20),
    sim=rep(1/dat$beta,20),
    est=c(rep(p$Smax,nrow(df)),
      rep(pac$Smax,nrow(df)),
      rep(ptva$Smax,nrow(df)),
      ptvb$Smax,
      ptvab$Smax,
      rep(phmma$Smax,nrow(df)),
      phmmb$Smax[phmmb$regime],
      phmmb$Smax%*%phmmb$probregime,
      phmm$Smax[phmm$regime],
      phmm$Smax%*%phmm$probregime,
      rep(b$Smax,nrow(df)),
      rep(bac$Smax,nrow(df)),
      rep(btva$Smax,nrow(df)),
      btvb$Smax,
      btvab$Smax, 
      rep(bhmma$Smax,nrow(df)),
      bhmmb$Smax_regime,
      bhmmb$Smax_wgt,
      bhmmab$Smax_regime,
      bhmmab$Smax_wgt),
     convergence=c(rep(c(p$model$convergence + p$conv_problem,
      pac$model$convergence + pac$conv_problem,
      ptva$model$convergence + ptva$conv_problem,
      ptvb$model$convergence + ptvb$conv_problem,
      ptvab$model$convergence + ptvab$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      as.numeric(abs(b$mcmcsummary["S_max","Rhat"]-1)>.1),
      as.numeric(abs(bac$mcmcsummary["S_max","Rhat"]-1)>.1),
      as.numeric(abs(btva$mcmcsummary["S_max","Rhat"]-1)>.1)
      ),each=nrow(df)),
      as.numeric(abs(btvb$mcmcsummary[grep("S_max",rownames(btvb$mcmcsummary)),"Rhat"]-1)>.1),
      as.numeric(abs(btvab$mcmcsummary[grep("S_max",rownames(btvab$mcmcsummary)),"Rhat"]-1)>.1),     
      rep(as.numeric(abs(bhmma$mcmcsummary[grep("S_max",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1),nrow(df)),
      as.numeric(abs(bhmmb$mcmcsummary[grep("S_max_t\\[",rownames(bhmmb$mcmcsummary)),"Rhat"]-1)>.1),
      as.numeric(abs(bhmmb$mcmcsummary[grep("S_max_wt\\[",rownames(bhmmb$mcmcsummary)),"Rhat"]-1)>.1),
      as.numeric(abs(bhmmab$mcmcsummary[grep("S_max_t\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1),
      as.numeric(abs(bhmmab$mcmcsummary[grep("S_max_wt\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1))
      )
   dfsmax$pbias<- ((dfsmax$est-dfsmax$sim)/dfsmax$sim)*100

   rmsesmax<-aggregate((dfsmax$est-dfsmax$sim)^2, 
    list(model=dfsmax$model,method=dfsmax$method), function(x)sqrt(mean(x)))
   rmsesmax$iteration<-u
   rmsesmax$parameter<-"Smax"
   rmsesmax$convergence <- c(p$model$convergence + p$conv_problem,
      pac$model$convergence + pac$conv_problem,
      ptva$model$convergence + ptva$conv_problem,
      ptvb$model$convergence + ptvb$conv_problem,
      ptvab$model$convergence + ptvab$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      as.numeric(abs(b$mcmcsummary["S_max","Rhat"]-1)>.1),
      as.numeric(abs(bac$mcmcsummary["S_max","Rhat"]-1)>.1),
      as.numeric(abs(btva$mcmcsummary["S_max","Rhat"]-1)>.1),
      sum(as.numeric(abs(btvb$mcmcsummary[grep("S_max",rownames(btvb$mcmcsummary)),"Rhat"]-1)>.1)),
      sum(as.numeric(abs(btvab$mcmcsummary[grep("S_max",rownames(btvab$mcmcsummary)),"Rhat"]-1)>.1)),
      as.numeric(abs(bhmma$mcmcsummary[grep("S_max",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1),
      sum(as.numeric(abs(bhmmb$mcmcsummary[grep("S_max_t\\[",rownames(bhmmb$mcmcsummary)),"Rhat"]-1)>.1)),
      sum(as.numeric(abs(bhmmb$mcmcsummary[grep("S_max_wt\\[",rownames(bhmmb$mcmcsummary)),"Rhat"]-1)>.1)),
      sum(as.numeric(abs(bhmmab$mcmcsummary[grep("S_max_t\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1)),
      sum(as.numeric(abs(bhmmab$mcmcsummary[grep("S_max_wt\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1)))
          
      
      

  #sigma
   dfsig<- data.frame(parameter="sigma",
    iteration=u,
    method=rep(c(rep("MLE",8),rep("MCMC",8)),each=nrow(df)),
    model=rep(c("simple","autocorr","rwa","rwb","rwab","hmma_regime",
      "hmmb_regime","hmmab_regime","simple","autocorr","rwa","rwb","rwab",
      "hmma_regime","hmmb_regime","hmmab_regime"),each=nrow(df)),
    by=rep(dat$year,16),
    sim=rep(dat$sigma,16),
    est=c(rep(p$sig,nrow(df)),
      rep(pac$sig,nrow(df)),
      rep(ptva$sig,nrow(df)),
      rep(ptvb$sig,nrow(df)),
      rep(ptvab$sig,nrow(df)),
      rep(phmma$sigma,nrow(df)),
      rep(phmmb$sigma,nrow(df)),
      rep(phmm$sigma,nrow(df)),
      rep(b$sigobs,nrow(df)),
      rep(bac$sigobs,nrow(df)),
      rep(btva$sigobs,nrow(df)),
      rep(btvb$sigobs,nrow(df)),
      rep(btvab$sigobs,nrow(df)),
      rep(bhmma$sigobs,nrow(df)),
      rep(bhmmb$sigobs,nrow(df)),
      rep(bhmmab$sigobs,nrow(df))),
     convergence=rep(c(p$model$convergence + p$conv_problem,
      pac$model$convergence + pac$conv_problem,
      ptva$model$convergence + ptva$conv_problem,
      ptvb$model$convergence + ptvb$conv_problem,
      ptvab$model$convergence + ptvab$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      as.numeric(abs(b$mcmcsummary["sigma_e","Rhat"]-1)>.1),
      as.numeric(abs(bac$mcmcsummary["sigma_e","Rhat"]-1)>.1),
      as.numeric(abs(btva$mcmcsummary["sigma_e","Rhat"]-1)>.1),
      as.numeric(abs(btvb$mcmcsummary["sigma_e","Rhat"]-1)>.1),
      as.numeric(abs(btvab$mcmcsummary["sigma_e","Rhat"]-1)>.1),
      as.numeric(abs(bhmma$mcmcsummary["sigma","Rhat"]-1)>.1),
      as.numeric(abs(bhmmb$mcmcsummary["sigma","Rhat"]-1)>.1),
      as.numeric(abs(bhmmab$mcmcsummary["sigma","Rhat"]-1)>.1)
      ),each=nrow(df)))
   dfsig$pbias<- ((dfsig$est-dfsig$sim)/dfsig$sim)*100

   rmsesig<-aggregate((dfsig$est-dfsig$sim)^2, list(model=dfsig$model,method=dfsig$method), function(x)sqrt(mean(x)))
   rmsesig$iteration<-u
   rmsesig$parameter<-"sigma"
   rmsesig$convergence <- c(p$model$convergence + p$conv_problem,
      pac$model$convergence + pac$conv_problem,
      ptva$model$convergence + ptva$conv_problem,
      ptvb$model$convergence + ptvb$conv_problem,
      ptvab$model$convergence + ptvab$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      as.numeric(abs(b$mcmcsummary["sigma_e","Rhat"]-1)>.1),
      as.numeric(abs(bac$mcmcsummary["sigma_e","Rhat"]-1)>.1),
      as.numeric(abs(btva$mcmcsummary["sigma_e","Rhat"]-1)>.1),
      as.numeric(abs(btvb$mcmcsummary["sigma_e","Rhat"]-1)>.1),
      as.numeric(abs(btvab$mcmcsummary["sigma_e","Rhat"]-1)>.1),
      as.numeric(abs(bhmma$mcmcsummary["sigma","Rhat"]-1)>.1),
      as.numeric(abs(bhmmb$mcmcsummary["sigma","Rhat"]-1)>.1),
      as.numeric(abs(bhmmab$mcmcsummary["sigma","Rhat"]-1)>.1)
      )
  
 
       
  #Smsy
  smsysim<-smsySolver(dat$alpha,dat$beta)

 
  dfsmsy<- data.frame(parameter="smsy",
    iteration=u,
    method=rep(c(rep("MLE",11),rep("MCMC",5)),each=nrow(df)),
    model=rep(c("simple","autocorr","rwa","rwb","rwab","hmma_regime","hmma_average",
      "hmmb_regime","hmmb_average","hmmab_regime","hmmab_average",
      "simple_b", "autocorr_b", "rwa_b","rwb_b","rwab_b"),each=nrow(df)),
    by=rep(dat$year,16),
    sim=rep(smsysim,16),
    est=c(rep(p$Smsy,nrow(df)),
          rep(pac$Smsy,nrow(df)),
          ptva$Smsy,
          ptvb$Smsy,
          ptvab$Smsy,
          phmma$Smsy[phmma$regime],
          smsySolver(c(phmma$alpha%*%phmma$probregime),phmma$beta),
          phmmb$Smsy[phmmb$regime],
          smsySolver(phmmb$alpha,phmmb$beta%*%phmmb$probregime),
          phmm$Smsy[phmm$regime],
          smsySolver(phmm$alpha%*%phmm$probregime,phmm$beta%*%phmm$probregime),
          rep(b$Smsy,nrow(df)),
          rep(bac$Smsy,nrow(df)),
          btva$Smsy,
          btvb$Smsy,
          btvab$Smsy
          ),    
     convergence=c(rep(c(p$model$convergence + p$conv_problem,
      pac$model$convergence + pac$conv_problem,
      ptva$model$convergence + ptva$conv_problem,
      ptvb$model$convergence + ptvb$conv_problem,
      ptvab$model$convergence + ptvab$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      as.numeric(abs(b$mcmcsummary["S_msy","Rhat"]-1)>.1),
      as.numeric(abs(bac$mcmcsummary["S_msy","Rhat"]-1)>.1)),each=nrow(df)),
      as.numeric(abs(btva$mcmcsummary[grep("S_msy",rownames(btva$mcmcsummary)),"Rhat"]-1)>.1),
      as.numeric(abs(btvb$mcmcsummary[grep("S_msy",rownames(btvb$mcmcsummary)),"Rhat"]-1)>.1),
      as.numeric(abs(btvab$mcmcsummary[grep("S_msy",rownames(btvab$mcmcsummary)),"Rhat"]-1)>.1))
  ) 
  
  dfsmsy$pbias<- ((dfsmsy$est-dfsmsy$sim)/dfsmsy$sim)*100

  rmsesmsy<-aggregate((dfsmsy$est-dfsmsy$sim)^2, list(model=dfsmsy$model,method=dfsmsy$method), function(x)sqrt(mean(x)))
  rmsesmsy$iteration<-u
  rmsesmsy$parameter<-"smsy"
  rmsesmsy$convergence <- c(p$model$convergence + p$conv_problem,
      pac$model$convergence + pac$conv_problem,
      ptva$model$convergence + ptva$conv_problem,
      ptvb$model$convergence + ptvb$conv_problem,
      ptvab$model$convergence + ptvab$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      as.numeric(abs(b$mcmcsummary["S_msy","Rhat"]-1)>.1),
      as.numeric(abs(bac$mcmcsummary["S_msy","Rhat"]-1)>.1),
      sum(as.numeric(abs(btva$mcmcsummary[grep("S_msy",rownames(btva$mcmcsummary)),"Rhat"]-1)>.1)),
      sum(as.numeric(abs(btvb$mcmcsummary[grep("S_msy",rownames(btvb$mcmcsummary)),"Rhat"]-1)>.1)),
      sum(as.numeric(abs(btvab$mcmcsummary[grep("S_msy",rownames(btvab$mcmcsummary)),"Rhat"]-1)>.1)))
      
  


  #Sgen
  #calc Bayesian Sgens
  
  sgen_b<-median(unlist(mapply(sGenSolver,a=b$samples[b$samples$parameters=="log_a","value"],
      b=b$samples[b$samples$parameters=="b","value"],
      Smsy=b$samples[b$samples$parameters=="S_msy","value"])))
  sgen_bac<-median(unlist(mapply(sGenSolver,a=bac$samples[bac$samples$parameters=="log_a","value"],
      b=bac$samples[bac$samples$parameters=="b","value"],
      Smsy=bac$samples[bac$samples$parameters=="S_msy","value"])))


  sgen_tva<-NULL
  sgen_tvb<-NULL
  sgen_tvab<-NULL
  for(j in seq_len(nrow(dat))){
    #tva
    atva<-btva$samples[grep(paste0("log_a\\[",j,"\\]"),btva$samples$parameters),]
    smsytva<-btva$samples[grep(paste0("S_msy\\[",j,"\\]"),btva$samples$parameters),]
    
    sgen_tva[j]<-median(unlist(mapply(sGenSolver,a=atva$value,
      b=btva$samples[btva$samples$parameters=="b","value"],
      Smsy=smsytva$value)))
  

  #tvb
    betatvb <- btvb$samples[grep(paste0("^b\\[",j,"\\]"),btvb$samples$parameters),]
    smsytvb<-btvb$samples[grep(paste0("S_msy\\[",j,"\\]"),btvb$samples$parameters),]

    sgen_tvb[j]<-median(unlist(mapply(sGenSolver,a=btvb$samples[btvb$samples$parameters=="log_a","value"],
      b=betatvb$value,Smsy=smsytvb$value)))

    
    #tvab    
    atvab<-btvab$samples[grep(paste0("log_a\\[",j,"\\]"),btvab$samples$parameters),]
    betatvab<-btvab$samples[grep(paste0("^b\\[",j,"\\]"),btvab$samples$parameters),]
    smsytvab<-btvab$samples[grep(paste0("S_msy\\[",j,"\\]"),btvab$samples$parameters),]

    sgen_tvab[j]<-median(unlist(mapply(sGenSolver,a=atvab$value,
      b=betatvab$value,Smsy=smsytvab$value)))

    #hmma
    
  }
 




  dfsgen <- data.frame(parameter="sgen",
    iteration=u,
    method=rep(c(rep("MLE",11),rep("MCMC",5)),each=nrow(df)),
    model=rep(c("simple","autocorr","rwa","rwb","rwab","hmma_regime","hmma_average",
      "hmmb_regime","hmmb_average","hmmab_regime","hmmab_average",
      "simple_b","autocorr_b","rwa_b","rwb_b","rwab_b"),each=nrow(df)),
    by=rep(dat$year,16),
    sim=rep(unlist(mapply(sGenSolver,a=dat$alpha,Smsy=smsysim, b=dat$beta)),12),
    est=c(unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="simple"],Smsy=dfsmsy$est[dfsmsy$model=="simple"], b=1/dfsmax$est[dfsmax$model=="simple"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="autocorr"],Smsy=dfsmsy$est[dfsmsy$model=="autocorr"], b=1/dfsmax$est[dfsmax$model=="autocorr"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="rwa"],Smsy=dfsmsy$est[dfsmsy$model=="rwa"], b=1/dfsmax$est[dfsmax$model=="rwa"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="rwb"],Smsy=dfsmsy$est[dfsmsy$model=="rwb"], b=1/dfsmax$est[dfsmax$model=="rwb"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="rwab"],Smsy=dfsmsy$est[dfsmsy$model=="rwab"], b=1/dfsmax$est[dfsmax$model=="rwab"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="hmma_regime"],Smsy=dfsmsy$est[dfsmsy$model=="hmma_regime"], b=1/dfsmax$est[dfsmax$model=="hmma_regime"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="hmma_average"],Smsy=dfsmsy$est[dfsmsy$model=="hmma_average"], b=1/dfsmax$est[dfsmax$model=="hmma_regime"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="hmmb_regime"],Smsy=dfsmsy$est[dfsmsy$model=="hmmb_regime"], b=1/dfsmax$est[dfsmax$model=="hmmb_regime"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="hmmb_regime"],Smsy=dfsmsy$est[dfsmsy$model=="hmmb_average"], b=1/dfsmax$est[dfsmax$model=="hmmb_average"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="hmmab_regime"],Smsy=dfsmsy$est[dfsmsy$model=="hmmab_regime"], b=1/dfsmax$est[dfsmax$model=="hmmab_regime"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="hmmab_average"],Smsy=dfsmsy$est[dfsmsy$model=="hmmab_average"], b=1/dfsmax$est[dfsmax$model=="hmmab_average"])),
        sgen_b,
        sgen_bac,
        sgen_tva,
        sgen_tvb,
        sgen_tvab
        ),
     convergence=c(rep(c(p$model$convergence + p$conv_problem,
      pac$model$convergence + pac$conv_problem,
      ptva$model$convergence + ptvab$conv_problem,
      ptvb$model$convergence + ptvab$conv_problem,
      ptvab$model$convergence + ptvab$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      sum(as.numeric(abs(b$mcmcsummary["log_a","Rhat"]-1)>.1),
          as.numeric(abs(b$mcmcsummary["b","Rhat"]-1)>.1)),
      sum(as.numeric(abs(bac$mcmcsummary["log_a","Rhat"]-1)>.1),
          as.numeric(abs(bac$mcmcsummary["b","Rhat"]-1)>.1)))
      ,each=nrow(df)),
      (as.numeric(abs(btva$mcmcsummary[grep("log_a\\[",rownames(btva$mcmcsummary)),"Rhat"]-1)>.1)+
          as.numeric(abs(btva$mcmcsummary["b","Rhat"]-1)>.1)),
      (as.numeric(abs(btvb$mcmcsummary["log_a","Rhat"]-1)>.1)+
          as.numeric(abs(btvb$mcmcsummary[grep("^b\\[",rownames(btvb$mcmcsummary)),"Rhat"]-1)>.1)),
      (as.numeric(abs(btvab$mcmcsummary[grep("log_a\\[",rownames(btvab$mcmcsummary)),"Rhat"]-1)>.1)+
          as.numeric(abs(btvab$mcmcsummary[grep("^b\\[",rownames(btvab$mcmcsummary)),"Rhat"]-1)>.1))))

   

   rmsesgen<-aggregate((dfsgen$est-dfsgen$sim)^2, list(model=dfsgen$model,method=dfsgen$method), function(x)sqrt(mean(x)))
   rmsesgen$iteration<-u
   rmsesgen$parameter<-"sgen"
   rmsesgen$convergence <- c(p$model$convergence + p$conv_problem,
      pac$model$convergence + pac$conv_problem,
      ptva$model$convergence + ptva$conv_problem,
      ptvb$model$convergence + ptvb$conv_problem,
      ptvab$model$convergence + ptvab$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      sum(as.numeric(abs(b$mcmcsummary["log_a","Rhat"]-1)>.1),
          as.numeric(abs(b$mcmcsummary["b","Rhat"]-1)>.1)),
      sum(as.numeric(abs(bac$mcmcsummary["log_a","Rhat"]-1)>.1),
          as.numeric(abs(bac$mcmcsummary["b","Rhat"]-1)>.1)),
      sum(as.numeric(abs(btva$mcmcsummary["log_a","Rhat"]-1)>.1),
          as.numeric(abs(btva$mcmcsummary["b","Rhat"]-1)>.1)),
      sum(as.numeric(abs(btvb$mcmcsummary["log_a","Rhat"]-1)>.1),
          as.numeric(abs(btvb$mcmcsummary["b","Rhat"]-1)>.1)),
      sum(as.numeric(abs(btvab$mcmcsummary["log_a","Rhat"]-1)>.1),
          as.numeric(abs(btvab$mcmcsummary["b","Rhat"]-1)>.1))
      )
  

  #umsy
  #1-gsl::lambert_W0(exp(1 - ptva$alpha))) /ptva$beta

    dfumsy<- data.frame(parameter="umsy",
    iteration=u,
    method=rep(c(rep("MLE",10),rep("MCMC",1)),each=nrow(df)),
    model=rep(c("simple","autocorr","rwa","rwb","rwab","hmma_regime","hmma_average",
      "hmmb_regime","hmmab_regime","hmmab_average","simple_b"),each=nrow(df)),
    by=rep(dat$year,11),
    sim=rep((.5 * dat$alpha - 0.07 * dat$alpha^2),11),
    est=c(rep((1 - gsl::lambert_W0(exp(1 - p$alpha))), nrow(df)),
          rep((1 - gsl::lambert_W0(exp(1 - pac$alpha))), nrow(df)),
          (1 - gsl::lambert_W0(exp(1 - ptva$alpha))),
          rep((1 - gsl::lambert_W0(exp(1 - ptvb$alpha))), nrow(df)),
          (1 - gsl::lambert_W0(exp(1 - ptvab$alpha))),
          (1 - gsl::lambert_W0(exp(1 - phmma$alpha[phmma$regime]))),
          (1 - gsl::lambert_W0(exp(1 - c(phmma$alpha%*%phmma$probregime)))),
          rep((1 - gsl::lambert_W0(exp(1 - phmmb$alpha))), nrow(df)),
          (1 - gsl::lambert_W0(exp(1 - phmm$alpha[phmm$regime]))),
          (1 - gsl::lambert_W0(exp(1 - c(phmm$alpha%*%phmm$probregime)))),
          rep((1 - gsl::lambert_W0(exp(1 - b$alpha))),nrow(df)),
          rep((1 - gsl::lambert_W0(exp(1 - bac$alpha))),nrow(df))),
    #c(rep(.5 * p$alpha - 0.07 * p$alpha^2,nrow(df)),
    #      rep(.5 * pac$alpha - 0.07 * pac$alpha^2,nrow(df)),
    #      .5 * ptva$alpha - 0.07 * ptva$alpha^2,
    #      rep(.5 * ptvb$alpha - 0.07 * ptvb$alpha^2,nrow(df)),
    #      .5 * ptvab$alpha - 0.07 * ptvab$alpha^2,
    #      .5 * phmma$alpha[phmma$regime] - 0.07 * phmma$alpha[phmma$regime]^2,
    #      .5 * c(phmma$alpha%*%phmma$probregime) - 0.07 * c(phmma$alpha%*%phmma$probregime)^2,
    #      rep(.5 * phmmb$alpha - 0.07 * phmmb$alpha^2,nrow(df)),
    #      .5 * phmm$alpha[phmm$regime] - 0.07 * phmm$alpha[phmm$regime]^2,
    #      .5 * c(phmm$alpha%*%phmm$probregime) - 0.07 * c(phmm$alpha%*%phmm$probregime)^2),
     convergence=rep(c(p$model$convergence + p$conv_problem,
      pac$model$convergence + pac$conv_problem,
      ptva$model$convergence+ ptva$conv_problem,
      ptvb$model$convergence+ ptvb$conv_problem,
      ptvab$model$convergence + ptvab$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      as.numeric(abs(b$mcmcsummary["log_a","Rhat"]-1)>.1),
      as.numeric(abs(bac$mcmcsummary["log_a","Rhat"]-1)>.1)
      ),each=nrow(df))
    )

    dfumsy$pbias<- ((dfumsy$est-dfumsy$sim)/dfumsy$sim)*100

    rmseumsy<-aggregate((dfumsy$est-dfumsy$sim)^2, list(model=dfumsy$model,method=dfumsy$method), function(x)sqrt(mean(x)))
    rmseumsy$iteration<-u
    rmseumsy$parameter<-"umsy"
    rmseumsy$convergence <- c(p$model$convergence + p$conv_problem,
      pac$model$convergence + pac$conv_problem,
      ptva$model$convergence + ptva$conv_problem,
      ptvb$model$convergence + ptvb$conv_problem,
      ptvab$model$convergence + ptvab$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmma$model$convergence + phmma$conv_problem,
      phmmb$model$convergence + phmmb$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      phmm$model$convergence + phmm$conv_problem,
      as.numeric(abs(b$mcmcsummary["log_a","Rhat"]-1)>.1),
      as.numeric(abs(bac$mcmcsummary["log_a","Rhat"]-1)>.1)
      )
  

   rmse[[u]]<-rbind(rmsea,rmseb,rmsesig,rmsesmsy,rmsesgen,rmseumsy)
   simest[[u]]<-rbind(dfa,dfb,dfsig,dfsmsy,dfsgen,dfumsy)
}
   


  

  
#=================================
#plots

dfpbias<-do.call("rbind",simest)
dfpbias<- dfpbias[dfpbias$convergence==0,]
dfpbias$model <- factor(dfpbias$model, levels=c("simple","simple_b", "autocorr", "rwa", "rwb",
           "rwab", "hmma_regime", "hmma_average", "hmmb_regime", "hmmb_average", 
           "hmmab_regime", "hmmab_average"))
dfpbias$method <- factor(dfpbias$method, levels=c("MLE","MCMC"))

head(dfpbias)

fig <- ggplot(dfpbias,aes(x=model,y=pbias)) +
geom_boxplot(aes(fill=method)) +
coord_cartesian(ylim = c(-100,100))+
geom_hline(yintercept=0) +
theme_bw(14)+ #theme(legend.position="none")+
facet_wrap(~parameter, scales="free_y")+
scale_fill_viridis_d(begin=.3, end=.9) +
stat_summary(fun.data = give.n, geom = "text", hjust = 0.5,
    vjust = -2)+
theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
fig


dfrmse<-do.call("rbind",rmse)
head(dfrmse)
dfrmse<- dfrmse[dfrmse$convergence==0,]
dfrmse$model <- factor(dfrmse$model, levels=c("simple", "autocorr", "rwa", "rwb",
           "rwab", "hmma_regime", "hmma_average", "hmmb_regime", "hmmb_average", "hmmab_regime", "hmmab_average"))



fig2 <- ggplot(dfrmse, aes(x=model,y=x)) +
geom_boxplot() +
theme_bw(14)+ theme(legend.position="none")+
facet_wrap(~parameter, scales="free_y")+
scale_colour_viridis_d() +
ylab("RMSE")+
stat_summary(fun.data = give.n, geom = "text", hjust = 0.5,
    vjust = -2)+
theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
fig2



for(u in unique(simData$iteration)){
   dat<-simData[simData$iteration==u,]
  dat<-dat[dat$year>(max(dat$year)-46),]

  dat <- dat[!is.na(dat$obsRecruits),]
  length(unique(dat$year))
  df <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))

#=======================
  #lfo comparison
  lfostatic<-tmb_mod_lfo_cv(data=df,tv.par='static')
  lfoac <- tmb_mod_lfo_cv(data=df,tv.par='staticAC')
  lfoalpha <- tmb_mod_lfo_cv(data=df,tv.par='alpha', siglfo="obs")
  lfobeta <- tmb_mod_lfo_cv(data=df,tv.par='beta', siglfo="obs")
  lfohmm <- tmb_mod_lfo_cv(data=df,tv.par='HMM')
  lfohmma <- tmb_mod_lfo_cv(data=df,tv.par='HMM_a')
  lfohmmb <- tmb_mod_lfo_cv(data=df,tv.par='HMM_b')

  lfodf[u,] <- c(sum(lfostatic), sum(lfoac), 
    sum(lfoalpha$lastparam), sum(lfoalpha$last3paramavg), sum(lfoalpha$last5paramavg), 
    sum(lfobeta$lastparam), sum(lfobeta$last3paramavg), sum(lfobeta$last5paramavg),
    sum(lfohmm$regime_pick),sum(lfohmm$regime_average),
    sum(lfohmma$regime_pick),sum(lfohmma$regime_average),
    sum(lfohmmb$regime_pick),sum(lfohmmb$regime_average))


}


head(lfodf)
dimnames(lfodf)[[2]]
lfo<-apply(lfodf,1,which.max)

lfochoice<-data.frame(
  chsnmod=dimnames(lfodf)[[2]][apply(lfodf,1,which.max)])

lfochoice$chsnmod<-factor(lfochoice$chsnmod, levels=c("simple", "autocorr", "rwa_lastparam", "rwa_last3paramavg", "rwa_last5paramavg",
 "rwb_lastparam", "rwb_last3paramavg", "rwb_last5paramavg", "hmm_regime_pick", "hmm_regime_average", 
 "hmma_regime_pick", "hmma_regime_average", "hmmb_regime_pick", "hmmb_regime_average"))



ggplot(lfochoice) +  
 geom_bar(aes(chsnmod))+
theme_bw(14)




for(i in 1:11){#seq_len(nrow(simPar))){

  genericRecoverySim(simPar=simPar[i,], cuPar=cuPar, catchDat=NULL, srDat=srDat,
            variableCU=FALSE, ricPars=ricPars , larkPars=NULL,cuCustomCorrMat= corrmat,
            outDir="example", nTrials=2, makeSubDirs=TRUE, random=FALSE)
  #simPar=simPar[1,];catchDat=NULL;variableCU=FALSE;larkPars=NULL;cuCustomCorrMat= corrmat;
  #outDir="example"; nTrials=1; makeSubDirs=TRUE; random=FALSE;larkPars=NULL;
  #erCorrMat=NULL;uniqueProd=TRUE;uniqueSurv=FALSE

  simData[[i]] <- readRDS(here("example","SamSimOutputs","simData", simPar$nameOM[i],simPar$scenario[i],
                         paste(simPar$nameOM[i],"_", simPar$nameMP[i], "_", "CUsrDat.RData",sep="")))$srDatout
  

  if(plotscn ==TRUE){
    df<-simData[[i]] 
    df<-df[df$iteration==1,]
    stackcu1<-cbind(df[,-c(9,10,11)],stack(df[,9:11]))

    p[[i]] <- ggplot(stackcu1) +
      geom_line(aes(x=year,y=values, col=as.factor(CU)))+
      geom_point(aes(x=year,y=values, col=as.factor(CU)))+
      theme_bw(14)+ theme(legend.position="none")+
      facet_wrap(~ind, scales="free_y")+
      scale_colour_viridis_d() +
      labs(title = simPar$nameOM[i])

    
  }


  nyr <- max(unique(simData[[i]]$year))

  simData[[i]] <- simData[[i]] %>% 
    filter(CU == 1, year %in% (nyr-50+1):nyr) %>% 
    filter(!is.na(obsRecruits))%>%
    mutate() %>% 
    rename(byr=year, spwn=spawners, rec=recruits, alpha_true=alpha, beta_true=beta) %>% 
    mutate(scenario = scenNames[i], alpha=99., beta=99., alpha_se=99., beta_se=99.) %>% # cols for output
    select(scenario, everything())  #reorder cols

  if(i==1){
    dlm_Out <- simData[[i]]
  }else{
    dlm_Out <-rbind(dlm_Out, simData[[i]])
  }
  
}

ggsave(
      filename = "../plots/scenarios.pdf", 
      plot = marrangeGrob(p, nrow=1, ncol=1), 
      width = 12, height = 5
    )
  

output<- readRDS(paste0("../example/SamSimOutputs/simData/", simPar$nameOM[1],"/",simPar$scenario[1],"/",
                         simPar$nameOM[1],"CUsrDat.RData",sep=""))$srDatout
  

output1<-output[output$iteration==1,]

output1<-cbind(output1[,1:3],stack(output1[,12:14]))

ggplot(output1) +
      geom_line(aes(x=year,y=values, col=as.factor(ind)))+
      geom_point(aes(x=year,y=values, col=as.factor(ind)))+
      theme_bw(14)+ 
      facet_wrap(~CU, scales="free_y")+
      scale_colour_viridis_d() +
      labs(title = simPar$nameOM[i])

estNames <- c("1_Stat", "2_Alpha_vary", "3_Beta_vary", "4_Alpha_Beta_vary")
  dlm_Out_stat <- dlm_Out_alpha <- dlm_Out_beta <- dlm_Out_alphabeta <- dlm_Out

iter <- unique(simData[[1]]$iteration)
nsc <- length(scenNames)

#check to see if models have already been run (this file is in the Google drive)
dlm_filename <- paste0("../example/estimation",   "model_estimates_all_combos.csv")

if(file.exists(dlm_filename)){
  dlm_out_all_combo <- readr::read_csv(dlm_filename)
} else{
  for(j in 1:nsc){
    for(i in seq_along(iter)){
      
      dat <- dlm_Out %>% 
        dplyr::filter(scenario == scenNames[j], iteration==i) %>% 
        dplyr::select(-c(alpha,beta,alpha_se,beta_se)) #need to remove these for fitting
      
      # alpha and beta fixed in estimation model
      dlm_model_stat <- fitDLM(data = dat,
                          alpha_vary = FALSE,
                          beta_vary = FALSE)
      
      dlm_Out_stat[which(dlm_Out_stat$scenario==scenNames[j] & dlm_Out_stat$iteration==i),] <- dlm_model_stat$results
      
      # alpha varies in estimation model
      dlm_model_alpha <- fitDLM(data = dat,
                          alpha_vary = TRUE,
                          beta_vary = FALSE)
      
      dlm_Out_alpha[which(dlm_Out_alpha$scenario==scenNames[j] & dlm_Out_alpha$iteration==i),] <- dlm_model_alpha$results
      
      # beta varies in estimation model
      dlm_model_beta <- fitDLM(data = dat,
                          alpha_vary = FALSE,
                          beta_vary = TRUE)
      
      dlm_Out_beta[which(dlm_Out_beta$scenario==scenNames[j] & dlm_Out_beta$iteration==i),] <- dlm_model_beta$results
      
      # alpha and beta vary in estimation model
      dlm_model_alphabeta <- fitDLM(data = dat,
                          alpha_vary = TRUE,
                          beta_vary = TRUE)
      
      dlm_Out_alphabeta[which(dlm_Out_alphabeta$scenario==scenNames[j] & dlm_Out_alphabeta$iteration==i),] <- dlm_model_alphabeta$results
   
    }#end j
  }# end i
}  

 # Now append the estimation model name to each dataframe
dlm_Out_stat <- dlm_Out_stat %>% 
  mutate(estModel = estNames[1])
dlm_Out_alpha <- dlm_Out_alpha %>% 
  mutate(estModel = estNames[2])
dlm_Out_beta <- dlm_Out_beta %>% 
  mutate(estModel = estNames[3])
dlm_Out_alphabeta <- dlm_Out_alphabeta %>% 
  mutate(estModel = estNames[4])

# Now make one gigantic dataframe (not sure if we want this)
dlm_out_all_combo <- rbind(dlm_Out_stat, dlm_Out_alpha, dlm_Out_beta,dlm_Out_alphabeta)

readr::write_csv(dlm_out_all_combo, here::here("example/estimation",   "model_estimates_all_combos.csv"))



dlm_out_all_combo$beta <- -dlm_out_all_combo$beta


allbias <- dlm_out_all_combo %>%
  group_by(scenario,estModel, iteration) %>%
  summarize(
    alpha_mpb=mean((alpha-alpha_true)/alpha_true)*100,
    beta_mpb=mean((beta-beta_true)/beta_true)*100) 
#allbiasstack<-cbind(allbias[,-c(4,5)],stack(allbias[,4:5]))


plotBias <- allbias %>% 
  filter(scenario %in% scenNames[1:5])

  ggplot(plotBias,aes(x=factor(scenario), y=beta_mpb, fill=estModel))+
   geom_boxplot(outlier.shape = NA)+
   #facet_wrap(vars(parameter))+
   xlab("Scenario") +
   ylab("Mean % bias")+
   ggtitle("beta")+
   theme_bw()+theme(axis.text=element_text(size=14),
                    axis.title=element_text(size=16),
                    title=element_text(size=16, face="bold"))

 ggplot(plotBias,aes(x=factor(scenario), y=alpha_mpb, fill=estModel))+
   geom_boxplot(outlier.shape = NA)+
   #facet_wrap(vars(parameter))+
   xlab("Scenario") +
   ylab("Mean % bias")+
   ggtitle("alpha")+
   theme_bw()+theme(axis.text=element_text(size=14),
                    axis.title=element_text(size=16),
                    title=element_text(size=16, face="bold"))




#==========================
#todo





# add sceatios for changes of age at spawners - later


