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
source("dlm-wrapper.R")
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


for(u in unique(simData$iteration)){
  #u=1
  dat<-simData[simData$iteration==u,]
  dat<-dat[dat$year>(max(dat$year)-46),]

  dat <- dat[!is.na(dat$obsRecruits),]
  length(unique(dat$year))
  df <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))
 

  p<-ricker_TMB(data=df)
  #p[1:5]
  pnp<-ricker_TMB(data=df, prior=0)
  #pnp[1:5]
  

  b<-ricker_stan(data=df,iter = 2000)  

  #a
   dfa<- data.frame(parameter="alpha",
    iteration=u,
    method=c(rep("MLE",2),rep("MCMC",1)),
    model=c("simple","simple_noprior","simple_bayes"),
    sim=rep(dat$alpha[nrow(dat)],3),
    est=c(p$alpha,pnp$alpha,b$alpha),
     convergence=c(p$model$convergence,
      pnp$model$convergence,
      as.numeric(abs(b$mcmcsummary["log_a","Rhat"]-1)>.1)
      ))
   dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
   
   rmsea<-aggregate((dfa$est-dfa$sim)^2, list(model=dfa$model,method=dfa$method), function(x)sqrt(mean(x)))
   rmsea$iteration<-u
   rmsea$parameter<-"alpha"
   rmsea$convergence <- c(p$model$convergence,
      pnp$model$convergence,
      as.numeric(abs(b$mcmcsummary["log_a","Rhat"]-1)>.1)
      )
  


  #b
  dfb<- data.frame(parameter="beta",
    iteration=u,
    method=c(rep("MLE",2),rep("MCMC",1)),
    model=c("simple","simple_noprior","simple_bayes"),
    sim=rep(dat$beta[nrow(dat)],3),
    est=c(p$beta, pnp$beta, b$beta),
    convergence=c(p$model$convergence,
      pnp$model$convergence,
      as.numeric(abs(b$mcmcsummary["b","Rhat"]-1)>.1))
      )
   dfb$pbias<- ((dfb$est-dfb$sim)/dfb$sim)*100

   rmseb<-aggregate((dfb$est-dfb$sim)^2, list(model=dfb$model,method=dfb$method), function(x)sqrt(mean(x)))
   rmseb$iteration<-u
   rmseb$parameter<-"beta"
   rmseb$convergence <- c(p$model$convergence,
      pnp$model$convergence,
      as.numeric(abs(b$mcmcsummary["b","Rhat"]-1)>.1)
      )
  

  #sigma
   dfsig<- data.frame(parameter="sigma",
    iteration=u,
    method=c(rep("MLE",2),rep("MCMC",1)),
    model=c("simple","simple_noprior","simple_bayes"),
    sim=rep(dat$sigma[nrow(dat)],3),
    est=c(p$sig, pnp$sig, b$sigobs),
    convergence=c(p$model$convergence,
      pnp$model$convergence,
      as.numeric(abs(b$mcmcsummary["sigma_e","Rhat"]-1)>.1)
      ))
   dfsig$pbias<- ((dfsig$est-dfsig$sim)/dfsig$sim)*100

   rmsesig<-aggregate((dfsig$est-dfsig$sim)^2, list(model=dfsig$model,method=dfsig$method), function(x)sqrt(mean(x)))
   rmsesig$iteration<-u
   rmsesig$parameter<-"sigma"
   rmsesig$convergence <- c(p$model$convergence,
      pnp$model$convergence,
      as.numeric(abs(b$mcmcsummary["sigma_e","Rhat"]-1)>.1)
      )
  
 
       
  #Smsy
  smsysim<-((1 - gsl::lambert_W0(exp(1 - dat$alpha))) /dat$beta)
 
  dfsmsy<- data.frame(parameter="smsy",
    iteration=u,
    method=c(rep("MLE",2),rep("MCMC",1)),
    model=c("simple","simple_noprior","simple_bayes"),
    sim=rep(smsysim[nrow(dat)],3),
    est=c((1 - gsl::lambert_W0(exp(1 - p$alpha))) /p$beta,
          (1 - gsl::lambert_W0(exp(1 - pnp$alpha))) /pnp$beta,
         (1 - gsl::lambert_W0(exp(1 - b$alpha))) /b$beta),
     convergence=c(p$model$convergence,
      pnp$model$convergence,
      sum(as.numeric(abs(b$mcmcsummary["log_a","Rhat"]-1)>.1),
          as.numeric(abs(b$mcmcsummary["log_b","Rhat"]-1)>.1))
      ))
  
  dfsmsy$pbias<- ((dfsmsy$est-dfsmsy$sim)/dfsmsy$sim)*100

  rmsesmsy<-aggregate((dfsmsy$est-dfsmsy$sim)^2, list(model=dfsmsy$model,method=dfsmsy$method), function(x)sqrt(mean(x)))
  rmsesmsy$iteration<-u
  rmsesmsy$parameter<-"smsy"
  rmsesmsy$convergence <- c(p$model$convergence,
      pnp$model$convergence,     
      sum(as.numeric(abs(b$mcmcsummary["log_a","Rhat"]-1)>.1),
          as.numeric(abs(b$mcmcsummary["log_b","Rhat"]-1)>.1))
      )
  


  #Sgen
  dfsgen<-data.frame(parameter="sgen",
    iteration=u,
    method=c(rep("MLE",2),rep("MCMC",1)),
    model=c("simple","simple_noprior","simple_bayes"),
    sim=unlist(mapply(sGenSolver,a=dat$alpha,Smsy=smsysim, b=dat$beta))[nrow(dat)],
    est=c(unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="simple"],Smsy=dfsmsy$est[dfsmsy$model=="simple"], b=dfb$est[dfb$model=="simple"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="simple_noprior"],Smsy=dfsmsy$est[dfsmsy$model=="simple_noprior"], b=dfb$est[dfb$model=="simple_noprior"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="simple_bayes"],Smsy=dfsmsy$est[dfsmsy$model=="simple_bayes"], b=dfb$est[dfb$model=="simple_bayes"]))
        ),
     convergence=c(p$model$convergence,
      pnp$model$convergence,
      sum(as.numeric(abs(b$mcmcsummary["log_a","Rhat"]-1)>.1),
          as.numeric(abs(b$mcmcsummary["log_b","Rhat"]-1)>.1))
      ))
   dfsgen$pbias<- ((dfsgen$est-dfsgen$sim)/dfsgen$sim)*100

   rmsesgen<-aggregate((dfsgen$est-dfsgen$sim)^2, list(model=dfsgen$model,method=dfsgen$method), function(x)sqrt(mean(x)))
   rmsesgen$iteration<-u
   rmsesgen$parameter<-"sgen"
   rmsesgen$convergence <- c(p$model$convergence,
      pnp$model$convergence,      
      sum(as.numeric(abs(b$mcmcsummary["log_a","Rhat"]-1)>.1),
          as.numeric(abs(b$mcmcsummary["log_b","Rhat"]-1)>.1))
      )
  

  #umsy
  #1-gsl::lambert_W0(exp(1 - ptva$alpha))) /ptva$beta

    dfumsy<- data.frame(parameter="umsy",
    iteration=u,
    method=c(rep("MLE",2),rep("MCMC",1)),
    model=c("simple","simple_noprior","simple_bayes"),
    sim=rep((.5 * dat$alpha - 0.07 * dat$alpha^2)[nrow(dat)],3),
    est=c((1 - gsl::lambert_W0(exp(1 - p$alpha))),
          (1 - gsl::lambert_W0(exp(1 - pnp$alpha))),
         (1 - gsl::lambert_W0(exp(1 - b$alpha)))),
    convergence=c(p$model$convergence,
      pnp$model$convergence,
      as.numeric(abs(b$mcmcsummary["log_a","Rhat"]-1)>.1)
      )
    )

    dfumsy$pbias<- ((dfumsy$est-dfumsy$sim)/dfumsy$sim)*100

    rmseumsy<-aggregate((dfumsy$est-dfumsy$sim)^2, list(model=dfumsy$model,method=dfumsy$method), function(x)sqrt(mean(x)))
    rmseumsy$iteration<-u
    rmseumsy$parameter<-"umsy"
    rmseumsy$convergence <- c(p$model$convergence,
      pnp$model$convergence,
      as.numeric(abs(b$mcmcsummary["log_a","Rhat"]-1)>.1)
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


