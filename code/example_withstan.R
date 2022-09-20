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
ricPars <- read.csv("../data/samsimIFcoho/cohoRickerpars.csv")

corrmat <- read.csv("../data/samsimIFcoho/cohoCorrMat.csv", header=F)

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simPar$species),sep = "_"))

## First check to ensure that a single scenario can be run (only a small number
# of trials necessary)
plotscn <- TRUE
p <- list()
simData <- list()





genericRecoverySim(simPar=simPar[1,], cuPar=cuPar, catchDat=NULL, srDat=NULL,
            variableCU=FALSE, ricPars=ricPars , larkPars=NULL,cuCustomCorrMat= corrmat,
            outDir="test", nTrials=100, makeSubDirs=TRUE, random=FALSE)

simData <- readRDS(here("test","SamSimOutputs","simData", simPar$nameOM[1],simPar$scenario[1],
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

  dat<-simData[simData$iteration==u,]
  dat <- dat[!is.na(dat$obsRecruits),]
  df <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))
 

  p<-rickerTMB(data=df)
  pac<-rickerTMB(data=df, AC=TRUE)
  ptva <- ricker_rwa_TMB(data=df)
  ptvb <- ricker_rwb_TMB(data=df)
  ptvab <- ricker_rwab_TMB(data=df)
  phmma <- ricker_HMM_TMB_a(data=df)
  phmmb <- ricker_HMM_TMB_b(data=df)
  phmm <- ricker_HMM_TMB(data=df)

  #a
   dfa<- data.frame(parameter="alpha",
    iteration=u,
    model=rep(c("simple","autocorr","rwa","rwb","rwab","hmma_regime","hmma_average",
      "hmmb_regime","hmmab_regime","hmmab_average"),each=nrow(df)),
    by=rep(dat$year,10),
    sim=rep(dat$alpha,10),
    est=c(rep(p$alpha,nrow(df)),rep(pac$alpha,nrow(df)),
      ptva$alpha,rep(ptvb$alpha,nrow(df)),ptvab$alpha,
      phmma$alpha[phmma$regime],c(phmma$alpha%*%phmma$probregime),
      rep(phmmb$alpha,nrow(df)),
      phmm$alpha[phmm$regime],phmm$alpha%*%phmm$probregime),
     convergence=rep(c(p$model$convergence,
      pac$model$convergence,
      ptva$model$convergence,
      ptvb$model$convergence,
      ptvab$model$convergence,
      phmma$model$convergence,
      phmma$model$convergence,
      phmmb$model$convergence,
      phmm$model$convergence,
      phmm$model$convergence
      ),each=nrow(df)))
   dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
   
   rmsea<-aggregate((dfa$est-dfa$sim)^2, list(model=dfa$model), function(x)sqrt(mean(x)))
   rmsea$iteration<-u
   rmsea$parameter<-"alpha"
   rmsea$convergence <- c(p$model$convergence,
      pac$model$convergence,
      ptva$model$convergence,
      ptvb$model$convergence,
      ptvab$model$convergence,
      phmma$model$convergence,
      phmma$model$convergence,
      phmmb$model$convergence,
      phmm$model$convergence,
      phmm$model$convergence
      )
  


  #b
   dfb<- data.frame(parameter="beta",
    iteration=u,
    model=rep(c("simple","autocorr","rwa","rwb","rwab","hmma_regime",
      "hmmb_regime","hmmb_average","hmmab_regime","hmmab_average"),each=nrow(df)),
    by=rep(dat$year,10),
    sim=rep(dat$beta,10),
    est=c(rep(p$beta,nrow(df)),rep(pac$beta,nrow(df)),
      rep(ptva$beta,nrow(df)),ptvb$beta,ptvab$beta,
      rep(phmma$beta,nrow(df)),
      phmmb$beta[phmmb$regime],phmmb$beta%*%phmmb$probregime,
      phmm$beta[phmm$regime],phmm$beta%*%phmm$probregime),
     convergence=rep(c(p$model$convergence,
      pac$model$convergence,
      ptva$model$convergence,
      ptvb$model$convergence,
      ptvab$model$convergence,
      phmma$model$convergence,
      phmmb$model$convergence,
      phmmb$model$convergence,
      phmm$model$convergence,
      phmm$model$convergence
      ),each=nrow(df)))
   dfb$pbias<- ((dfb$est-dfb$sim)/dfb$sim)*100

   rmseb<-aggregate((dfb$est-dfb$sim)^2, list(model=dfb$model), function(x)sqrt(mean(x)))
   rmseb$iteration<-u
   rmseb$parameter<-"beta"
   rmseb$convergence <- c(p$model$convergence,
      pac$model$convergence,
      ptva$model$convergence,
      ptvb$model$convergence,
      ptvab$model$convergence,
      phmma$model$convergence,
      phmmb$model$convergence,
      phmmb$model$convergence,
      phmm$model$convergence,
      phmm$model$convergence
      )
  

  #sigma

   dfsig<- data.frame(parameter="sigma",
    iteration=u,
    model=rep(c("simple","autocorr","rwa","rwb","rwab","hmma_regime",
      "hmmb_regime","hmmab_regime"),each=nrow(df)),
    by=rep(dat$year,8),
    sim=rep(dat$sigma,8),
    est=c(rep(p$sig,nrow(df)),rep(pac$sig,nrow(df)),
      rep(ptva$sig,nrow(df)),rep(ptvb$sig,nrow(df)),
      rep(ptvab$sig,nrow(df)),rep(phmma$sigma,nrow(df)),
      rep(phmmb$sigma,nrow(df)),rep(phmm$sigma,nrow(df))),
     convergence=rep(c(p$model$convergence,
      pac$model$convergence,
      ptva$model$convergence,
      ptvb$model$convergence,
      ptvab$model$convergence,
      phmma$model$convergence,
      phmmb$model$convergence,
      phmm$model$convergence
      ),each=nrow(df)))
   dfsig$pbias<- ((dfsig$est-dfsig$sim)/dfsig$sim)*100

   rmsesig<-aggregate((dfsig$est-dfsig$sim)^2, list(model=dfsig$model), function(x)sqrt(mean(x)))
   rmsesig$iteration<-u
   rmsesig$parameter<-"sigma"
   rmsesig$convergence <- c(p$model$convergence,
      pac$model$convergence,
      ptva$model$convergence,
      ptvb$model$convergence,
      ptvab$model$convergence,
      phmma$model$convergence,
      phmmb$model$convergence,
      phmm$model$convergence
      )
  
       
  #Smsy
  smsysim<-((1 - gsl::lambert_W0(exp(1 - dat$alpha))) /dat$beta)
 
  dfsmsy<- data.frame(parameter="smsy",
    iteration=u,
    model=rep(c("simple","autocorr","rwa","rwb","rwab","hmma_regime","hmma_average",
      "hmmb_regime","hmmb_average","hmmab_regime","hmmab_average"),each=nrow(df)),
    by=rep(dat$year,11),
    sim=rep(smsysim,11),
    est=c(rep((1 - gsl::lambert_W0(exp(1 - p$alpha))) /p$beta,nrow(df)),
          rep((1 - gsl::lambert_W0(exp(1 - pac$alpha))) /pac$beta,nrow(df)),
          (1 - gsl::lambert_W0(exp(1 - ptva$alpha))) /ptva$beta,
          (1 - gsl::lambert_W0(exp(1 - ptvb$alpha))) /ptvb$beta,
          (1 - gsl::lambert_W0(exp(1 - ptvab$alpha))) /ptvab$beta,
          (1 - gsl::lambert_W0(exp(1 - phmma$alpha[phmma$regime]))) /phmma$beta,
          (1 - gsl::lambert_W0(exp(1 - c(phmma$alpha%*%phmma$probregime)))) /phmma$beta,
          (1 - gsl::lambert_W0(exp(1 - phmmb$alpha))) /phmmb$beta[phmmb$regime],
          (1 - gsl::lambert_W0(exp(1 - phmmb$alpha))) /phmmb$beta%*%phmmb$probregime,
          (1 - gsl::lambert_W0(exp(1 - phmm$alpha[phmm$regime]))) /phmm$beta[phmm$regime],
          (1 - gsl::lambert_W0(exp(1 - c(phmm$alpha%*%phmm$probregime)))) /c(phmm$beta%*%phmm$probregime)),
     convergence=rep(c(p$model$convergence,
      pac$model$convergence,
      ptva$model$convergence,
      ptvb$model$convergence,
      ptvab$model$convergence,
      phmma$model$convergence,
      phmma$model$convergence,
      phmmb$model$convergence,
      phmmb$model$convergence,
      phmm$model$convergence,
      phmm$model$convergence
      ),each=nrow(df)))
  
  dfsmsy$pbias<- ((dfsmsy$est-dfsmsy$sim)/dfsmsy$sim)*100

  rmsesmsy<-aggregate((dfsmsy$est-dfsmsy$sim)^2, list(model=dfsmsy$model), function(x)sqrt(mean(x)))
  rmsesmsy$iteration<-u
  rmsesmsy$parameter<-"smsy"
  rmsesmsy$convergence <- c(p$model$convergence,
      pac$model$convergence,
      ptva$model$convergence,
      ptvb$model$convergence,
      ptvab$model$convergence,
      phmma$model$convergence,
      phmma$model$convergence,
      phmmb$model$convergence,
      phmmb$model$convergence,
      phmm$model$convergence,
      phmm$model$convergence
      )
  


  #Sgen
  dfsgen<-data.frame(parameter="sgen",
    iteration=u,
    model=rep(c("simple","autocorr","rwa","rwb","rwab","hmma_regime","hmma_average",
      "hmmb_regime","hmmb_average","hmmab_regime","hmmab_average"),each=nrow(df)),
    by=rep(dat$year,11),
    sim=rep(unlist(mapply(sGenSolver,a=dat$alpha,Smsy=smsysim, b=dat$beta)),11),
    est=c(unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="simple"],Smsy=dfsmsy$est[dfsmsy$model=="simple"], b=dfb$est[dfb$model=="simple"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="autocorr"],Smsy=dfsmsy$est[dfsmsy$model=="autocorr"], b=dfb$est[dfb$model=="autocorr"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="rwa"],Smsy=dfsmsy$est[dfsmsy$model=="rwa"], b=dfb$est[dfb$model=="rwa"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="rwb"],Smsy=dfsmsy$est[dfsmsy$model=="rwb"], b=dfb$est[dfb$model=="rwb"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="rwab"],Smsy=dfsmsy$est[dfsmsy$model=="rwab"], b=dfb$est[dfb$model=="rwab"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="hmma_regime"],Smsy=dfsmsy$est[dfsmsy$model=="hmma_regime"], b=dfb$est[dfb$model=="hmma_regime"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="hmma_average"],Smsy=dfsmsy$est[dfsmsy$model=="hmma_average"], b=dfb$est[dfb$model=="hmma_regime"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="hmmb_regime"],Smsy=dfsmsy$est[dfsmsy$model=="hmmb_regime"], b=dfb$est[dfb$model=="hmmb_regime"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="hmmb_regime"],Smsy=dfsmsy$est[dfsmsy$model=="hmmb_average"], b=dfb$est[dfb$model=="hmmb_average"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="hmmab_regime"],Smsy=dfsmsy$est[dfsmsy$model=="hmmab_regime"], b=dfb$est[dfb$model=="hmmab_regime"])),
        unlist(mapply(sGenSolver,a=dfa$est[dfa$model=="hmmab_average"],Smsy=dfsmsy$est[dfsmsy$model=="hmmab_average"], b=dfb$est[dfb$model=="hmmab_average"]))
        ),
     convergence=rep(c(p$model$convergence,
      pac$model$convergence,
      ptva$model$convergence,
      ptvb$model$convergence,
      ptvab$model$convergence,
      phmma$model$convergence,
      phmma$model$convergence,
      phmmb$model$convergence,
      phmmb$model$convergence,
      phmm$model$convergence,
      phmm$model$convergence
      ),each=nrow(df)))
   dfsgen$pbias<- ((dfsgen$est-dfsgen$sim)/dfsgen$sim)*100

   rmsesgen<-aggregate((dfsgen$est-dfsgen$sim)^2, list(model=dfsgen$model), function(x)sqrt(mean(x)))
   rmsesgen$iteration<-u
   rmsesgen$parameter<-"sgen"
   rmsesgen$convergence <- c(p$model$convergence,
      pac$model$convergence,
      ptva$model$convergence,
      ptvb$model$convergence,
      ptvab$model$convergence,
      phmma$model$convergence,
      phmma$model$convergence,
      phmmb$model$convergence,
      phmmb$model$convergence,
      phmm$model$convergence,
      phmm$model$convergence
      )
  

  #umsy
    dfumsy<- data.frame(parameter="umsy",
    iteration=u,
    model=rep(c("simple","autocorr","rwa","rwb","rwab","hmma_regime","hmma_average",
      "hmmb_regime","hmmab_regime","hmmab_average"),each=nrow(df)),
    by=rep(dat$year,10),
    sim=rep((.5 * dat$alpha - 0.07 * dat$alpha^2),10),
    est=c(rep(.5 * p$alpha - 0.07 * p$alpha^2,nrow(df)),
          rep(.5 * pac$alpha - 0.07 * pac$alpha^2,nrow(df)),
          .5 * ptva$alpha - 0.07 * ptva$alpha^2,
          rep(.5 * ptvb$alpha - 0.07 * ptvb$alpha^2,nrow(df)),
          .5 * ptvab$alpha - 0.07 * ptvab$alpha^2,
          .5 * phmma$alpha[phmma$regime] - 0.07 * phmma$alpha[phmma$regime]^2,
          .5 * c(phmma$alpha%*%phmma$probregime) - 0.07 * c(phmma$alpha%*%phmma$probregime)^2,
          rep(.5 * phmmb$alpha - 0.07 * phmmb$alpha^2,nrow(df)),
          .5 * phmm$alpha[phmm$regime] - 0.07 * phmm$alpha[phmm$regime]^2,
          .5 * c(phmm$alpha%*%phmm$probregime) - 0.07 * c(phmm$alpha%*%phmm$probregime)^2),
     convergence=rep(c(p$model$convergence,
      pac$model$convergence,
      ptva$model$convergence,
      ptvb$model$convergence,
      ptvab$model$convergence,
      phmma$model$convergence,
      phmma$model$convergence,
      phmmb$model$convergence,
      phmm$model$convergence,
      phmm$model$convergence
      ),each=nrow(df))
    )

    dfumsy$pbias<- ((dfumsy$est-dfumsy$sim)/dfumsy$sim)*100

    rmseumsy<-aggregate((dfumsy$est-dfumsy$sim)^2, list(model=dfumsy$model), function(x)sqrt(mean(x)))
    rmseumsy$iteration<-u
    rmseumsy$parameter<-"umsy"
    rmseumsy$convergence <- c(p$model$convergence,
      pac$model$convergence,
      ptva$model$convergence,
      ptvb$model$convergence,
      ptvab$model$convergence,
      phmma$model$convergence,
      phmma$model$convergence,
      phmmb$model$convergence,
      phmm$model$convergence,
      phmm$model$convergence
      )
  

   rmse[[u]]<-rbind(rmsea,rmseb,rmsesig,rmsesmsy,rmsesgen,rmseumsy)
   simest[[u]]<-rbind(dfa,dfb,dfsig,dfsmsy,dfsgen,dfumsy)

   


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

  
#=================================
#plots

dfpbias<-do.call("rbind",simest)
dfpbias<- dfpbias[dfpbias$convergence==0,]
dfpbias$model <- factor(dfpbias$model, levels=c("simple", "autocorr", "rwa", "rwb",
           "rwab", "hmma_regime", "hmma_average", "hmmb_regime", "hmmb_average", "hmmab_regime", "hmmab_average"))

summary(dfpbias)

fig <- ggplot(dfpbias,aes(x=model,y=pbias)) +
geom_boxplot() +
coord_cartesian(ylim = c(-100,100))+
geom_hline(yintercept=0) +
theme_bw(14)+ theme(legend.position="none")+
facet_wrap(~parameter, scales="free_y")+
scale_colour_viridis_d() +
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


