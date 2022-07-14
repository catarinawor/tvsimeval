#=============================================================
#Example use of samsim for time varying simulation evaluation
#using the coho data as it is compliant with the most recent
#samSim updates
#Catarina Wor
# March 2022 
#=============================================================

#TODO
#add true Sgen, smsy and umsy to data output
#polta of param scenarios
#plot of Er scenarios
#add estimation routines
 
if(!"KFfuncs" %in% rownames(installed.packages())){
  remotes::install_github("carrieholt/KF-funcs")
}

#install samsim 
#remotes::install_github("Pacific-salmon-assess/samSim", force=TRUE, ref="timevar")

library(samSim)
library(dplyr)
library(KFfuncs)
library(cmdstanr);library(loo)


source("dlm-wrapper.R")

compile("TMB/Ricker_tva_Smax_ratiovar.cpp")
dyn.load(dynlib("TMB/Ricker_tva_Smax_ratiovar"))


system(paste("cp stan/ricker_linear_varying_a.stan", paste0(cmdstan_path(),"/timevarmodels")))
system(paste("cp stan/ricker_linear_varying_b.stan", paste0(cmdstan_path(),"/timevarmodels")))
system(paste("cp stan/ricker_linear_varying_a_and_b.stan", paste0(cmdstan_path(),"/timevarmodels")))
system(paste("cp stan/ricker_linear_varying_a_GP.stan", paste0(cmdstan_path(),"/timevarmodels")))
#system(paste("cp stan/ricker_linear_varying_b_GP.stan", paste0(cmdstan_path(),"/stanmodels")))
#system(paste("cp stan/ricker_linear_varying_a_and_b_GP.stan", paste0(cmdstan_path(),"/stanmodels")))



#' simest
#' 
#' Run samSim and estimation models
#' 
#' 
#'@param simPar samSim input file, see ?samSim::simParexam ple for details
#'@param cuPar samSim input file
#'@param srDat samSim input file, stock recruit initial data to initilaize simulation
#'@param ricPars samSim input file, ricker stock recruitment parameters for simulation
#'@param corrmat samSim input file, correlation matrix among CUs
#'@param outDir head directory for output simulations
#'@param iteration number of simulation runs, must be >1
#' 
#' 
#' 
#'@return A list containing two objects:
#' * 
#' * simData simulation true and observed data 
#' * estData estimation results   
#' 
#' 
sim_est <- function(simPar, cuPar, srDat, ricPars, corrmat, outDir, iteration, simrun=TRUE){

  ## Store relevant object names to help run simulation 
  scenNames <- unique(simPar$scenario)
  dirNames <- sapply(scenNames, function(x) paste(x, unique(simPar$species),sep = "_"))

  plotscn <- TRUE
  p <- list()
  simData <- list()
  estData <- list()

  for(i in seq_len(nrow(simPar))){

    #i=1;iteration=2
    if(simrun){
      genericRecoverySim(simPar=simPar[i,], cuPar=cuPar, catchDat=NULL, srDat=srDat,
              variableCU=FALSE, ricPars=ricPars , larkPars=NULL,cuCustomCorrMat= corrmat,
              outDir=outDir, nTrials=iteration, makeSubDirs=TRUE, random=FALSE)
      #simPar=simPar[1,];catchDat=NULL;variableCU=FALSE;larkPars=NULL;cuCustomCorrMat= corrmat;
      #outDir="example"; nTrials=1; makeSubDirs=TRUE; random=FALSE;larkPars=NULL;
      #erCorrMat=NULL;uniqueProd=TRUE;uniqueSurv=FALSE
    }
    
    simData[[i]] <- readRDS(here(outDir,"SamSimOutputs","simData", simPar$nameOM[i],simPar$scenario[i],
                         paste(simPar$nameOM[i],"_", simPar$nameMP[i], "_", "CUsrDat.RData",sep="")))$srDatout
  
    nyr <- max(unique(simData[[i]]$year))

    simData[[i]] <- simData[[i]] %>% 
      filter(CU == 1, year %in% (nyr-simPar[i,"simYears"]+1):nyr) %>% 
      filter(!is.na(obsRecruits))%>%
      mutate() %>% 
      rename(byr=year, spwn=spawners, rec=recruits, alpha_true=alpha, beta_true=beta) %>% 
      mutate(scenario = scenNames[i]) %>% # cols for output
      select(scenario, everything()) 

    nobs<-length(unique(simData[[i]]$byr))

    #dlm
    dlmalpha_alpha <- matrix(NA,ncol=iteration,nrow=nobs)
    dlmalpha_beta <- matrix(NA,ncol=iteration,nrow=nobs)
    dlmalpha_sigma <- matrix(NA,ncol=iteration,nrow=nobs)
    dlmalpha_smsy <- matrix(NA,ncol=iteration,nrow=nobs)
    dlmalpha_umsy <- matrix(NA,ncol=iteration,nrow=nobs)
    dlmalpha_sgen <- matrix(NA,ncol=iteration,nrow=nobs)
    dlmalpha_AIC <- rep(NA,length=iteration)
    dlmalpha_BIC <- rep(NA,length=iteration)
    dlmalpha_nll <- rep(NA,length=iteration)
    dlmalpha_convergence <- rep(NA,length=iteration)
    
    dlmbeta_alpha <- matrix(NA,ncol=iteration,nrow=nobs)
    dlmbeta_beta <- matrix(NA,ncol=iteration,nrow=nobs)
    dlmbeta_sigma <- matrix(NA,ncol=iteration,nrow=nobs)
    dlmbeta_smsy <- matrix(NA,ncol=iteration,nrow=nobs)
    dlmbeta_umsy <- matrix(NA,ncol=iteration,nrow=nobs)
    dlmbeta_sgen <- matrix(NA,ncol=iteration,nrow=nobs)
    dlmbeta_AIC <- rep(NA,length=iteration)
    dlmbeta_BIC <- rep(NA,length=iteration)
    dlmbeta_nll <- rep(NA,length=iteration)
    dlmbeta_convergence <- rep(NA,length=iteration)
    
    dlmalphabeta_alpha <- matrix(NA,ncol=iteration,nrow=nobs)
    dlmalphabeta_beta <- matrix(NA,ncol=iteration,nrow=nobs)
    dlmalphabeta_sigma <- matrix(NA,ncol=iteration,nrow=nobs)
    dlmalphabeta_smsy <- matrix(NA,ncol=iteration,nrow=nobs)
    dlmalphabeta_umsy <- matrix(NA,ncol=iteration,nrow=nobs)
    dlmalphabeta_sgen <- matrix(NA,ncol=iteration,nrow=nobs)
    dlmalphabeta_AIC <- rep(NA,length=iteration)
    dlmalphabeta_BIC <- rep(NA,length=iteration)
    dlmalphabeta_nll <- rep(NA,length=iteration)
    dlmalphabeta_convergence <- rep(NA,length=iteration)
    

    #KF 
    tmbKF_alpha <- matrix(NA,ncol=iteration,nrow=nobs)
    tmbKF_beta <- matrix(NA,ncol=iteration,nrow=nobs)
    tmbKF_sigma <- matrix(NA,ncol=iteration,nrow=nobs)
    tmbKF_smsy <- matrix(NA,ncol=iteration,nrow=nobs)
    tmbKF_umsy <- matrix(NA,ncol=iteration,nrow=nobs)
    tmbKF_sgen <- matrix(NA,ncol=iteration,nrow=nobs)

    tmbKF_AIC <- rep(NA,length=iteration)
    tmbKF_BIC <- rep(NA,length=iteration)
    tmbKF_nll <- rep(NA,length=iteration)
    tmbKF_convergence <- rep(NA,length=iteration)

    #Recursive Bayes tmb  
    tmbRB_alpha <- matrix(NA,ncol=iteration,nrow=nobs)
    tmbRB_beta <- matrix(NA,ncol=iteration,nrow=nobs)
    tmbRB_sigma <- matrix(NA,ncol=iteration,nrow=nobs)
    tmbRB_smsy <- matrix(NA,ncol=iteration,nrow=nobs)
    tmbRB_umsy <- matrix(NA,ncol=iteration,nrow=nobs)
    tmbRB_sgen <- matrix(NA,ncol=iteration,nrow=nobs)

    tmbRB_AIC <- rep(NA,length=iteration)
    tmbRB_BIC <- rep(NA,length=iteration)
    tmbRB_nll <- rep(NA,length=iteration)
    tmbRB_convergence <- rep(NA,length=iteration)

    #Recursive Bayes alpha vary  stan
    stanRBalpha_alpha <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBalpha_beta <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBalpha_sigma <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBalpha_smsy <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBalpha_umsy <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBalpha_sgen <- matrix(NA,ncol=iteration,nrow=nobs)

    stanRBalpha_loo <- list()
    #stanRBalpha_AIC <- rep(NA,length=iteration)
    #stanRBalpha_BIC <- rep(NA,length=iteration)
    #stanRBalpha_nll <- rep(NA,length=iteration)
    
    stanRBalpha_convergence_alpha <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBalpha_convergence_beta <- rep(NA,length=iteration)
    stanRBalpha_convergence_sigma <- rep(NA,length=iteration)
    
    #Recursive Bayes beta vary  stan
    stanRBbeta_alpha <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBbeta_beta <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBbeta_sigma <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBbeta_smsy <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBbeta_umsy <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBbeta_sgen <- matrix(NA,ncol=iteration,nrow=nobs)
    
    stanRBbeta_loo <- list()
    
    stanRBbeta_convergence_alpha <- rep(NA,length=iteration)
    stanRBbeta_convergence_beta <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBbeta_convergence_sigma <- rep(NA,length=iteration)

    #Recursive Bayes alpha and beta vary  stan
    stanRBalphabeta_alpha <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBalphabeta_beta <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBalphabeta_sigma <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBalphabeta_smsy <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBalphabeta_umsy <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBalphabeta_sgen <- matrix(NA,ncol=iteration,nrow=nobs)

    stanRBalphabeta_loo <- list()
    
    stanRBalphabeta_convergence_alpha <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBalphabeta_convergence_beta <- matrix(NA,ncol=iteration,nrow=nobs)
    stanRBalphabeta_convergence_sigma <- rep(NA,length=iteration)

    #Recursive Bayes alpha and beta vary  stan
    stanGPalpha_alpha <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPalpha_beta <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPalpha_sigma <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPalpha_smsy <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPalpha_umsy <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPalpha_sgen <- matrix(NA,ncol=iteration,nrow=nobs)

    stanGPalpha_loo <- list()
    
    #beta gaussian
    stanGPbeta_convergence_alpha <- rep(NA,length=iteration)
    stanGPbeta_convergence_beta <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPbeta_convergence_sigma <- rep(NA,length=iteration)

    stanGPbeta_alpha <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPbeta_beta <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPbeta_sigma <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPbeta_smsy <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPbeta_umsy <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPbeta_sgen <- matrix(NA,ncol=iteration,nrow=nobs)

    stanGPbeta_loo <- list()
    
    stanGPbeta_convergence_alpha <- rep(NA,length=iteration)
    stanGPbeta_convergence_beta <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPbeta_convergence_sigma <- rep(NA,length=iteration)

    #alpha beta gaussian
    stanGPalphabeta_convergence_alpha <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPalphabeta_convergence_beta <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPalphabeta_convergence_sigma <- rep(NA,length=iteration)

    stanGPalphabeta_alpha <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPalphabeta_beta <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPalphabeta_sigma <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPalphabeta_smsy <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPalphabeta_umsy <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPalphabeta_sgen <- matrix(NA,ncol=iteration,nrow=nobs)

    stanGPalphabeta_loo <- list()
    
    stanGPalphabeta_convergence_alpha <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPalphabeta_convergence_beta <- matrix(NA,ncol=iteration,nrow=nobs)
    stanGPalphabeta_convergence_sigma <- rep(NA,length=iteration)

    #n=1
    for( n in seq_len(iteration)){
      print(paste("iteration",n))

      dat <- simData[[i]]  %>% 
        dplyr::filter( iteration==n, CU==1) %>% 
        dplyr::select(byr, obsSpawners,obsRecruits ) %>%
        dplyr::rename(spwn=obsSpawners, rec=obsRecruits)
 
      tryCatch({
        outdlm_alpha <- fitDLM(data = dat,
                            alpha_vary = TRUE,
                            beta_vary = FALSE)

        dlmalpha_alpha[,n] <- outdlm_alpha$results$alpha
        dlmalpha_beta[,n] <- -outdlm_alpha$results$beta
        dlmalpha_sigma[,n] <- rep(outdlm_alpha$sd.est[1],nrow(dat))
        dlmalpha_smsy[,n] <- outdlm_alpha$results$Smsy
        dlmalpha_umsy[,n] <- outdlm_alpha$results$Umsy
        dlmalpha_sgen[,n] <- outdlm_alpha$results$Sgen

        dlmalpha_AIC[n] <- outdlm_alpha$AICc
        dlmalpha_BIC[n] <- outdlm_alpha$BIC
        dlmalpha_nll[n] <- outdlm_alpha$nll
        dlmalpha_convergence[n] <- outdlm_alpha$convergence

      },  error = function(e) {
        print(paste("dlm_alpha failed in iteration",n))
      })
  
      tryCatch({
        # beta varies in estimation model
        outdlm_beta <- fitDLM(data = dat,
                            alpha_vary = FALSE,
                            beta_vary = TRUE)

        dlmbeta_alpha[,n] <- outdlm_beta$results$alpha
        dlmbeta_beta[,n] <- -outdlm_beta$results$beta
        dlmbeta_sigma[,n] <- rep(outdlm_beta$sd.est[1],nrow(dat))
        dlmbeta_smsy[,n] <- outdlm_beta$results$Smsy
        dlmbeta_umsy[,n] <- outdlm_beta$results$Umsy
        dlmbeta_sgen[,n] <- outdlm_beta$results$Sgen

        dlmbeta_AIC[n] <- outdlm_beta$AICc
        dlmbeta_BIC[n] <- outdlm_beta$BIC
        dlmbeta_nll[n] <- outdlm_beta$nll
        dlmbeta_convergence[n] <- outdlm_beta$convergence


      },  error = function(e) {
        print(paste("dlm_beta failed in iteration",n))
      })
    
      tryCatch({
        # alpha and beta vary in estimation model
        outdlm_alphabeta <- fitDLM(data = dat,
                            alpha_vary = TRUE,
                            beta_vary = TRUE)
        
        dlmalphabeta_alpha[,n] <- outdlm_alphabeta$results$alpha
        dlmalphabeta_beta[,n] <- -outdlm_alphabeta$results$beta
        dlmalphabeta_sigma[,n] <- rep(outdlm_alphabeta$sd.est[1],nrow(dat))
        dlmalphabeta_smsy[,n] <- outdlm_alphabeta$results$Smsy
        dlmalphabeta_umsy[,n] <- outdlm_alphabeta$results$Umsy
        dlmalphabeta_sgen[,n] <- outdlm_alphabeta$results$Sgen

        dlmalphabeta_AIC[n] <- outdlm_alphabeta$AICc
        dlmalphabeta_BIC[n] <- outdlm_alphabeta$BIC
        dlmalphabeta_nll[n] <- outdlm_alphabeta$nll
        dlmalphabeta_convergence[n] <- outdlm_alphabeta$convergence

      },  error = function(e) {
        print(paste("dlm_alphabeta failed in iteration",n))
      })

      #add
      #Kalman Filter Holt TMB - avary
      tryCatch({
        # alpha and beta vary in estimation model
        s <- data.frame(S=dat$spwn,
          R=dat$rec)
        rekf <- kfTMB(data=s, silent = FALSE, control = TMBcontrol())
        kfrep <- summary(sdreport(rekf$tmb_obj))
        
        tmbKF_alpha[,n] <- kfrep[which(rownames(kfrep)=="smoothemeana"),1] 
        tmbKF_beta[,n] <- -rep(kfrep[which(rownames(kfrep)=="b"),1], nrow(dat))
        tmbKF_sigma[,n] <- rep(kfrep[which(rownames(kfrep)=="sige"),1], nrow(dat))
        tmbKF_smsy[,n] <- (1 - gsl::lambert_W0(exp(1 - tmbKF_alpha[,n])))/tmbKF_beta[,n]
        tmbKF_umsy[,n] <- .5 * tmbKF_alpha[,n] - 0.07 * (tmbKF_alpha[,n])^2
        tmbKF_sgen[,n] <- unlist(mapply(sGenSolverdlm, a=tmbKF_alpha[,n],Smsy=tmbKF_smsy[,n], b=tmbKF_beta[,n]))

        tmbKF_nll[n]<-rekf$model$objective
        nparskftmb<-length(rekf$model$par)

        tmbKF_AIC[n] <- 2*tmbKF_nll[n] + 2*nparskftmb +(2*nparskftmb*(nparskftmb+1)/(nrow(dat)-nparskftmb-1))
        tmbKF_BIC[n] <- nparskftmb*log(nrow(dat)) +2*tmbKF_nll[n]
        tmbKF_convergence[n] <- rekf$model$convergence


      },  error = function(e) {
        print(paste("KF TMB failed in iteration",n))
      })

      #Recursive Bayes TMB - avary
      tryCatch({

        srm <- lm(log(dat$rec/dat$spwn) ~ dat$spwn)
        SRdata<-list(obs_logRS=log(dat$rec/dat$spwn),obs_S=dat$spwn, prbeta1=1.5,
        prbeta2=1.5)
  
        #Model 1 - TMB
        parameters<- list(
          alphao=srm$coefficients[[1]],
          logSmax = log(1/ifelse(-srm$coefficients[[2]]<0,1e-08,-srm$coefficients[2])),
          rho=.5,
          logvarphi=0,
          alpha=rep(srm$coefficients[1],length(dat$rec))
        )    

        obj <- MakeADFun(SRdata,parameters,DLL="Ricker_tva_Smax_ratiovar",random="alpha")#,lower = -Inf, upper = Inf)
        newtonOption(obj, smartsearch=FALSE)

        opt <- nlminb(obj$par,obj$fn,obj$gr)
  
        sdrep <- summary(sdreport(obj))

        Smsyrb <- (1 - gsl::lambert_W0(exp(1 - sdrep[rownames(sdrep)=="alpha",1]))) /sdrep[rownames(sdrep)=="beta",1]
    
        
        tmbRB_alpha[,n] <- sdrep[which(rownames(sdrep)=="alpha"),1] 
        tmbRB_beta[,n] <- rep(sdrep[which(rownames(sdrep)=="beta"),1], nrow(dat))
        tmbRB_sigma[,n] <- rep(sdrep[which(rownames(sdrep)=="sig"),1], nrow(dat))
        tmbRB_smsy[,n] <- (1 - gsl::lambert_W0(exp(1 - tmbRB_alpha[,n])))/tmbRB_beta[,n]
        tmbRB_umsy[,n] <- sdrep[which(rownames(sdrep)=="umsy"),1] 
        tmbRB_sgen[,n] <- unlist(mapply(sGenSolverdlm, a=tmbRB_alpha[,n],Smsy=tmbRB_smsy[,n], b=tmbRB_beta[,n]))

        tmbRB_nll[n]<-opt$objective
        nparsrbtmb<-length(opt$par)

        tmbRB_AIC[n] <- 2*tmbRB_nll[n] + 2*nparsrbtmb +(2*nparsrbtmb*(nparsrbtmb+1)/(nrow(dat)-nparsrbtmb-1))
        tmbRB_BIC[n] <- nparsrbtmb*log(nrow(dat)) +2*tmbRB_nll[n]
        tmbRB_convergence[n] <- opt$convergence


      },  error = function(e) {
        print(paste("KF TMB failed in iteration",n))
      })

      #Recursive Bayes stan - a,b and abvary
      tryCatch({

        file1 <- file.path(cmdstan_path(),'timevarmodels',"ricker_linear_varying_a.stan")

        mod <- cmdstan_model(file1)
        #mcmc
        fit<- mod$sample(
            data =list(R_S = log(dat$rec/dat$spwn),
                                                 N=nrow(dat),
                                                 TT=as.numeric(factor(seq_len(nrow(dat)))),
                                                 S=c(dat$spwn)),
            seed = 123, 
            chains = 6, 
            parallel_chains = 6,
            iter_warmup = 500,
            iter_sampling = 1000,
            refresh = 500,
            adapt_delta = 0.99,
            max_treedepth = 20 # print update every 500 iters
          )

        stanRBalpha_loo[[n]] <- fit$loo(cores = 2)
        
        params1<- fit$draws(format='df',variables=c('log_a','b','log_b','sigma_a','sigma_e'))
        parssummary<-summary(params1)
        
       
        stanRBalpha_alpha[,n] <- parssummary[grep("log_a",parssummary$variable),]$median
        stanRBalpha_beta[,n] <- rep(parssummary[parssummary$variable=="b",]$median, nrow(dat))
        stanRBalpha_sigma[,n] <- rep(parssummary[parssummary$variable=="sigma_e",]$median, nrow(dat))

        #params_stan_rb<-rstan::extract(stan_rb)
        
        logastanrb <- params1[,grep("log_a",names(params1))]
        umsystanrb <- .5 *  logastanrb - 0.07 * ( logastanrb)^2
        Smsystanrb <- matrix(NA, nrow=nrow(umsystanrb), ncol=ncol(umsystanrb))
        Sgenstanrb <- matrix(NA, nrow=nrow(umsystanrb), ncol=ncol(umsystanrb))
        
        for(j in 1:ncol(Smsystanrb)){
          Smsystanrb[,j] <- (1 - gsl::lambert_W0(exp(1 - logastanrb[[j]] )))/params1$b
          Sgenstanrb[,j] <- unlist(mapply(sGenSolverdlm,a=logastanrb[[j]],
             Smsy=Smsystanrb[,j], b=params1$b))
        }
    
        
        stanRBalpha_smsy[,n] <- apply(Smsystanrb,2,median)
        stanRBalpha_umsy[,n] <- apply(umsystanrb,2,median) 
        stanRBalpha_sgen[,n] <- apply(Sgenstanrb,2,median) 
       
        
        #remove mle - these seem t0 change depending on seed 
        #also lack of convergence is frequent
        #fit_mle <- mod$optimize(data = list(R_S = log(dat$rec/dat$spwn),
        #                                         N=nrow(dat),
        #                                         TT=as.numeric(factor(seq_len(nrow(dat)))),
        #                                         S=c(dat$spwn)),
        #                                          seed = 12234) 
        
        #s <- fit_mle$return_codes()

        #fit_mle$summary()$variable
        #fit_mle$lp()
        #nllstanRBalpha <-fit_mle$mle()
        #nparsstanRBalpha <-length(opt$par)
        
        #stanRBalpha_AIC[n] <- 2*nllrbtmb + 2*nparsrbtmb +(2*nparsrbtmb*(nparsrbtmb+1)/(nrow(dat)-nparsrbtmb-1))
        #stanRBalpha_BIC[n] <- nparsrbtmb*log(nrow(dat)) +2*nllrbtmb
        stanRBalpha_convergence_alpha[,n] <- as.numeric(abs(parssummary[grep("log_a",parssummary$variable),]$rhat-1)>.1)
        stanRBalpha_convergence_beta[n] <- as.numeric(abs(parssummary[parssummary$variable=="b",]$rhat-1)>.1)
        stanRBalpha_convergence_sigma[n] <- as.numeric(abs(parssummary[parssummary$variable=="sigma_e",]$rhat-1)>.1)

        

      },  error = function(e) {
        print(paste("stan RBa failed in iteration",n))
      })
      
      tryCatch({

        file2 <- file.path(cmdstan_path(),'timevarmodels',"ricker_linear_varying_b.stan")

        mod2 <- cmdstan_model(file2)
        #mcmc
        fit2<- mod2$sample(
            data =list(R_S = log(dat$rec/dat$spwn),
                                                 N=nrow(dat),
                                                 TT=as.numeric(factor(seq_len(nrow(dat)))),
                                                 S=c(dat$spwn)),
            seed = 123, 
            chains = 6, 
            parallel_chains = 6,
            iter_warmup = 500,
            iter_sampling = 1000,
            refresh = 500,
            adapt_delta = 0.99,
            max_treedepth = 20 # print update every 500 iters
          )

        stanRBbeta_loo[[n]] <- fit2$loo(cores = 2)
        
        params2<- fit2$draws(format='df',variables=c('log_a','b','log_b','sigma_b','sigma_e'))
        parssummary2<-summary(params2)
        
        stanRBbeta_alpha[,n] <- rep(parssummary2[grep("log_a",parssummary2$variable),]$median, nrow(dat))
        stanRBbeta_beta[,n] <- parssummary2[grep("^b",parssummary2$variable),]$median
        stanRBbeta_sigma[,n] <- rep(parssummary2[parssummary2$variable=="sigma_e",]$median, nrow(dat))

        #params_stan_rb<-rstan::extract(stan_rb)
        stanRBbeta_umsy[,n] <- rep(median(.5 *  params2$log_a - 0.07 * ( params2$log_a)^2), nrow(dat))
        
        logastanrb <- params2[,grep("^b",names(params2))]
        Smsystanrb <- matrix(NA, nrow=nrow(logastanrb), ncol=ncol(logastanrb))
        Sgenstanrb <- matrix(NA, nrow=nrow(logastanrb), ncol=ncol(logastanrb))
        
        for(j in 1:ncol(Smsystanrb)){
          Smsystanrb[,j] <- (1 - gsl::lambert_W0(exp(1 - params2$log_a )))/logastanrb[[j]]
          Sgenstanrb[,j] <- unlist(mapply(sGenSolverdlm,a=params2$log_a,
             Smsy=Smsystanrb[,j], b=logastanrb[[j]]))
        }
    
        
        stanRBbeta_smsy[,n] <- apply(Smsystanrb,2,median)
        stanRBbeta_sgen[,n] <- apply(Sgenstanrb,2,median)


        stanRBbeta_convergence_alpha[n] <- as.numeric(abs(parssummary2[parssummary2$variable=="b",]$rhat-1)>.1)        
        stanRBbeta_convergence_beta[,n] <- as.numeric(abs(parssummary2[grep("^b",parssummary2$variable),]$rhat-1)>.1)
        stanRBbeta_convergence_sigma[n] <- as.numeric(abs(parssummary2[parssummary2$variable=="sigma_e",]$rhat-1)>.1)
   
        
      },  error = function(e) {
        print(paste("stan RB b failed in iteration",n))
      })


      tryCatch({

        file3 <- file.path(cmdstan_path(),'timevarmodels',"ricker_linear_varying_a_and_b.stan")

        mod3 <- cmdstan_model(file3)
        #mcmc
        fit3<- mod3$sample(
            data =list(R_S = log(dat$rec/dat$spwn),
                                                 N=nrow(dat),
                                                 TT=as.numeric(factor(seq_len(nrow(dat)))),
                                                 S=c(dat$spwn)),
            seed = 123, 
            chains = 6, 
            parallel_chains = 6,
            iter_warmup = 500,
            iter_sampling = 2000,
            refresh = 500,
            adapt_delta = 0.99,
            max_treedepth = 20 # print update every 500 iters
          )

        stanRBalphabeta_loo[[n]] <- fit3$loo(cores = 2)
        
        params3<- fit3$draws(format='df',variables=c('log_a','b','log_b','sigma_b','sigma_e'))
        parssummary3<-summary(params3)
        
        stanRBalphabeta_alpha[,n] <- parssummary3[grep("log_a",parssummary3$variable),]$median
        stanRBalphabeta_beta[,n] <- parssummary3[grep("^b",parssummary3$variable),]$median
        stanRBalphabeta_sigma[,n] <- rep(parssummary3[parssummary3$variable=="sigma_e",]$median, nrow(dat))

        #params_stan_rb<-rstan::extract(stan_rb)
        
        logastanrba <- params3[,grep("log_a",names(params3))]
        logastanrbb <- params3[,grep("^b",names(params3))]
        umsystanrb <- .5 *  logastanrba - 0.07 * ( logastanrba)^2
         
        Smsystanrb <- matrix(NA, nrow=nrow(logastanrba), ncol=ncol(logastanrba))
        Sgenstanrb <- matrix(NA, nrow=nrow(logastanrba), ncol=ncol(logastanrba))
        
        for(j in 1:ncol(Smsystanrb)){
          Smsystanrb[,j] <- (1 - gsl::lambert_W0(exp(1 - logastanrba[[j]])))/logastanrbb[[j]]
          Sgenstanrb[,j] <- unlist(mapply(sGenSolverdlm,a=logastanrba[[j]],
             Smsy=Smsystanrb[,j], b=logastanrbb[[j]]))
        }
    
        stanRBalphabeta_umsy[,n] <- apply(umsystanrb,2,median)
        stanRBalphabeta_smsy[,n] <- apply(Smsystanrb,2,median)
        stanRBalphabeta_sgen[,n] <- apply(Sgenstanrb,2,median)


        stanRBalphabeta_convergence_alpha[,n] <- as.numeric(abs(parssummary3[grep("log_a",parssummary3$variable),]$rhat-1)>.1)       
        stanRBalphabeta_convergence_beta[,n] <- as.numeric(abs(parssummary3[grep("^b",parssummary3$variable),]$rhat-1)>.1)
        stanRBalphabeta_convergence_sigma[n] <- as.numeric(abs(parssummary3[parssummary3$variable=="sigma_e",]$rhat-1)>.1)
   
        
      },  error = function(e) {
        print(paste("stan RB b failed in iteration",n))
      })

      #Gaussian process stan - a, b and and b vary
      tryCatch({

        fileGP1 <- file.path(cmdstan_path(),'timevarmodels',"ricker_linear_varying_a_GP.stan")

        modGP1 <- cmdstan_model(fileGP1)
        #mcmc
        fitGP1<- modGP1$sample(
            data =list(N=nrow(dat),
                      TT=as.numeric(factor(seq_len(nrow(dat)))),
                      R_S = log(dat$rec/dat$spwn),
                      S=c(dat$spwn)),
            seed = 123, 
            chains = 6, 
            parallel_chains = 6,
            iter_warmup = 500,
            iter_sampling = 2000,
            refresh = 500,
            adapt_delta = 0.99,
            max_treedepth = 20 # print update every 500 iters
          )
        
        stanGPalpha_loo[[n]] <- fitGP1$loo(cores = 2)
        paramsGP1<- fitGP1$draws(format='df',variables=c('log_a','b','log_b','sigma_e'))
        parssummaryGP1<-summary(paramsGP1)
        
        stanGPalpha_alpha[,n] <- parssummaryGP1[grep("log_a",parssummaryGP1$variable),]$median
        stanGPalpha_beta[,n] <- rep(parssummaryGP1[grep("^b",parssummaryGP1$variable),]$median, nrow(dat))
        stanGPalpha_sigma[,n] <- rep(parssummaryGP1[parssummaryGP1$variable=="sigma_e",]$median, nrow(dat))

        #params_stan_rb<-rstan::extract(stan_rb)
        
        logastangpa <- paramsGP1[,grep("log_a",names(paramsGP1))]
        umsystangp <- .5 *  logastangpa - 0.07 * ( logastangpa)^2
         
        Smsystangp <- matrix(NA, nrow=nrow(logastangpa), ncol=ncol(logastangpa))
        Sgenstangp <- matrix(NA, nrow=nrow(logastangpa), ncol=ncol(logastangpa))
        
        for(j in 1:ncol(Smsystangp)){
          Smsystangp[,j] <- (1 - gsl::lambert_W0(exp(1 - logastangpa[[j]] )))/paramsGP1$b
          Sgenstangp[,j] <- unlist(mapply(sGenSolverdlm,a=logastangpa[[j]],
             Smsy=Smsystangp[,j], b=paramsGP1$b))
        }
    
        stanGPalpha_umsy[,n] <- apply(umsystangp,2,median)
        stanGPalpha_smsy[,n] <- apply(Smsystangp,2,median)
        stanGPalpha_sgen[,n] <- apply(Sgenstangp,2,median)

      
        stanGPalpha_convergence_alpha[,n] <- as.numeric(abs(parssummaryGP1[grep("log_a",parssummaryGP1$variable),]$rhat-1)>.1)       
        stanGPalpha_convergence_beta[n] <- as.numeric(abs(parssummaryGP1[grep("^b",parssummaryGP1$variable),]$rhat-1)>.1)
        stanGPalpha_convergence_sigma[n] <- as.numeric(abs(parssummaryGP1[parssummaryGP1$variable=="sigma_e",]$rhat-1)>.1)
   
        
      },  error = function(e) {
        print(paste("stan RB b failed in iteration",n))
      })

      tryCatch({

        fileGP2 <- file.path(cmdstan_path(),'timevarmodels',"ricker_linear_varying_b_GP.stan")

        modGP2 <- cmdstan_model(fileGP2)
        #mcmc
        fitGP2<- modGP2$sample(
            data =list(N=nrow(dat),
                      TT=as.numeric(factor(seq_len(nrow(dat)))),
                      R_S = log(dat$rec/dat$spwn),
                      S=c(dat$spwn)),
            seed = 123, 
            chains = 6, 
            parallel_chains = 6,
            iter_warmup = 500,
            iter_sampling = 2000,
            refresh = 500,
            adapt_delta = 0.99,
            max_treedepth = 20 # print update every 500 iters
          )
        
        stanGPbeta_loo[[n]] <- fitGP2$loo(cores = 2)
        paramsGP2<- fitGP2$draws(format='df',variables=c('log_a','b','log_b','sigma_e'))
        parssummaryGP2<-summary(paramsGP2)
        
        stanGPbeta_alpha[,n] <- parssummaryGP2[grep("log_a",parssummaryGP2$variable),]$median
        stanGPbeta_beta[,n] <- rep(parssummaryGP2[grep("^b",parssummaryGP2$variable),]$median, nrow(dat))
        stanGPbeta_sigma[,n] <- rep(parssummaryGP2[parssummaryGP2$variable=="sigma_e",]$median, nrow(dat))
        #parei aqui
        #params_stan_rb<-rstan::extract(stan_rb)
        
        #logastangpb <- paramsG2[,grep("^b",names(paramsGP2))]
        #umsystangp <- .5 *  logastangpa - 0.07 * ( logastangpa)^2
        # 
        #Smsystangp <- matrix(NA, nrow=nrow(logastangpa), ncol=ncol(logastangpa))
        #Sgenstangp <- matrix(NA, nrow=nrow(logastangpa), ncol=ncol(logastangpa))
        #
        #for(j in 1:ncol(Smsystangp)){
        #  Smsystangp[,j] <- (1 - gsl::lambert_W0(exp(1 - logastangpa[[j]] )))/paramsGP1$b
        #  Sgenstangp[,j] <- unlist(mapply(sGenSolverdlm,a=logastangpa[[j]],
        #     Smsy=Smsystangp[,j], b=paramsGP1$b))
        #}
    
        #stanGPalpha_umsy[,n] <- apply(umsystangp,2,median)
        #stanGPalpha_smsy[,n] <- apply(Smsystangp,2,median)
        #stanGPalpha_sgen[,n] <- apply(Sgenstangp,2,median)
     
        #stanGPalpha_convergence_alpha[,n] <- as.numeric(abs(parssummaryGP1[grep("log_a",parssummaryGP1$variable),]$rhat-1)>.1)       
        #stanGPalpha_convergence_beta[n] <- as.numeric(abs(parssummaryGP1[grep("^b",parssummaryGP1$variable),]$rhat-1)>.1)
        #stanGPalpha_convergence_sigma[n] <- as.numeric(abs(parssummaryGP1[parssummaryGP1$variable=="sigma_e",]$rhat-1)>.1)
   
        
      },  error = function(e) {
        print(paste("stan RB b failed in iteration",n))
      })



      #need to add
      #HM regime shift model jags a,b and and and b
      #HM regime shift stan  a,b and and and b

      


    }

    dlm_alpha <- list(alpha=dlmalpha_alpha,
    beta=dlmalpha_beta, 
    sigma=dlmalpha_sigma,
    smsy=dlmalpha_smsy,
    umsy=dlmalpha_umsy,
    sgen=dlmalpha_sgen,
    
    nll=dlmalpha_nll,
    AIC=dlmalpha_AIC, 
    BIC=dlmalpha_BIC, 
    convergence=dlmalpha_convergence )

    dlm_beta <- list(alpha=dlmbeta_alpha,
    beta=dlmbeta_beta, 
    sigma=dlmbeta_sigma,
    smsy=dlmbeta_smsy,
    umsy=dlmbeta_umsy,
    sgen=dlmbeta_sgen,

    nll=dlmbeta_nll, 
    AIC=dlmbeta_AIC, 
    BIC=dlmbeta_BIC, 
    convergence=dlmbeta_convergence )

    dlm_alphabeta <- list(alpha=dlmalphabeta_alpha,
    beta=dlmalphabeta_beta, 
    sigma=dlmalphabeta_sigma,
    smsy=dlmalphabeta_smsy,
    umsy=dlmalphabeta_umsy,
    sgen=dlmalphabeta_sgen,
    nll= dlmalphabeta_nll,
    AIC=dlmalphabeta_AIC, 
    BIC=dlmalphabeta_BIC, 
    convergence=dlmalphabeta_convergence )

    tmbKF <- list(alpha= tmbKF_alpha,
      beta=tmbKF_beta,
      sigma=tmbKF_sigma,
      smsy=tmbKF_smsy,
      umsy=tmbKF_umsy,
      sgen=tmbKF_sgen,
      nll=tmbKF_nll,
      AIC=tmbKF_AIC,
      BIC=tmbKF_BIC,
      convergence= tmbKF_convergence)

    tmbRB <- list(alpha= tmbRB_alpha,
      beta=tmbRB_beta,
      sigma=tmbRB_sigma,
      smsy=tmbRB_smsy,
      umsy=tmbRB_umsy,
      sgen=tmbRB_sgen,
      nll=tmbRB_nll,
      AIC=tmbRB_AIC,
      BIC=tmbRB_BIC,
      convergence= tmbRB_convergence)

    stanRBalpha <- list(alpha= stanRBalpha_alpha,
      beta= stanRBalpha_beta,
      sigma=stanRBalpha_sigma,
      smsy=stanRBalpha_smsy,
      umsy=stanRBalpha_umsy,
      sgen=stanRBalpha_sgen,
      loo=stanRBalpha_loo,
      convergence_alpha= stanRBalpha_convergence_alpha,
      convergence_beta= stanRBalpha_convergence_beta,
      convergence_sigma= stanRBalpha_convergence_sigma
      )
    
    stanRBbeta <- list(alpha= stanRBbeta_alpha,
                        beta= stanRBbeta_beta,
                        sigma=stanRBbeta_sigma,
                        smsy=stanRBbeta_smsy,
                        umsy=stanRBbeta_umsy,
                        sgen=stanRBbeta_sgen,
                        loo=stanRBbeta_loo,
                        
                        convergence_alpha= stanRBbeta_convergence_alpha,
                        convergence_beta= stanRBbeta_convergence_beta,
                        convergence_sigma= stanRBbeta_convergence_sigma
    )

    stanRBalphabeta <- list(alpha= stanRBalphabeta_alpha,
                        beta= stanRBalphabeta_beta,
                        sigma=stanRBalphabeta_sigma,
                        smsy=stanRBalphabeta_smsy,
                        umsy=stanRBalphabeta_umsy,
                        sgen=stanRBalphabeta_sgen,
                        loo=stanRBalphabeta_loo,
                        
                        convergence_alpha= stanRBalphabeta_convergence_alpha,
                        convergence_beta= stanRBalphabeta_convergence_beta,
                        convergence_sigma= stanRBalphabeta_convergence_sigma
    )

    stanGPalpha <- list(alpha= stanGPalpha_alpha,
                        beta= stanGPalpha_beta,
                        sigma=stanGPalpha_sigma,
                        smsy=stanGPalpha_smsy,
                        umsy=stanGPalpha_umsy,
                        sgen=stanGPalpha_sgen,
                        loo=stanGPalpha_loo,
                        
                        convergence_alpha= stanGPalpha_convergence_alpha,
                        convergence_beta= stanGPalpha_convergence_beta,
                        convergence_sigma= stanGPalpha_convergence_sigma
    )

    stanGPbeta <- list(alpha= stanGPbeta_alpha,
                        beta= stanGPbeta_beta,
                        sigma=stanGPbeta_sigma,
                        smsy=stanGPbeta_smsy,
                        umsy=stanGPbeta_umsy,
                        sgen=stanGPbeta_sgen,
                        loo=stanGPbeta_loo,
                        
                        convergence_alpha= stanGPbeta_convergence_alpha,
                        convergence_beta= stanGPbeta_convergence_beta,
                        convergence_sigma= stanGPbeta_convergence_sigma
    )

    stanGPalphabeta <- list(alpha= stanGPalphabeta_alpha,
                        beta= stanGPalphabeta_beta,
                        sigma=stanGPalphabeta_sigma,
                        smsy=stanGPalphabeta_smsy,
                        umsy=stanGPalphabeta_umsy,
                        sgen=stanGPalphabeta_sgen,
                        loo=stanGPalphabeta_loo,
                        
                        convergence_alpha= stanGPalphabeta_convergence_alpha,
                        convergence_beta= stanGPalphabeta_convergence_beta,
                        convergence_sigma= stanGPalphabeta_convergence_sigma
    )

    
    
    estData[[i]]<-list(dlm_alpha=dlm_alpha, 
      dlm_beta=dlm_beta, 
      dlm_alphabeta=dlm_alphabeta,
      tmbKF=tmbKF,
      tmbRB=tmbRB,
      stanRBalpha=stanRBalpha,
      stanRBbeta=stanRBbeta,
      stanRBalphabeta=stanRBalphabeta)
  }


  return(list(simData=simData,
  estData=estData))

}



 
#' calcpbias
#' 
#' Calculate percent bias in main recruitment parameters and management quantities
#' 
#' 
#'@param simData output of sim_est function
#'@param estData output of sim_est function
#'
#' 
#' 
#' 
#'@return a data frame for plotting percent bias of recruitment parameters and management quantities
#' 
calcpbias<-function(simData, estData){
  
  #simData<-out$simData
  #estData<-out$estData
  allpb<- list()

  for(i in seq_len(length(estData))){


    #pbias alpha
    nsimyr <- length(unique(simData[[i]]$byr))
    #alpha_sim <- simData[[i]]$alpha_true
    #beta_sim <- simData[[i]]$beta_true
    #sig_sim <- simData[[i]]$sigma
    #smsy_sim <- simData[[i]]$sMSY
    #umsy_sim <- simData[[i]]$uMSY
    #sgen_sim <- simData[[i]]$sGen
    
    #need o add the anyna condition here. 
    length(simData[[i]][["alpha_true"]])

    params <- c( "alpha","beta","sigma","smsy","umsy","sgen")
    simparams <- c( "alpha_true","beta_true","sigma","sMSY","uMSY","sGen")
    #alpha
    
    pbdflist <- list()

    for(m in seq_along(params)){
      o<-pbiasdfdlm(estx=estData[[i]]$dlm_alpha, 
        simx=simData[[i]][[simparams[m]]], 
        estmodel="dlm a vary", 
        parameter=params[m], 
        nsimyr=nsimyr)
      oo<-pbiasdfdlm(estx=estData[[i]]$dlm_beta, 
        simx=simData[[i]][[simparams[m]]], 
        estmodel="dlm b vary", 
        parameter=params[m], 
        nsimyr=nsimyr)
      ooo<-pbiasdfdlm(estx=estData[[i]]$dlm_alphabeta, 
        simx=simData[[i]][[simparams[m]]], 
        estmodel="dlm a and b vary", 
        parameter=params[m],    
        nsimyr=nsimyr)
      o4<- pbiasdfdlm(estx=estData[[i]]$tmbKF, 
        simx=simData[[i]][[simparams[m]]], 
        estmodel="TMB Kalman Filter a vary", 
        parameter=params[m],    
        nsimyr=nsimyr)
      o5<- pbiasdfdlm(estx=estData[[i]]$tmbRB, 
        simx=simData[[i]][[simparams[m]]], 
        estmodel="TMB Recursive Bayes a vary", 
        parameter=params[m],    
        nsimyr=nsimyr)
    
      o6<- pbiasdfdlm(estx=estData[[i]]$stanRBalpha, 
                      simx=simData[[i]][[simparams[m]]], 
                      estmodel="stan Recursive Bayes a vary", 
                      parameter=params[m],    
                      nsimyr=nsimyr)
      o7<- pbiasdfdlm(estx=estData[[i]]$stanRBbeta, 
                      simx=simData[[i]][[simparams[m]]], 
                      estmodel="stan Recursive Bayes b vary", 
                      parameter=params[m],    
                      nsimyr=nsimyr)
      
      
      pbdflist[[m]]<-rbind(o,oo,ooo,o4,o5,o6)
    }
  
    pbdf<-do.call("rbind",pbdflist)
    pbdf$scenario<-simData[[i]]$scenario 
   
    allpb[[i]] <- pbdf
   
  }

  dfall<-do.call("rbind",allpb)

  return(dfall)

}




 
#' pbiasdfdlm
#' 
#' Aux function to calculate percent bias in main recruitment parameters and management quantities for dlm estimates
#' 
#' 
#'@param estx list containing estimated quantities
#'@param simx list containing simulated quantities
#'@param estmodel string to identify the estimation model
#'@param nsimyr number of simulation years
#' 
#' 
#' 
#'@return a data frame for plotting percent bias of recruitment parameters and management quantities
#' 
pbiasdfdlm <- function(estx, simx, estmodel, parameter, nsimyr){
  
 

  dlma<-reshape::melt(estx[[parameter]], varnames=c('byr', 'iteration'))
  
    dlma$pbias<-(dlma$value-simx)/simx
    dlma$estmodel<- estmodel
    dlma$parameter<- parameter
    dlma$convergence <- rep(estx[["convergence"]], each=nsimyr)
    dlma$true_value <- simx
  
   

   return(dlma)
}