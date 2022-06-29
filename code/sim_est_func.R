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
 

#install samsim 
#remotes::install_github("Pacific-salmon-assess/samSim", force=TRUE)

library(samSim)
source("dlm-wrapper.R")



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
sim_est <- function(simPar, cuPar, srDat, ricPars, corrmat, outDir, iteration){

  ## Store relevant object names to help run simulation 
  scenNames <- unique(simPar$scenario)
  dirNames <- sapply(scenNames, function(x) paste(x, unique(simPar$species),sep = "_"))

  plotscn <- TRUE
  p <- list()
  simData <- list()
  estData <- list()



  for(i in seq_len(nrow(simPar))){

    genericRecoverySim(simPar=simPar[i,], cuPar=cuPar, catchDat=NULL, srDat=srDat,
              variableCU=FALSE, ricPars=ricPars , larkPars=NULL,cuCustomCorrMat= corrmat,
              outDir="example", nTrials=iteration, makeSubDirs=TRUE, random=FALSE)
    #simPar=simPar[1,];catchDat=NULL;variableCU=FALSE;larkPars=NULL;cuCustomCorrMat= corrmat;
    #outDir="example"; nTrials=1; makeSubDirs=TRUE; random=FALSE;larkPars=NULL;
    #erCorrMat=NULL;uniqueProd=TRUE;uniqueSurv=FALSE

    simData[[i]] <- readRDS(here("example","SamSimOutputs","simData", simPar$nameOM[i],simPar$scenario[i],
                         paste(simPar$nameOM[i],"_", simPar$nameMP[i], "_", "CUsrDat.RData",sep="")))$srDatout
  


    nyr <- max(unique(simData[[i]]$year))

    simData[[i]] <- simData[[i]] %>% 
      filter(CU == 1, year %in% (nyr-50+1):nyr) %>% 
      filter(!is.na(obsRecruits))%>%
      mutate() %>% 
      rename(byr=year, spwn=spawners, rec=recruits, alpha_true=alpha, beta_true=beta) %>% 
      mutate(scenario = scenNames[i]) %>% # cols for output
      select(scenario, everything())  

    dlm_alpha <- list()
    dlm_beta <- list() 
    dlm_alphabeta <- list() 

    for( n in seq_len(iteration)){

      dat <- simData[[i]]  %>% 
        dplyr::filter( iteration==n, CU==1) %>% 
        dplyr::select(byr, obsSpawners,obsRecruits ) %>%
        dplyr::rename(spwn=obsSpawners, rec=obsRecruits)

       
      dlm_alpha[[n]] <- NA
      dlm_beta[[n]] <- NA
      dlm_alphabeta[[n]] <- NA
      tryCatch({
         dlm_alpha[[n]] <-fitDLM(data = dat,
                            alpha_vary = TRUE,
                            beta_vary = FALSE)
      },  error = function(e) {
        print(paste("dlm_alpha failed in iteration",n))
      })
  
      tryCatch({
        # beta varies in estimation model
        dlm_beta[[n]] <- fitDLM(data = dat,
                            alpha_vary = FALSE,
                            beta_vary = TRUE)
      },  error = function(e) {
        print(paste("dlm_beta failed in iteration",n))
      })
    
      tryCatch({
        # alpha and beta vary in estimation model
        dlm_alphabeta[[n]] <- fitDLM(data = dat,
                            alpha_vary = TRUE,
                            beta_vary = TRUE)
      },  error = function(e) {
        print(paste("dlm_alphabeta failed in iteration",n))
      })

    }


    estData[[i]]<-list(dlm_alpha=dlm_alpha, dlm_beta=dlm_beta, dlm_alphabeta=dlm_alphabeta)
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

  for(i in seq_len(length(estData))){


    #pbias alpha
    alpha_sim <- simData[[i]]$alpha_true
    beta_sim <- simData[[i]]$beta_true
    umsy_sim <- simData[[i]]$uMSY
    smsy_sim <- simData[[i]]$sMSY
    sgen_sim <- simData[[i]]$sGen
    sig_sim <- simData[[i]]$sigma
    
    
    alpha_dlmavary <- do.call("rbind",lapply(estData[[i]]$dlm_alpha, function(x){x$results}))$alpha
    alpha_dlmbvary <- do.call("rbind",lapply(estData[[i]]$dlm_beta, function(x){x$results}))$alpha
    alpha_dlmabvary <- do.call("rbind",lapply(estData[[i]]$dlm_alphabeta, function(x){x$results}))$alpha


    convdlm <- c(do.call("c",lapply(estData[[i]]$dlm_alpha, function(x){rep(x$convergence,nrow(x$results))})),
        do.call("c",lapply(estData[[i]]$dlm_beta, function(x){rep(x$convergence,nrow(x$results))})),
        do.call("c",lapply(estData[[i]]$dlm_alphabeta, function(x){rep(x$convergence,nrow(x$results))})))
     
  
    
    pbdfa <- data.frame(iteration=rep(simData[[i]]$iteration,3),
      estmodel=rep(c("dlm a vary","dlm b vary","dlm a and b vary"), each=length(alpha_sim)),
      parameter= "alpha",
      pbias=c((alpha_dlmavary-alpha_sim)/alpha_sim,
        (alpha_dlmbvary -alpha_sim)/alpha_sim,
        (alpha_dlmabvary-alpha_sim)/alpha_sim),
      convergence= convdlm
    )


    #pbias beta
    beta_dlmavary <- -do.call("rbind",lapply(estData[[i]]$dlm_alpha, function(x){x$results}))$beta
    beta_dlmbvary <- -do.call("rbind",lapply(estData[[i]]$dlm_beta, function(x){x$results}))$beta
    beta_dlmabvary <- -do.call("rbind",lapply(estData[[i]]$dlm_alphabeta, function(x){x$results}))$beta
 

    pbdfb <- data.frame(iteration=rep(simData[[i]]$iteration,3),
      estmodel=rep(c("dlm a vary","dlm b vary","dlm a and b vary"), each=length(beta_sim)),
      parameter= "beta",
      pbias=c((beta_dlmavary-beta_sim)/beta_sim,
        (beta_dlmbvary -beta_sim)/beta_sim,
        (beta_dlmabvary-beta_sim)/beta_sim),
      convergence= convdlm
    )

    #pbias smsy
    smsy_dlmavary <- do.call("rbind",lapply(estData[[i]]$dlm_alpha, function(x){x$results}))$Smsy
    smsy_dlmbvary <- do.call("rbind",lapply(estData[[i]]$dlm_beta, function(x){x$results}))$Smsy
    smsy_dlmabvary <- do.call("rbind",lapply(estData[[i]]$dlm_alphabeta, function(x){x$results}))$Smsy
 

    pbdfsmsy <- data.frame(iteration=rep(simData[[i]]$iteration,3),
      estmodel=rep(c("dlm a vary","dlm b vary","dlm a and b vary"), each=length(smsy_sim)),
      parameter= "Smsy",
      pbias=c((smsy_dlmavary-smsy_sim)/smsy_sim,
        (smsy_dlmbvary -smsy_sim)/smsy_sim,
        (smsy_dlmabvary-smsy_sim)/smsy_sim),
      convergence= convdlm)


    #pbias umsy
    umsy_dlmavary <- do.call("rbind",lapply(estData[[i]]$dlm_alpha, function(x){x$results}))$Umsy
    umsy_dlmbvary <- do.call("rbind",lapply(estData[[i]]$dlm_beta, function(x){x$results}))$Umsy
    umsy_dlmabvary <- do.call("rbind",lapply(estData[[i]]$dlm_alphabeta, function(x){x$results}))$Umsy
 

    pbdfumsy <- data.frame(iteration=rep(simData[[i]]$iteration,3),
      estmodel=rep(c("dlm a vary","dlm b vary","dlm a and b vary"), each=length(beta_sim)),
      parameter= "Umsy",
      pbias=c((umsy_dlmavary-umsy_sim)/Umsy_sim,
        (umsy_dlmbvary -umsy_sim)/umsy_sim,
        (umsy_dlmabvary-umsy_sim)/umsy_sim),
      convergence= convdlm)


    #pbias sgen
    sgen_dlmavary <- do.call("rbind",lapply(estData[[i]]$dlm_alpha, function(x){x$results}))$Sgen
    sgen_dlmbvary <- do.call("rbind",lapply(estData[[i]]$dlm_beta, function(x){x$results}))$Sgen
    sgen_dlmabvary <- do.call("rbind",lapply(estData[[i]]$dlm_alphabeta, function(x){x$results}))$Sgen
 

    pbdfsgen <- data.frame(iteration=rep(simData[[i]]$iteration,3),
      estmodel=rep(c("dlm a vary","dlm b vary","dlm a and b vary"), each=length(sgen_sim)),
      parameter= "Sgen",
      pbias=c((sgen_dlmavary-sgen_sim)/sgen_sim,
        (sgen_dlmbvary -sgen_sim)/sgen_sim,
        (sgen_dlmabvary-sgen_sim)/sgen_sim),
      convergence= convdlm)

    #pbias sig
    sig_dlmavary <- do.call("c",lapply(estData[[i]]$dlm_alpha, function(x){rep(x$sd.est[1],nrow(x$result))}))
    sig_dlmbvary <- do.call("c",lapply(estData[[i]]$dlm_beta, function(x){rep(x$sd.est[1],nrow(x$result))}))
    sgen_dlmabvary <- do.call("c",lapply(estData[[i]]$dlm_alphabeta,  function(x){rep(x$sd.est[1],nrow(x$result))}))
    
    pbdfsig <- data.frame(iteration=rep(simData[[i]]$iteration,3),
      estmodel=rep(c("dlm a vary","dlm b vary","dlm a and b vary"), each=length(sig_sim)),
      parameter= "sigma",
      pbias=c((sig_dlmavary-sigma_sim)/sigma_sim,
        (sig_dlmavary -sigma_sim)/sigma_sim,
        (sig_dlmavary-sigma_sim)/sigma_sim),
      convergence= convdlm)
 



    #pbias siga -- not really a parameter of interest.
     
     df<- rbind(pbdfa,
      pbdfb,
      pbdfsmsy,
      pbdfumsy,
      pbdfsgen,
      pbdfsig)  

    df$scn<-simData[[i]]$scenario 

   

  }

 dfall<-do.call("rbind",df)

 return(df)


}

