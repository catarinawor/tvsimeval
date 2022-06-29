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



out<-sim_est(simPar=simPar, cuPar=cuPar, srDat=srDat, ricPars=ricPars, corrmat=corrmat, outDir="example", iteration=2)

df<-calcpbias(simData=out$simData,estData= out$estData)

dfc<-df[df$convergence==0,]

ggplot(dfc,aes(x=estmodel, y=pbias)) +
geom_boxplot()+
coord_cartesian(ylim = c(-1,1))+
geom_hline(yintercept=0) +
theme_bw(14)+
#stat_summary(fun.data = give.n, geom = "text", hjust = 0.5,
#    vjust = 0.9)+
facet_wrap(~parameter)



## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simPar$species),sep = "_"))

## First check to ensure that a single scenario can be run (only a small number
# of trials necessary)
plotscn <- TRUE
p <- list()
simData <- list()
estData <- list()

iteration<-2

for(i in seq_len(nrow(simPar))){

  genericRecoverySim(simPar=simPar[i,], cuPar=cuPar, catchDat=NULL, srDat=srDat,
            variableCU=FALSE, ricPars=ricPars , larkPars=NULL,cuCustomCorrMat= corrmat,
            outDir="example", nTrials=iteration, makeSubDirs=TRUE, random=FALSE)
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
    select(scenario, everything())  

  dlm_alpha <- list()
  dlm_beta <- list() 
  dlm_alphabeta <- list() 

  for( n in seq_len(iteration)){

    dat <- simData[[i]]  %>% 
        dplyr::filter( iteration==n, CU==1) %>% 
        dplyr::select(byr, obsSpawners,obsRecruits ) %>%
        dplyr::rename(spwn=obsSpawners, rec=obsRecruits)

    #lnRS <- log(dat$rec/dat$spwn)
    #lm(lnRS ~dat$spwn)
    #plot(dat$spwn,lnRS)    
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

 

#write routine to compare
#parameter bias
# management quantities biases
#BIC and AIC summaries. 






ggsave(
      filename = "../plots/scenarios.pdf", 
      plot = marrangeGrob(p, nrow=1, ncol=1), 
      width = 12, height = 5
    )
  



#===================================================
#estimation routines

iter <- unique(simData[[1]]$iteration)
nsc <- length(scenNames)

























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


