#=============================================================
#Example use of samsim for time varying simulation evaluation
#using the Harrison Chinook data as to test
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




#here::here()
## Load relevant input data
# Simulation run parameters describing different scenarios
simPar <- read.csv("../data/samsimHarCk/harcnkSimPars.csv")
# CU-specific parameters
cuPar <- read.csv("../data/samsimHarCk/harcnkCUPars.csv")



## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simPar$species),sep = "_"))

## First check to ensure that a single scenario can be run (only a small number
# of trials necessary)


genericRecoverySim(simPar=simPar[1,], cuPar=cuPar, catchDat=NULL, srDat=NULL,
            variableCU=FALSE, ricPars=NULL , larkPars=NULL,cuCustomCorrMat= NULL,
            outDir="testing", nTrials=2, makeSubDirs=TRUE, random=FALSE)
#simPar=simPar[1,]; cuPar=cuPar; catchDat=NULL; srDat=NULL;
#variableCU=FALSE; ricPars=NULL; larkPars=NULL; cuCustomCorrMat= NULL;
#outDir="testing"; nTrials=2; makeSubDirs=TRUE; random=FALSE; uniqueProd=TRUE; uniqueSurv=FALSE




simData <- readRDS(here("testing","SamSimOutputs","simData", simPar$nameOM,simPar$scenario,
                         paste(simPar$nameOM,"_", simPar$nameMP, "_", "CUsrDat.RData",sep="")))$srDatout
  

df<-simData
df<-df[df$iteration==1,]
   

p <- ggplot(df) +
      geom_point(aes(x=spawners,y=recruits, col=year))+
      geom_text(aes(x=spawners,y=recruits,label=year, col=year))+
      theme_bw(14)+ theme(legend.position="none")+
      scale_colour_viridis_c() +
      labs(title = simPar$nameOM)
p


#TODO


#plug in the estimation routines

lm(log(df$recruits/df$spawners)~df$spawners)



for(i in 1:11){#seq_len(nrow(simPar))){

  genericRecoverySim(simPar=simPar[i,], cuPar=cuPar, catchDat=NULL, srDat=NULL,
            variableCU=FALSE, ricPars=NULL , larkPars=NULL,cuCustomCorrMat= NULL,
            outDir="testing", nTrials=2, makeSubDirs=TRUE, random=FALSE; uniqueSurv=FALSE)
  #simPar=simPar[1,];catchDat=NULL;variableCU=FALSE;larkPars=NULL;cuCustomCorrMat= corrmat;
  #outDir="example"; nTrials=1; makeSubDirs=TRUE; random=FALSE;larkPars=NULL;
  #erCorrMat=NULL;uniqueProd=TRUE;uniqueSurv=FALSE; uniqueProd=TRUE; uniqueSurv=FALSE

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

harsr<-read.csv("../data/samsimHarCk/HARSR.csv")


shar<-harsr$Natural.Origin.Spawners
lrshar<-log((harsr$AEQ_Recruitment..age.2.5.)/harsr$Natural.Origin.Spawners)

lmhar<-lm(lrshar~shar)
names(lmhar)
alpha<-lmhar$coefficients[[1]]
b<--lmhar$coefficients[[2]]

#179313.1
Seq<-alpha/b
S <- NULL
R <- NULL



A<- matrix(NA,ncol=5, nrow=40)
Aprime<- matrix(NA,ncol=5, nrow=40)
RET<- matrix(NA,ncol=5, nrow=40)


ma<- c(.0,.3,.5,.8,1)
M<-c(.4,.3,.2,.1,.1)
#M<-rep(0,5)

x<-rep(2,5)
 
aeq<-NULL
aeq[5]<-1
aeq[4]<-ma[4]+((1-ma[4])*(1-M[4])*aeq[5])
aeq[3]<-ma[3]+((1-ma[3])*(1-M[3])*aeq[4])
aeq[2]<-ma[2]+((1-ma[2])*(1-M[2])*aeq[3])
aeq[1]<-ma[1]+((1-ma[1])*(1-M[1])*aeq[2])

(x*aeq)/aeq



lxo<-NULL
#lxr<-NULL
 lxo[1]<-1
 #lxr[1]<-1*ma[1]
for(i in 2:5){
  lxo[i]<-lxo[i-1]*(1-ma[i-1])*(1-M[i-1])
  #lxr[i]<-lxo[i-1]*(1-M[i-1])*(1-ma[i-1])*ma[i]
 }
lxr<-lxo*ma
sum(alpha*lxr/b*lxr)




Seq


#hack if parameters were estimated with Aeq values (inflate rec to account for M and m)

for(y in 1:2){
S[y] <- Seq
R[y] <- Seq/sum(lxo*ma)

A[y,] <- lxo*R[y]
Aprime[y,]<- A[y,]
RET[y,]<- Aprime[y,]*ma

}

sum(RET[1,])

apply(RET,1,sum)


for(y in 3:40){
R[y] <- (exp(alpha)*S[y-2]*exp(-b*S[y-2]))/sum(lxo*ma)
A[y,1]<- R[y]
#fishery - all ages fully selected and ER=.3
Aprime[y,1]<- A[y,1]#*(1-.3)
RET[y,1]<- Aprime[y,1]*(ma[1])

  for(a in 2:5){
    A[y,a]<-  (Aprime[y-1,a-1] - RET[y-1,a-1])*(1-M[a-1])
    #fishery - all ages fully selected and ER=.3
    Aprime[y,a]<- A[y,a]#*(1-.3)
    RET[y,a]<- Aprime[y,a]*(ma[a])
  }

  S[y]<-sum(RET[y,])
}







#Using botsford equilibrium recruitment
#optional, adjust alpha if paremeter wer estimated based on aeq values
alpha1<-log(exp(alpha)/sum(lxo*ma))

S <- NULL
R <- NULL
A<- matrix(NA,ncol=5, nrow=40)
Aprime<- matrix(NA,ncol=5, nrow=40)
RET<- matrix(NA,ncol=5, nrow=40)


for(y in 1:2){
R[y] <- (alpha1+log(sum(lxr)))/(b*sum(lxr))

A[y,] <- lxo*R[y]
Aprime[y,]<- A[y,]
RET[y,]<- Aprime[y,]*ma
S[y]<-sum(RET[y,])
}

sum(RET[1,])

apply(RET,1,sum)


for(y in 3:40){
R[y] <- (exp(alpha1)*S[y-2]*exp(-b*S[y-2]))
A[y,1]<- R[y]
#fishery - all ages fully selected and ER=.3
Aprime[y,1]<- A[y,1]#*(1-.3)
RET[y,1]<- Aprime[y,1]*(ma[1])

  for(a in 2:5){
    A[y,a]<-  (Aprime[y-1,a-1] - RET[y-1,a-1])*(1-M[a-1])
    #fishery - all ages fully selected and ER=.3
    Aprime[y,a]<- A[y,a]#*(1-.3)
    RET[y,a]<- Aprime[y,a]*(ma[a])
  }

  S[y]<-sum(RET[y,])
}

