#==========================================
#utility functions
#==========================================




give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
} 




sumpair <- function(x,k=2){
  #x<-as.numeric(abs(bhmma$mcmcsummary[grep("^gamma\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1)

  y<-matrix(x,ncol=k)
  return(apply(y,1,sum))
  
  # experiment with the multiplier to find the perfect position
} 

