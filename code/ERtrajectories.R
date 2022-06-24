#=============================================================
#Visual inspection of Empirical exploitation rates to inform simulation scenarios
#
#Catarina Wor
# April 2022 
#=============================================================

library(here)
library(ggplot2)
#read in data

# chum and Pink datasets from Mike Malik
chum <- read.csv(here("data","HRchumpink","chum_catch_v1.1.csv"))
pink <- read.csv(here("data","HRchumpink","pink_catch_v1.1.csv"))


cphr<-list()

for(i in seq_along(unique(chum$region))){
  ch1 <- chum[chum$region==unique(chum$region)[i],]

  pc <- ggplot(ch1)+
    geom_point(aes(x=return.yr,y=harvest.rate ))+
    geom_smooth(aes(x=return.yr,y=harvest.rate), method="lm")+
    facet_wrap(~stock)+
    coord_cartesian(ylim=c(0,1)) +ggtitle(paste("Chum stocks in",unique(chum$region)[i]))+
    theme_bw(14)
 
 cphr[[i]]<-pc
}


ggsave(
   filename = "../plots/hr_explore/plot_chumhr.pdf", 
   plot = marrangeGrob(cphr, nrow=1, ncol=1), 
   width = 15, height = 9
)


head(pink)


pphr<-list()

for(i in seq_along(unique(pink$region))){
  ph1 <- pink[pink$region==unique(pink$region)[i],]

  pc <- ggplot(ph1)+
    geom_point(aes(x=return.yr,y=harvest.rate ))+
    geom_smooth(aes(x=return.yr,y=harvest.rate), method="lm")+
    facet_wrap(~stock)+
    coord_cartesian(ylim=c(0,1)) +ggtitle(paste("Pink stocks in",unique(pink$region)[i]))+
    theme_bw(14)
 
 pphr[[i]]<-pc
}


ggsave(
   filename = "../plots/hr_explore/plot_pinkhr.pdf", 
   plot = marrangeGrob(pphr, nrow=1, ncol=1), 
   width = 15, height = 9
)