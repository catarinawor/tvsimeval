#samSim updates
#Catarina Wor
#June 2022 
#=============================================================

#install samsim 
#remotes::install_github("Pacific-salmon-assess/samSim", force=TRUE)


#or locally

#setwd("C:\\Users\\worc\\Documents\\timevar\\samSim\\R")


devtools::document()
devtools::load_all()
#devtools::build()

#detach("package:samSim", unload=TRUE)
#devtools::build(pkg = "C:/Users/worc/Documents/LRP/samSim/.")
#devtools::install_local("C:/Users/worc/Documents/LRP/samSim_0.0.0.9.tar.gz")

#library(samSim)
?simParexample

?synchList


usethis::use_data(cuParexample,  internal = TRUE)


