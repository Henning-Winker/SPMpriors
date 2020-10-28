
#*********************************************************************
# SPMpriors: R package for deriving Schaefer and Pella-Tomlinson 
# from FishLife (Thorson 2020) through a MVN Age-Structured Monte-Carlo
# simulaiton approach (Winker et al. 2020)
# Developed by Henning Winker
# JRC-EC, Ispra, 2020
#*********************************************************************


#********************************
# Installation Instructions
#********************************
install.packages('devtools')
devtools::install_github("james-thorson/FishLife")
devtools::install_github("henning-winker/SPMpriors")
#********************************

library(SPMpriors)
library(FishLife)

#><> North Sea Flounder 
# Get MVN stock par replicates from FishLife, while tuning Loo, Lm and h
stk = flmvn_traits(Genus="Platichthys",Species="flesus",Loo=c(41,0.1),Lm=c(21,0.1),h=c(0.6,0.9),Plot=T,savepng = F)

stk$traits

# the r prior can be used Schaefer SPM, but should not be applied in
# Pella-Tomlison SPMs. For the latter it is more appropriate to approximate
# r and shape from an Age-Structured Equilibrium Model (Winker et al. 2020)

# The default assumption is that length at first capture = Lm
fl2asem(stk,mc=1000,plot.progress = T)

#><> White Angler Fish
stk = flmvn_traits(Genus="Lophius",Species="piscatorius",tmax=c(20,0.2),h=c(0.6,0.9),Plot=T,savepng = F)

# assume Lc = Lm
fl2asem(stk,mc=1000,plot.progress = T)

# what if Lc < Lm
fl2asem(stk,mc=1000,Lc=30,plot.progress = T)






