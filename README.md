

## SPMpriors
##### Henning Winker (Joint Research Centre of the European Commission, Ispra, 2020)

This package aims to provide a tool box for generating priors for stock assessments from FishLife (Thorson et al. 2020) unsing Multivariate-Normal (MVN) Monte-Carlo simulations, while allowing to impose stock specific prior tuning criteria. The main focus is to translate stock parameters into surplus production model priors r, shape, K. In addition other useful quantaties, such as Generation time, steepness h, sigmaR or recruitment auto-correlation, can be objectively obtained. It is recommended to only use the r priors from FishLife in a Schafer model. To translate the MVN parameters into Pella-Tomlinson type parameterization an age-structured equilibrium approach is proposed (c.f. Winker et al. 2020), which require an approximation of the length at first capture a additional input (coming soon).  

Installing SPMpriors requires the `library(devtools)` and `library(FishLife)` prior to the installation of `library(SPMpriors)` using the following steps:
<br/>

`install.packages('devtools')`

`devtools::install_github("james-thorson/FishLife")`

`install_github("henning-winker/SPMpriors")`


`library(SPMpriors)`

<br/>


## Example 
Get Flounder traits tuned  to North Sea Loo and Lm 

`stk = flmvn_traits(Genus="Platichthys",Species="flesus",Loo=c(41,0.1),Lm=c(21,0.1),h=c(0.6,0.9),Plot=T)`

<img src="https://github.com/Henning-Winker/SPMpriors/blob/main/Example/Platichthys.flesus.fl_ellipse.png" width = "800" >

<br/>

<img src="https://github.com/Henning-Winker/SPMpriors/blob/main/Example/Platichthys.flesus.fl_stktraits.png" width = "800" >




Note that Leslie-matrix r prior from FishLife should can be used Schaefer SPM, but should not be applied in a Pella-Tomlison formaltion where the shape parameter differs from n = 2. 

For a Pella-Tomlison model it is more appropriate to approximate r and shape from an Age-Structured Equilibrium Model (Winker et al. 2020) The default assumption is that length at first capture = Lm:

`fl2asem(stk,mc=1000,plot.progress = T)`

Another example could be White Angler Fish

`stk = flmvn_traits(Genus="Lophius",Species="piscatorius",tmax=c(20,0.2),h=c(0.6,0.9),Plot=T,savepng = F)`

First assume Lc = Lm

`fl2asem(stk,mc=1000,plot.progress = T)`

..but what happens if Lc < Lm?

`fl2asem(stk,mc=1000,Lc=30)`

more plotting functions coming soon....



`fl2asem(stk,Lc=20)` ...coming soon


