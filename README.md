

## SPMpriors
##### Henning Winker (EC-JRC, 2020)

This package aims to provide a tool box for generating priors for stock assessments from FishLife (Thorson et al. 2020) unsing Multivariate-Normal (MVN) Monte-Carlo simulations, while allowing to impose stock specific prior tuning criteria. The main focus is to translate stock parameters into surplus production model priors r, shape, K. In addition other useful quantaties, such as Generation time, steepness h, sigmaR or recruitment auto-correlation, can be objectively obtained. It is recommended to only use the r priors from FishLife in a Schafer model. To translate the MVN parameters into Pella-Tomlinson type parameterization an age-structured equilibrium approach is proposed (c.f. Winker et al. 2020), which require an approximation of the length at first capture a additional input (coming soon).  

Installing SPMpriors requires the `library(devtools)` and `library(FishLife)` prior to the installation of `library(SPMpriors)` using the following steps:
<br/>

`install.packages('devtools')`

`devtools::install_github("james-thorson/FishLife")`

`install_github("henning-winker/SPMpriors")`


`library(SPMpriors)`

<br/>


## Example 
Get Flounder traits

`stk = mvn_traits(Genus="Platichthys",Species="flesus",h=c(0.6,0.9),Plot=T)`

check

`stk$traits`

Tune to North Sea Loo and Lm 

`stkLinfK = mvn_traits(Genus="Platichthys",Species="flesus",Loo=c(41,0.1),Lm=c(21,0.1),h=c(0.6,0.9),Plot=T)`

<img src="https://github.com/Henning-Winker/JARA/blob/master/JARAplotting/trj2x1.png" width = "800" >

<br/>

<img src="https://github.com/Henning-Winker/JARA/blob/master/JARAplotting/trj2x1.png" width = "800" >


Note: Leslie r should be only used in a Schaefer model 

Age-Structure Equilibrium Model will follow for Pella-Tomlinson estimates

`aspm_priors(stkLinfK$mvnsim,Lc)` ...coming soon


