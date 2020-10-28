#' fl2asem()
#'
#' Translates MVN stocks pars from FishLife (Thorson et al. 2020) 
#' into Pella-Tomlison SPM priors r, shape, bmsy/k (JABBA/SPicT type) 
#' through an Age-Structured Equilibrium Model (Winker et al. 2020) 
#' @param flstk output object from flmvn_traits()   
#' @param t0 theoretical age a length zero (default -0.5)
#' @param aW parameter of the lenght-weight relation (default 0.01)
#' @param bW parameter of the lenght-weight relation (default 3.04)
#' @param dLm width of ogive (default=0.1*Lm)   
#' @param tmin minimum age (default = 0)
#' @param Lc length at first capture (default = Lm)
#' @param dLc width of ogive (default = 0.1*Lc)
#' @param Fvec vector of instantanious fishing fortality
#' @param R0 unfished recruitment default = 1
#' @param plot.progress plot simulated production function
#' @return Pella-Tomlison Prior 
#' @export
#' @author Henning Winker (JRC-EC)

 
fl2asem <- function(flstk,mc=1000,t0=-0.5,aW=0.01,bW=3.04,dLm=NULL,tmin=0,tmax,Lc=NULL,dLc=NULL,Fvec=NULL,R0=1,plot.progress=FALSE){
  mvnstk = flstk$mvnstk
  steps = floor(nrow(mvnstk)/mc)
  thinning = seq(1,1000*steps,steps)
  mvnstk = mvnstk[thinning,]
  if(is.null(Fvec)) Fvec = c(0,exp(seq(-5,1.5,0.005)))
  
  pt_prior = NULL
  
  if(plot.progress){
    mu = flstk$traits
    hat = asem(Loo=mu[mu$trait=="Loo",4],K=mu[mu$trait=="K",4],t0=t0,aW=aW,bW=bW,Lm=mu[mu$trait=="Lm",4],dLm=dLm,tmin=0,tmax=round(mu[mu$trait=="tmax",4],0),h=mu[mu$trait=="h",4],M=mu[mu$trait=="M",4],Lc=Lc,dLc=dLc,Fvec=Fvec,R0=R0)
    plot(hat$eqdyn$bk,hat$eqdyn$yield,lwd=3,type="l",col=4,xlab="b/k",ylab="relative yield",ylim=c(0,max(hat$eqdyn$yield)*1.6))
  } 
  for(i in 1:nrow(mvnstk)){
  mv = mvnstk[i,]
  if(is.null(dLm)) dLm = 0.1*Lm
  if(is.null(Lc)) Lc = Lm 
  if(is.null(dS)) dS = 0.1*Lc
  mv$tmax = max(round(mv$tmax,0),2)
  pti = asem(Loo=mv$Loo,K=mv$K,t0=t0,aW=0.01,bW=3.04,Lm=mv$Lm,dLm=dLm,tmin=0,tmax=mv$tmax,h=mv$h,M=mv$M,Lc=Lc,dLc=dLc,Fvec=Fvec,R0=R0)
  if(plot.progress) lines(pti$eqdyn$bk,pti$eqdyn$yield,col=gray(runif(1,0.5,0.8),0.4))
  pt_prior = rbind(pt_prior,pti$spm)   
  }
  if(plot.progress){
    lines(hat$eqdyn$bk,hat$eqdyn$yield,lwd=2,col=4) 
    lines(rep(hat$spm$bmsyk,2),c(0,hat$spm["msy"]),lty=2,col=4)
  }  
  mu = exp(apply(log(pt_prior),2,median))
  logsd =  apply(log(pt_prior),2,sd)
  priors= rbind(mu,logsd)[,-3]
  return(priors)
}
  


