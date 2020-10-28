 
fl2asem <- function(mvnstk,mc=1000,t0=-0.5,aW=0.01,bW=3.04,dLm=NULL,tmin=0,tmax,Lc=NULL,dLc=NULL,Fvec=NULL,R0=1,Plot=FALSE){
  steps = floor(nrow(mvnstk)/mc)
  thinning = seq(1,1000*steps,steps)
  mvnstk = mvnstk[thinning,]
  if(is.null(Fvec)) Fvec = c(0,exp(seq(-5,1.5,0.005)))
  
  pt_prior = NULL
  
  for(i in 1:nrow(mvnstk)){
  mv = mvnstk[i,]
  if(is.null(dLm)) dLm = 0.1*Lm
  if(is.null(Lc)) Lc = Lm 
  if(is.null(dS)) dS = 0.1*Lc
  mv$tmax = max(round(mv$tmax,0),2)
  pti = asem(Loo=mv$Loo,K=mv$K,t0=t0,aW=0.01,bW=3.04,Lm=mv$Lm,dLm=dLm,tmin=0,tmax=mv$tmax,h=mv$h,M=mv$M,Lc=Lc,dLc=dLc,Fvec=Fvec,R0=R0)
  pt_prior = rbind(pt_prior,pti$spm)   
  }
  
}
  