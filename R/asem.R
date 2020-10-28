#' asem()
#'
#' Age-Structured Equilibrium Model to compute ratio Fmsy= MSY/Bmsy and Bmsy/K
#' 
#' @param Loo  infinitive length (cm)  
#' @param K  brody  growth coefficient K   
#' @param t0 theoretical age a length zero (default -0.5)
#' @param aW parameter of the lenght-weight relation (default 0.01)
#' @param bW parameter of the lenght-weight relation (default 3.04)
#' @param Lm length-at-50-maturity (cm)   
#' @param dLm width of ogive (default=0.1*Lm)   
#' @param tmin minimum age
#' @param tmax maximum age 
#' @param h steepness of B&H SSR 
#' @param M natural mortality
#' @param Lc length at first capture
#' @param dLc width of ogive (default = 0.1*Lc)
#' @param Fvec vector of instantanious fishing fortality
#' @param R0 unfished recruitment default = 1
#' @return spm pars, equilibrium dynamics, life history table 
#' @export
#' @author Henning Winker (JRC-EC)

asem <- function(Loo,K,t0=-0.5,aW=0.01,bW=3.04,Lm,dLm=NULL,tmin=0,tmax,h,M,Lc=NULL,dLc=NULL,Fvec=NULL,R0=1){
  if(is.null(Fvec)) Fvec = c(0,exp(seq(-5,1.5,0.005)))
  tmax= round(tmax,0)
  age = tmin:tmax
  nages = length(age)
  ntilda  = mat.or.vec(length(Fvec),nages)  
  L = Lhalf = W = Whalf = mat.or.vec(nages,1)
  
  
  #Length-at-age
  L = Loo*(1-exp(-K*(age-t0)))
  Lhalf = Loo*(1-exp(-K*(age+0.5-t0)))
  # Weight-at-age
  W = aW*L^bW 
  Whalf = aW*Lhalf^bW 
  
  # Maturity-at-age
  if(is.null(dLm)) dLm = 0.1*Lm
  mat = 1/(1+exp(-(L-Lm)/dLm)) 
  
  #Selectivity-at-age (logistic)
  if(is.null(Lc)) Lc = Lm 
  if(is.null(dLc)) dLc = 0.1*Lc
  sel = 1/(1+exp(-(L-Lc)/dLc))
  
  
  # compute unFvecshed Spawning biomass per recruit (SBR0)
  n0 = 0
  for (t in 1:nages){
    if(t==1) n0[t] <- 1 
    if(t>1) n0[t] <- n0[t-1]*exp(-M) 
    if(t==nages) n0[t] <- n0[t]#/(1-exp(-M)) 
  }
  
  SBR0 = sum(n0*W*mat)
  EBR0 = sum(n0*W*sel) 
  
  # compute Spawning Biomass per rectuit as a function of F (SBR)
  for (t in 1:nages){
    if(t==1) ntilda[,t] <- 1
    if(t>1) ntilda[,t] <- ntilda[,t-1]*exp(-(M+sel[t-1]*Fvec)) 
    if(t==nages) ntilda[,t] <- ntilda[,t]#/(1-exp(-(M+sel[t]*Fvec))) 
  }
  
  ntilda[ntilda < 0] = 0 
  
  #get spawner biomass per recruit
  SBR= apply(ntilda %*% diag(W) %*% diag(mat),1,sum)
  # get Z matrix              
  Z = M+t(sel%*%t(Fvec))
  # Yield-per-recruit matrix
  YPR = apply(ntilda %*% diag(W) %*% diag(sel)*Fvec/Z*(1-exp(-Z)),1,sum)
  #get Exploitable Biomass per recruit
  EBR = apply(ntilda %*% diag(W) %*% diag(sel),1,sum)
  
  # equilibrium recruitment
  RF = R0*(4*h*SBR-(1-h)*SBR0)/(SBR*(5*h-1))
  recruits = RF 
  if(min(recruits) < 0) recruits[min(which(RF<0)):length(RF)] = 0
  
  # get spawner biomass depletion
  SBtoSB0 = SBR*recruits/SBR0
  EBtoEB0 = EBR*recruits/EBR0 #><>NEW
  
  lhtab = data.frame(L_a=L,W_a=W,Mat_a=mat,S_a=sel,n0=n0)
  
  eqf = data.frame(Fvec,bk=EBR*recruits/EBR0,yield=YPR*recruits,recruits=recruits, b=EBR*recruits)
  
  msy = max(eqf$yield)
  bmsyk = eqf$bk[eqf$yield==msy]
  fmsy = msy/eqf$b[eqf$yield==msy]
  # Find shape for  SBmsytoK 
  rshape = seq(0.1,5,0.001)
  check.shape =((rshape)^(-1/(rshape-1))-bmsyk)^2 #><> NEW (EBmsyEB0)
  shape = rshape[check.shape==min(check.shape)]
  r =fmsy*(shape-1)/(1-1/shape)
  spm = data.frame(r,shape,msy,fmsy,bmsyk)
  # return Spawning Biomass as fraction of unfished levels
  return(list(spm=spm,eqdyn=eqf,lhtab=lhtab))
}


  