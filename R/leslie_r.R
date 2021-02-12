#' leslie_r()
#'
#' Generate r from Leslie Matrix assuming BevHolt SSR (Afte MacAllister)
#'  
#' @param Loo  infinitive length (cm)  
#' @param K  brody  growth coefficient K   
#' @param t0 theoretical age a length zero (default -0.5)
#' @param aW parameter of the lenght-weight relation (default 0.01)
#' @param bW parameter of the lenght-weight relation (default 3.04)
#' @param mat length- c(Length,0) or age-at-50-maturity c(Age,1) (default is length, i.e. c(Lenght,0))    
#' @param dmat width of ogive (default=0.1*Mat)   
#' @param tmin minimum age
#' @param tmax maximum age 
#' @param M natural mortality
#' @param h steepness of BevHolt SSR 
#' @return r, G (generation time) 
#' @export
#' @author Henning Winker (JRC-EC)

leslie_r <- function(Loo=80,K=0.2,t0=-0.5,aW=0.01,bW=3.04,mat=c(35,0),dmat=NULL,minage=0,maxage=12,M=0.25,h=0.7){
  age = minage:maxage
  nages = length(age)
  
  #Length-at-age
  L = Loo*(1-exp(-K*(age-t0)))
  # Weight-at-age
  W = aW*L^bW 
  # Maturity
  if(is.null(dmat)){
    dmat = 0.1*mat[1]
  }
  if(mat[2]==0){  
  maturity = 1/(1+exp(-(L-mat[1])/dmat)) 
  } else {
  maturity = 1/(1+exp(-(age-mat[1])/dmat)) 
  }  
  
  # mean Spawner weight at age
  SWt = W*maturity
  
  # compute unfished Spawning biomass per recruit (SBR0)
  n0 = 0
  for (t in 1:nages)
  {
    if(t==1) n0[t] <- 1
    if(t>1) n0[t] <- n0[t-1]*exp(-M) 
    if(t==nages) n0[t] <- n0[t]/(1-exp(-M)) 
  }
  
  n0[n0 < 0] = 0
  # Unfished spawning biomass per recruit
  SBR0 = sum(n0*W*maturity)
  # Reproductive output Rs for bonyfish
  Rs = 4*h/(SBR0*(1-h))
  
  # Make Leslie matrix 
  L.Mat=mat.or.vec(nages,nages)
  
  L.Mat[1,] = Rs*SWt # for sharks e.g.  mat*N.pups_t 
  
  #fill rest of Matrix with Survival
  for(i  in 2:nages)
  {
    L.Mat[i,(i-1)] = exp(-M) 
  }
  
  # Net reproductive rate
  NR = sum(n0*SWt) 
  
  
  
  # return intrinsic rate of population increase r and generation GT
  return(list(r=log(as.numeric(eigen(L.Mat)$values[1])),GT = sum(age*n0*SWt)/NR))
  
}

