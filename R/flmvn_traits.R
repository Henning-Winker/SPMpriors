#' flmvn_traits()
#'
#' Generate stock parameters with tuning from FishLife (Thorson et al. 2000) using MVN
#'
#' @param Genus genus name must match FishBase 
#' @param Species species name must match FishBase 
#' @param Loo optional tuning prior dist for Linf c(mu, cv)  
#' @param K optional tuning prior dist for brody K c(mu, cv)  
#' @param tmax optional tuning prior dist for maximum age c(mu, cv)  
#' @param tm optional tuning prior dist for age maturity c(mu, cv)  
#' @param Lm  optional tuning prior dist for length maturity c(mu, cv)  
#' @param h optional uniform tuning prior dist for steepness h c(min, max)  
#' @param nmc number Monte-Carlo MVN draws
#' @param upper.quant upper quatile for CIs ,
#' @param Plot option TRUE/FALSE to switch of plots
#' @param savepng option TRUE/FALSE to save plots as png 
#' @param PlotPath directory to save plot if savepng=TRUE
#' @return stock parameter list, mvn covar, Monte-Carlo matrix of stock pars 
#' @export
#' @author Henning Winker (JRC-EC)

flmvn_traits <- function(Genus="Rhabdosargus",Species="globiceps",Loo = NULL,K=NULL,tmax=NULL,tm=NULL,M=NULL,Lm=NULL,h=NULL,nmc = 2*10^5,upper.quant=0.9,Plot=TRUE,savepng=FALSE,PlotPath=getwd()){
  parms=c("Loo","K","Lm","tm","tmax","M","logitbound_h","ln_margsd","rho","ln_r","ln_G")
  
  if(!Plot){
    PlotPath=tempdir()
    savepng = TRUE
  } 
  # Search Taxa
  taxa = FishLife::Search_species(Genus=Genus,Species = Species,add_ancestors=TRUE)$match_taxonomy
  tax = strsplit(taxa[[1]], "_")
  # Predict LifeLife Traits  
  
  if(savepng){
  png(file = paste0(PlotPath,"/",Genus,".",Species,".fl_ellipse.png"), width = 7, height = 8, 
      res = 200, units = "in")
  }
  predfl =FishLife::Plot_taxa(taxa,mfrow=c(3,2))
  if(savepng) dev.off()
  
  mu= predfl[[1]]$Mean_pred
  covs = predfl[[1]]$Cov_pred
  
  priors = parms
  priors[parms=="logitbound_h"] = "h"
  priors[parms=="ln_r"] = "r"
  priors[parms=="ln_G"] = "G"
  priors[parms=="ln_margsd"] = "sigR"
  
  pars= c(priors)
  priors = priors[c(1:6)]
  mvn = mvtnorm::rmvnorm(nmc,mean =  mu,sigma = covs,  method=c("eigen"))
  mcparms = data.frame(mvn)
  
  nprior = 1
  prd = NULL
  for(i in 1:6){
    getpr = get(priors[i])  
    
    # pdf of prior
    if(is.null(getpr)==FALSE){
      nprior = nprior+1 
      prd = cbind(prd,dnorm(mcparms[,parms[i]],log(getpr[1]),getpr[2]))
    }
  } # end of loop
  
  if(is.null(prd)==FALSE){
    prp = apply(prd,1,sum)/max(apply(prd,1,sum)) # Convert to probability  
    # subsample conditioned given prp
    rand = runif(nrow(mcparms))
    sel = which(rand < prp)
    subparms = mcparms[sel,]
  } else {
    subparms = mcparms
  }
  # subset h range
  if(is.null(h)){
    if(tax[[1]][1]=="Elasmobranchii"){h = c(0.2,0.6)} else {h=c(0.3,0.9)}
  }
  subparms = subparms[subparms$h>h[1] & subparms$h<h[2],]
  # Get means and vars from subset
  mu.upd = apply(subparms,2,mean)
  covs.upd = covs 
  diag(covs.upd) = apply(subparms,2,var)
  updfl = predfl
  updfl[[1]]$Mean_pred = mu.upd # Here only updating on species (stock level) 
  updfl[[1]]$Cov_pred = covs.upd # Here only updating on species (stock level)
  
  # Plot
  if(savepng){
  png(file = paste0(PlotPath,"/",Genus,".",Species,".fl_stktraits.png"), width = 7, height = 8, 
      res = 200, units = "in")}
  
  Par = list(mfrow=c(4,3),mai=c(0.5,0.5,0,.1),omi = c(0.2,0.2,0.1,0) + 0.1,mgp=c(2,0.5,0), tck = -0.02,cex=0.7)
  par(Par)
  for(i in 1:11){
    if(pars[i]!="h" | pars[i]!="rho"){
      xi = sort(rlnorm(1000,mu[parms[i]],sqrt(diag(covs)[parms[i]])))  
      pdf.upd =dlnorm(xi,mu.upd[parms[i]],sqrt(diag(covs.upd)[parms[i]]))
      pdf.fl = dlnorm(xi,mu[parms[i]],sqrt(diag(covs)[parms[i]])) 
    } 
    if(pars[i]=="h") {
      xi = sort((rnorm(1000,mu[parms[i]],sqrt(diag(covs)[parms[i]]))))  
      pdf.upd =dnorm(xi,mu.upd[parms[i]],sqrt(diag(covs.upd)[parms[i]]))
      pdf.fl = dnorm(xi,mu[parms[i]],sqrt(diag(covs)[parms[i]]))  
      xi = from_logith(xi)
    } 
    if(pars[i]=="rho"){
      xi = sort(rnorm(1000,mu[parms[i]],sqrt(diag(covs)[parms[i]])))  
      pdf.upd =dnorm(xi,mu.upd[parms[i]],sqrt(diag(covs.upd)[parms[i]]))
      pdf.fl = dnorm(xi,mu[parms[i]],sqrt(diag(covs)[parms[i]]))
    }
    if(i==12){
      
    }  
    
    plot(0,0,type="h",xlab=paste(pars[i]),ylab="Density",xlim=quantile(xi,c(0.01,0.99)),ylim=range(c(pdf.fl,pdf.upd),c(0.01)))
    polygon(c(xi,rev(xi)),c(pdf.fl,rep(0,1000)),col="grey")  
    polygon(c(xi,rev(xi)),c(pdf.upd,rep(0,1000)),col=rgb(1,0,0,0.5),border=rgb(1,0,0,0.5))  
    
  }
  plot(0,0,type="n",axes=F,ylab="",xlab="")
  legend("topright",c("FishLife","Stock"),pch=22,pt.cex = 2,pt.bg=c("grey",rgb(1,0,0,0.5)),bty="n",cex=1.1)
  if(savepng) dev.off()
  
  
  quant = qt(upper.quant,1000)
  # Summarize results
  trait.upd = subparms[,paste(parms)]
  trait.upd[,c(1:6,8,10:11)] = exp(trait.upd[,c(1:6,8,10:11)])
  trait.upd[,7] = from_logith(trait.upd[,7])
  CV.upd=apply(trait.upd,2,sd)/apply(trait.upd,2,mean)
  
  trait.fl = mcparms[,paste(parms)]
  trait.fl[,c(1:6,8,10:11)] = exp(trait.fl[,c(1:6,8,10:11)])
  trait.fl[,7] = from_logith(trait.fl[,7])
  CV.fl=apply(trait.fl,2,sd)/apply(trait.fl,2,mean)
  
  sd.fl = sqrt(diag(covs[parms,parms]))
  mean.fl= mu[parms]
  mu.fl = ifelse(parms!="logitbound_h",exp(mean.fl),from_logith(mean.fl))
  
  sd.upd = sqrt(diag(covs.upd[parms,parms]))
  mean.upd= mu.upd[parms]
  mu.upd = ifelse(parms!="logitbound_h",exp(mean.upd),from_logith(mean.upd))
  lcl.upd = ifelse(parms!="logitbound_h",exp(mean.upd-quant*sd.fl),from_logith( mean.upd-quant*sd.upd))
  ucl.upd = ifelse(parms!="logitbound_h",exp(mean.upd+quant*sd.fl),from_logith(mean.upd+quant*sd.upd))
  out = data.frame(trait=pars,mu.fl,CV.fl,mu.upd,CV.upd,lcl.upd,ucl.upd,upper.quant=rep(upper.quant,length(mu.fl)))
  out[,-1] = round(out[,-1],4)
  colnames(out) = c("trait","mu.sp","cv.sp","mu.stk","cv.stk","lc.stk","uc.stk","upper.quant")
  parname = names(trait.upd)
  parname[c(7,8,10,11)] = c("h","sigR","r","G") 
  colnames(trait.upd) = parname
  rownames(out) <- 1:nrow(out)
  return(list(traits=out,vcm=updfl,mvnstk = trait.upd))
  
} # End of function 

