#' from_logith
#'
#' convert steepness from logit
#' @param logit_h logit converted steepness h
#' @return steepness h 
#' @export
from_logith <- function(logit_h){
  out=  0.2001 + 0.7998*1/(1+exp(-logit_h))
  out}

#' to_logith
#'
#' convert steepness to logit
#' @param logit_h logit converted steepness h
#' @return logit transformed steepness h 
#' @export
to_logith <- function(h){
  x = seq(-5,5,0.0001)
  h_i = 0.2001 + 0.7998*1/(1+exp(-x))
  out=mean(x[(h_i-h)^2==min((h_i-h)^2)])
  out}
