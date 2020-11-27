#' @title Mean proportion of energy of wavelet coefficients
#'
#' @param object an object of class \code{\link{Wavelet2D}}.
#'
#' @return a matrix K×L which contains the mean proportion of energy of each coefficient,
#'         with K×L is the dimensions of the wavelet basis
#'
#' @method MeanPoe Wavelet2D
#'
#' @keywords internal
#' @export
MeanPoe.Wavelet2D <- function(object){
  beta2 <- object**2 # coefficients square
  SumBeta2 <- apply(beta2, 3, sum) # denominator (total energy of each map)

  # compute energy of wavelet coefficients for each image/map
  ncoeff <- dim(beta2)[1:2] # number of coefficients
  poe <- matrix(nrow=ncoeff[1],ncol=ncoeff[2]) # mean proportion of energy
  for(k in 1:ncoeff[1]){
    for(l in 1:ncoeff[2]){
      poe[k,l]<- mean(beta2[k,l,]/SumBeta2)
    }# end for l
  }# end for k


  attr(poe,"type")<-attr(object,"type")
  attr(poe,"J")<-attr(object,"J") # depth of the wavelet decomposition
  attr(poe,"wf")<-attr(object,"wf") # wavelet type
  attr(poe,"boundary")<-attr(object,"boundary")

  class(poe)<-"MeanPoe.Wavelet2D"
  return(poe)
}# end MeanPoe.Wavelet2D



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title Mean proportion of energy of B-splines control points
#'
#' @param object control point outputs of  \code{\link{coef.OrthoNormalBsplines2D}} .
#'
#' @return a matrix K×L which contains the mean proportion of energy of each coefficient,
#'         with K×L is the dimensions of the B-splines basis
#'
#' @method MeanPoe coef_OrthoNormalBsplines2D
#'
#' @keywords internal
#' @export
MeanPoe.coef_OrthoNormalBsplines2D <- function(object){
  beta2 <- object**2 # coefficients square
  SumBeta2 <- apply(beta2, 3, sum) # denominator (total energy of each map)

  # compute energy of wavelet coefficients for each image/map
  ncoeff <- dim(beta2)[1:2] # number of coefficients
  poe <- matrix(nrow=ncoeff[1],ncol=ncoeff[2]) # mean proportion of energy
  for(k in 1:ncoeff[1]){
    for(l in 1:ncoeff[2]){
      poe[k,l]<- mean(beta2[k,l,]/SumBeta2)
    }# end for l
  }# end for k

  class(poe)<-"MeanPoe.coef_OrthoNormalBsplines2D"
  return(poe)
}# end MeanPoe.Wavelet2D
