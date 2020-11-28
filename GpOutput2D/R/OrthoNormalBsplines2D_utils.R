#' @title Get control points of two-dimensional B-splines basis
#'
#' @param object an object of class \code{\link{OrthoNormalBsplines2D}}
#' with functions of the B-splines basis.
#' @param y two-dimensional data to approximate on B-splines basis.
#' @param ... other arguments.
#'
#' @return a matrix with control points on two-dimensional B-splines basis.
#'
#' @method coef OrthoNormalBsplines2D
#'
#' @seealso \code{\link{OrthoNormalBsplines2D}} \code{\link{Inverse2D.OrthoNormalBsplines2D}}
#' \code{\link{Inverse2D}}
#'
#' @export
coef.OrthoNormalBsplines2D<-function(object, y,...){

  # basis dimensions
  d<-dim(object)

  # data set size
  dy<-dim(y)
  n1 <- dy[1]; n2 <- dy[2]; n<-dy[3]

  # basis dimensions
  K<-d[2]; L<-d[4]

  # y in matrix form
  y_mat<-matrix(y,ncol=n)

  ###############
  # coefficients
  ###############

  # In matrix form
  basis <- attr(object,"matrix") # basis functions
  coeff_mat <- t(basis)%*%y_mat

  # Array
  coeff <- array(coeff_mat,dim=c(K,L,n))
  attr(coeff,"matrix")<-coeff_mat


  class(coeff)<-"coef_OrthoNormalBsplines2D"

  return(coeff)

}# end coef



##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
#' @title Inverse transform of control points of two-dimensional B-splines basis
#'
#' @param object an object of class \code{\link{OrthoNormalBsplines2D}}
#' with functions of the B-splines basis.
#' @param coeff  matrix with control points on two-dimensional B-splines basis.
#' @param ... other arguments.
#'
#' @return an array with two-dimensional approximated data.
#'
#' @method Inverse2D OrthoNormalBsplines2D
#'
#' @seealso \code{\link{coef.OrthoNormalBsplines2D}} \code{\link{OrthoNormalBsplines2D}}
#' \code{\link{Inverse2D}}
#'
#' @export
Inverse2D.OrthoNormalBsplines2D<-function(object,coeff,...){
  # basis dimensions
  d<-dim(object)

  # data set size
  n<-dim(coeff)[3]

  coeffm <-attr(coeff,"matrix") # coefficients on matrix

  # objects in matrix form
  basis <- attr(object,"matrix") # basis functions
  Y <-basis%*%coeffm

  # approximated data
  hatY <- array(Y,dim=c(d[1],d[3],n))

  return(hatY)

} # end Inverse2D
