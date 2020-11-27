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

    # object in matrix form

    y_mat<-matrix(y,ncol=n)
    coeff <- array(dim=c(K,L,n))

    i <- 0
    for(k in 1:K){
      for(l in 1:L){
        i<-i+1
        Phi_kl <- matrix(object[,k,,l],nrow=1)
        coeff[k,l,] <-  Phi_kl %*%y_mat
      }# end for l
    }# end for k



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

  coeffm <-matrix(coeff,ncol=n) # coefficients on matrix
  # K <-nrow(coeffm)# total number of coefficients
  # objectm <- array(object,dim=c((n1<-d[1]),(n2<-d[2]),K)) # object on 3D array

  # approximated data
  hatY <- array(dim=c((n1<-d[1]),(n2<-d[3]),n))

  t <- Sys.time()
  for(i in 1:n1){
    for(j in 1:n2){
      # hatY[i,j,]<-sapply(1:n,function(l){sum(coeffm[,l]*objectm[i,j,])})
      Phi_ij <- matrix(object[i,,j,],nrow=1)
      hatY[i,j,]<-Phi_ij%*%coeffm
    }# end for j
  }# end for i
  t <- Sys.time()-t

  return(hatY)

} # end Inverse2D
