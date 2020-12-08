#' @title Creating a two-dimensional othonormal B-splines basis
#'
#' @param x,y two vectors of argument values at which the B-spline basis functions are to be evaluated.
#' @param x.knots,y.knots The full set of knots used to define the basis functions.
#' @param order Order of the spline fit (degree= order-1). The default is 2.
#' @param ortho a character string which specify the orthogonalization method.
#' If it is "Redd", the \code{\link{orthogonalsplinebasis}} package is used.
#' If it is "GS", the Gram-Schmidt is used. The default is "Redd".
#' @param expand_knots a boolean. If it is TRUE, knots are expanded
#' for appropriate number of knots in bsplines (see \code{\link{expand.knots}}). The default is TRUE.
#' @param ... other arguments (see \code{\link{SplineBasis}}).
#'
#' @return an array with two dimensional orthonormal B-splines.
#' The two first dimensions correspond to the data dimensions.
#' The last two dimensions correspond to the basis dimensions.
#'
#' @import orthogonalsplinebasis
#' @importFrom pracma gramSchmidt
#'
#' @seealso \code{\link{coef.OrthoNormalBsplines2D}} \code{\link{Inverse2D.OrthoNormalBsplines2D}}
#' \code{\link{Inverse2D}}
#'
#' @examples
#'
#' #####################################
#' #
#' #  Create two-dimensional data set
#' #
#' ####################################
#'
#' # inputs of the Campbell2D function
#' x1<-rep(-1,8); x2<-rep(5,8); x3<-c(5,3,1,-1,5,3,1,-1)
#' X <- rbind(x1,x2,x3)
#'
#' # spatial domain of the Campbell2D output
#' nz<-64 # root of the size of the spatial domain
#' z<-seq(-90,90,length=nz)
#'
#' # Campbell2D function
#' Y = Campbell2D(X,z,z)
#'
#' ###########################################
#' #
#' #  Create two-dimensional B-splines basis
#' #
#' ###########################################
#'
#' # knots for B-splines basis
#' K<-35
#' z.knots <- seq(-90,90,length=K)
#'
#' # Generating a two-dimensional othonormal B-splines basis
#' OPhi <-OrthoNormalBsplines2D(z,z,z.knots,z.knots,ortho="GS")
#'
#' ########################
#' #
#' #  Get control points
#' #
#' ########################
#'
#' coeff_Y<-coef(OPhi,Y)
#'
#' #######################
#' #
#' #  Approximation of Y
#' #
#' ########################
#'
#' hatY <-Inverse2D(OPhi,coeff_Y)
#'
#'
#' @export
OrthoNormalBsplines2D <-function(x,y,x.knots,y.knots,order=2,ortho="GS",expand_knots=TRUE,...){

  # Dimensions
  nx<-length(x); ny<-length(y) # data dimensions

  # Expands knots for appropriate number of knots
  if(isTRUE(expand_knots)){
    x.knots <-expand.knots(x.knots,order=order)
    y.knots <-expand.knots(y.knots,order=order)
  }

  if(!ortho %in% c("Redd","GS")){
    stop("ortho must be 'Redd' or 'GS'.")
  }# end if

  # orthogonolized B-splines
  if(ortho=="Redd"){
    # orthogonolized B-splines basis on x and y
    OPhi1 <-OBasis(x.knots, order=order,...) # x
    OPhi2 <-OBasis(y.knots, order=order,...) # y

    # Orthogonal B-splines on x and y
    OPhi1 <- evaluate(OPhi1,x)
    OPhi2 <- evaluate(OPhi2,y)

    OPhi1[which(is.na(OPhi1))]<-0
    OPhi2[which(is.na(OPhi2))]<-0
  }

  if(ortho=="GS"){
    # B-splines basis on x and y
    Phi1 <-SplineBasis(x.knots, order=order,...) # x
    Phi2 <-SplineBasis(y.knots, order=order,...) # y

    # B-splines on x and y
    Phi1 <- evaluate(Phi1,x)
    Phi2 <- evaluate(Phi2,y)

    Phi1[which(is.na(Phi1))]<-0
    Phi2[which(is.na(Phi2))]<-0

    # Gram-Schmidt orthogonalization
    OPhi1 <-gramSchmidt(Phi1)$Q
    OPhi2 <-gramSchmidt(Phi2)$Q

    # memory deallocation
    rm(list=c("Phi1","Phi2"))
  }


  # basis dimensions
  K <- ncol(OPhi1); L <-ncol(OPhi2)

  # normalization
  norm_eucl1<- sqrt(apply(OPhi1**2,2,sum)) # norm
  norm_eucl2<- sqrt(apply(OPhi2**2,2,sum)) # norm

  OPhi1 <- sapply(1:K, function(k){OPhi1[,k]/norm_eucl1[k]})
  OPhi2 <- sapply(1:L, function(l){OPhi2[,l]/norm_eucl2[l]})

  # memory deallocation
  rm(list=c("norm_eucl1","norm_eucl2"))

  # tensor product
  OPhi_2D <-OPhi1%o%OPhi2

  # memory deallocation
  rm(list=c("OPhi1","OPhi2"))

  ## matrix (to compute coefficients)
  new.OPhi <- aperm(OPhi_2D, c(1,3,2,4))
  Phi_mat <-matrix(new.OPhi,nrow=nx*ny,ncol=K*L)
  attr(OPhi_2D,"matrix")<-Phi_mat

  class(OPhi_2D)<-"OrthoNormalBsplines2D"

  return(OPhi_2D)
}# end OrthoNormalBsplines2D
