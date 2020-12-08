#' @title Two dimensional functional principal component analysis (FPCA)
#' by using B-splines basis
#'
#' @param x a three dimensional array, which contains two-dimensional data
#' (images, maps, etc.). The two first dimensions correspond to data dimensions,
#' which are denoted M and N. The third one is the size of the data set.
#' @param z1,z2 two vectors of argument values at which the B-spline basis functions are to be evaluated.
#' @param z1.knots,z2.knots The full set of knots used to define the basis functions.
#' @param norder Order of the spline fit (degree= order-1). The default is 2.
#' @param ortho a character string which specify the orthogonalization method.
#' If it is "Redd", the \code{\link{orthogonalsplinebasis}} package is used.
#' If it is "GS", the Gram-Schmidt is used. The default is "GS".
#' @param expand_knots a boolean. If it is TRUE, knots are expanded
#' for appropriate number of knots in bsplines (see \code{\link{expand.knots}}). The default is FALSE.
#' @param p a value which fixes the total mean proportion of energy (or mean spatial variance).
#'          The number of coefficients (\code{ncoeff}) used for PCA is calibrated according to its value.
#'          The default is 1. If a value is given in \code{ncoeff}, \code{p} is not used.
#' @param ncoeff a value which fixes the number of coefficients used for PCA.
#'               The default is \code{NULL}. If it is \code{NULL}, the number of coefficients is calibrated by using the parameter \code{p}.
#' @param center a logical value indicating whether the coefficients should be shifted to be zero centered. The default is \code{TRUE}.
#' @param scale. a logical value indicating whether the coefficients should be scaled to have unit variance before the analysis takes place.
#'              The default is \code{FALSE}.
#' @param ... arguments of \code{\link{prcomp}}, which can fix characteristics of PCA.
#'
#'
#' @return an object of class "Fpca2d" which is a list with the following items :
#'  \itemize{
#'    \item \code{sdev} : the standard deviations of the principal components
#'                        (i.e., the square roots of the eigenvalues of the covariance operator).
#'    \item \code{EigFct} : a three dimensional array with the eigenfunctions.
#'                          The two first dimensions correspond to data (images, maps, etc.) dimensions.
#'                          The third one is the maximal number of principal component used.
#'    \item \code{x} : if \code{retx} is \code{TRUE} (see \code{\link{prcomp}}).
#'                     Scores (coordinates in the eigen basis) of the data are returned.
#'    \item \code{center}, \code{scale} : a logical value which indicates if data coefficients are respectively centered or scaled.
#'  }
#'
#' @seealso \code{\link{Fpca2d}}
#'
#' @importFrom stats prcomp coef
#'
#' @export
Fpca2d.Bsplines<-function(x,
                          z1,z2,z1.knots,z2.knots,norder=2, ortho="GS",expand_knots=FALSE,# B-splines  parameters
                          p=1,ncoeff=NULL, # number of coefficients for PCA
                          center=TRUE,scale.=FALSE,...){

  d<-dim(x) # dimensions of x

  res <- list(); class(res)<-"Fpca2d" # return scores, sdev etc.

  attr(res,"method")<-"Bsplines" # method of Fpca2d (basis)
  attr(res,"knots")<-data.frame(x.knots=z1.knots,y.knots=z2.knots) # knots
  attr(res,"coordinates")<-data.frame(x=z1,y=z2) # knots
  attr(res,"order")<-order # B-splines order

  ###############################
  ### B-splines decomposition ###
  ###############################
  OPhi <- OrthoNormalBsplines2D(x=z1,y=z2,x.knots=z1.knots,y.knots=z2.knots,
                                order=norder,ortho=ortho, expand_knots=expand_knots)
  attr(res,"SplinesBasis")<-OPhi

  # coefficients
  w<-coef(OPhi,x)
  attr(res,"coeff")<-w

  #################################
  ### Mean proportion of energy ###
  #################################
  SumW <- apply(w,3,sum)
  if(0%in% SumW){
    # simulations where there are only zero values
    idx.w <-which(SumW!=0) # indices where values are not all at zero
    w.poe <- w[,,idx.w]

    class(w.poe)<-"coef_OrthoNormalBsplines2D"

    poe <- MeanPoe(w.poe)
  }else{
    poe <- MeanPoe(w)
  }
  attr(res,"mean_poe")<-poe
  poe<-as.vector(poe)

  #######################################
  ### poe indices in decreasing order ###
  #######################################
  idx <-order(poe,decreasing = TRUE)

  ############################
  ### coefficients for PCA ###
  ############################

  #------------------ if is.null(ncoeff) -----------------------
  # if ncoeff is not given, then determinate ncoeff value
  if(is.null(ncoeff)){
    CumPoe <- cumsum(poe[idx]) # cumsum poe

    # if p==1 (<=> total energy)
    if(p==1){
      # choose the minimum number of coefficients such as CumPoe=1
      ncoeff <- min(which(CumPoe==1))
    }else{
      # choose the number of coefficients such as CumPoe<=p
      ncoeff <- length(which(CumPoe<=p))
    }# end ifelse (p==1)
  }else{
    # coefficients for PCA

    if(length(poe)<ncoeff){
      warning("The size of the B-splines basis is less than ncoeff. ncoeff has been replaced by
            the total size of the projection.")
      ncoeff <-length(poe)
    }
    if(length(poe)>ncoeff){
      warning("The size of the B-splines basis is more than ncoeff. ncoeff has been replaced by
            the total size of the projection.")
      ncoeff <-length(poe)
    }

  }# end if is.null(ncoeff)
  #-------------------------------------------------------------





  idx_pca <- idx[1:ncoeff] # coefficient indices for PCA
  w_mat <-matrix(w,ncol=d[3])
  x_pca <- w_mat[idx_pca,]

  ###########
  ### PCA ###
  ###########
  pca <- prcomp(t(x_pca),center=center,scale.=scale.,...)

  ######################
  ### eigenfunctions ###
  ######################
  rot <- pca$rotation # rotation matrix
  nPC <- ncol(rot) # number of principal components

  ## wavelet coefficients of the eigenfunctions ##
  d.basis <-dim(OPhi)
  coeff_eig_mat <-matrix(rep(0,nPC*d.basis[2]*d.basis[4]),ncol=nPC) # eigenfunctions
  coeff_eig_mat[idx_pca,]<-rot
  coeff_eig <- array(coeff_eig_mat,dim=c(d.basis[c(2,4)],nPC))

  attr(coeff_eig,"matrix")<-coeff_eig_mat

  ## eigenfunctions ##
  eigen_fct<-Inverse2D(OPhi,coeff_eig)

  ##############
  ### return ###
  ##############
  res$sdev<-pca$sdev
  res$EigFct<- eigen_fct # eigenfunctions
  res$x <- pca$x # scores
  res$center <-center
  res$scale <-pca$scale

  # informations of the pca step
  attr(res,"pca")<-pca
  attr(res,"ncoeff")<-ncoeff # number of coefficients
  p.res <-sum(poe[idx_pca])
  attr(res,"total_poe")<-p.res

  return(res)

}# end Fpca2d.Bsplines

