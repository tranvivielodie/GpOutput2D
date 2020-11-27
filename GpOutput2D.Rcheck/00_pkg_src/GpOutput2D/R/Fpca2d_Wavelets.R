#' @title Two dimensional functional principal component analysis (FPCA)
#'        by using wavelet basis
#'
#' @description performs FPCA on two-dimensional functional data (data set of images, maps, ...).
#'              The implementation of FPCA is based on wavelet basis.
#'
#' @param x a three dimensional array, which contains two-dimensional data (images, maps, etc.).
#'          The two first dimensions correspond to data dimensions, which are denoted M and N. The third one is the size of the data set.
#' @param wf a character string which specifies the wavelet filter (see \code{\link{dwt.2d}}).
#' @param J depth of the wavelet decomposition, must be a number less than or equal to log(min(M,N),2).
#' @param boundary a character string which specifies how boundaries are treated. Only "periodic" is currently implemented (see \code{\link{dwt.2d}}).
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
#' @importFrom stats prcomp
#'
#' @export
Fpca2d.Wavelets<-function(x,
                         wf,J,boundary="periodic", # wavelet parameters
                         p=1,ncoeff=NULL, # number of coefficients for PCA
                         center=TRUE,scale.=FALSE,...){ # PCA parameters

  d<-dim(x) # dimensions of x

  res <- list(); class(res)<-"Fpca2d" # return scores, sdev etc.

  attr(res,"method")<-"wavelet" # method of Fpca2d (basis)
  attr(res,"wf")<-wf # wavelet type
  attr(res,"J")<-J # depth of decomposition
  attr(res,"boundary")<-boundary # boundary method for wavelets

  #############################
  ### Wavelet decomposition ###
  #############################
  w <- Wavelet2D(x,wf=wf,J=J,boundary=boundary)
  attr(res,"coeff")<-w
  attr(res,"type")<-attr(w,"type") # tensor type

  #################################
  ### Mean proportion of energy ###
  #################################


  SumW <- apply(w,3,sum)
  if(0%in% SumW){
    # simulations where there are only zero values
    idx.w <-which(SumW!=0) # indices where values are not all at zero
    w.poe <- w[,,idx.w]
    attr(w.poe,"type")<-attr(w,"type")
    attr(w.poe,"J")<-attr(w,"J")
    attr(w.poe,"wf")<-attr(w,"wf")
    attr(w.poe,"boundary")<-attr(w,"boundary")
    attr(w.poe,"dim")<-dim(w.poe)

    class(w.poe)<-"Wavelet2D"
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
  }# end if is.null(ncoeff)
  #-------------------------------------------------------------

  # coefficients for PCA
  idx_pca <- idx[1:ncoeff] # coefficient indices for PCA
  x_pca <- matrix(w,ncol=d[3])[idx_pca,]

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
  coeff_eig <-matrix(rep(0,nPC*d[1]*d[2]),ncol=nPC) # eigenfunctions
  coeff_eig[idx_pca,]<-rot
  coeff_eig <- array(coeff_eig,dim=c(d[1:2],nPC))

  ## coefficients on Wavelet2D object ##
  attr( coeff_eig,"type")<-attr(w,"type")
  attr( coeff_eig,"J")<-attr(w,"J")  # depth of wavelet decomposition
  attr( coeff_eig,"wf")<-attr(w,"wf") #wavelet type
  attr( coeff_eig,"boundary")<-attr(w,"boundary")
  class(coeff_eig)<-"Wavelet2D"

  ## eigenfunctions ##
  eigen_fct<-Inverse2D(coeff_eig)

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
}# end Fpca2d_wavelet
