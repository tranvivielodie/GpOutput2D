#' @title Inverse transform of wavelet or orthonormal B-splines basis.
#'
#' @description embeds the coefficients of "Wavelet2D" or "OrthonormalBsplines2D" basis
#' onto two-dimensional domain.
#' @author Tran Vi-vi Elodie PERRIN
#'
#' @param object an object of class \code{\link{Wavelet2D}} or \code{\link{OrthoNormalBsplines2D}}.
#' @param ... if \code{object} is an object of class \code{OrthoNormalBsplines2D},
#' the matrix of control points must be given (see \code{\link{Inverse2D.OrthoNormalBsplines2D}}).
#'
#' @return a three dimensional array.
#' The first two dimensions correspond to maps dimensions.
#' The third dimension corresponds to the size of the data set.
#'
#' @seealso \code{\link{Inverse2D.OrthoNormalBsplines2D}} \code{\link{Inverse2D.Wavelet2D}}
#' @export
Inverse2D <- function(object,...){
  UseMethod("Inverse2D",object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Mean proportion of energy
#'
#' @description Wrapper function to compute mean proportion of energy for given coefficents.
#' The valid objects are "Wavelet2D" and "CoefOrthonormalBsplines2D".
#'
#' @param object an object of class \code{\link{Wavelet2D}} or \code{CoefOrthonormalBsplines2D}
#' which correspond to coefficients from wavelets or Bsplines basis.
#'
#' @return a matrix K×L which contains the mean proportion of energy of each coefficient,
#'         with K×L is the dimensions of the functional basis.
#'
#'
#'
#' @author Tran Vi-vi Elodie PERRIN
#'
#' @keywords internal
#' @export
MeanPoe<-function(object){
  UseMethod("MeanPoe",object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Two-dimensional Functional Principal Component Analysis (FPCA).
#'
#' @description A wrapper function to perform FPCA on two-dimensional data (images, maps, etc.),
#'              given a projection method. For instance, the valid projections are
#'              two methods based on projection method (here, "Wavelets" and "Bsplines").
#'
#' @param method a character string which specifies the method used to implement FPCA.
#' @param ... further parameters for the projection method: see \code{\link{Fpca2d.Wavelets}}, for using wavelets, and
#' \code{Fpca2d.Bsplines}, for using orthonormal B-splines.
#'
#' @return an object of class "Fpca2d" which is a list with the following items :
#'  \itemize{
#'    \item \code{sdev} : the standard deviations of the principal components
#'                        (i.e., the square roots of the eigenvalues of the covariance operator).
#'    \item \code{EigFct} : a three dimensional array with the eigenfunctions.
#'                          The two first dimensions correspond to data (images, maps, etc.) dimensions.
#'                          The third one is the number of modeled principal component.
#'    \item \code{x} : if \code{retx} is \code{TRUE} (see \code{\link{prcomp}}).
#'                     Data scores (coordinates in the eigen basis) are returned.
#'    \item \code{center}, \code{scale} : a logical value which indicates if data coefficients are respectively centered or scaled.
#'  }
#'
#' @author  Tran Vi-vi Elodie PERRIN
#'
#' @import lhs
#' @import DiceDesign
#'
#' @seealso \code{\link{Fpca2d.Wavelets}}
#'
#' @examples
#'
#' ################################
#' ### two-dimensional data set ###
#' ################################
#' n<-200 # size of the learning sample
#' nz<-64; z <- seq(-90,90,length=nz) # spatial domain
#'
#' ### inputs of Campbell2D ###
#' library(lhs)
#' library(DiceDesign)
#'
#' x <- maximinLHS(n=n,k=8)
#' X <-maximinSA_LHS(x)$design
#' X<-X*6 -1
#'
#' # Campbell2D
#' Y <- Campbell2D(X,z,z)
#'
#' ############
#' ### FPCA ###
#' ############
#'
#'
#' ### using wavelet basis ###
#' fpca_w<- Fpca2d(Y,method="Wavelets",
#'                 wf="d4", J=1, # wavelet parameters
#'                 ncoeff=1200, rank.=5) # FPCA configuration
#'
#' plot(fpca_w,type="eigenfunctions")
#' plot(fpca_w)
#'
#' \dontrun{
#' ### using B-splines basis ###
#'
#' # knots for B-splines basis
#' K<-35
#' z.knots <- seq(-90,90,length=K)
#'
#' fpca_Bs<- Fpca2d(Y,method="Bsplines",
#'                 z1=z,z2=z,z1.knots=z.knots,z2.knots=z.knots, # wavelet parameters
#'                 ncoeff=1225, rank.=2) # FPCA configuration
#'
#' plot(fpca_Bs,type="eigenfunctions")
#' plot(fpca_Bs)
#' }
#'
#' @export
Fpca2d <- function(method = c("Bsplines", "Wavelets"), ...) {
  basis <- match.arg(method)
  fun <- get(paste("Fpca2d.", basis, sep = ""))
  model <- fun(...)
  return(model)
}
