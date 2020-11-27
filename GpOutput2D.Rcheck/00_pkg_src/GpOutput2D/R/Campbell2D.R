#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# This script contains examples of functions with spatial output.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Example of function called Campbell2D
#'
#' @param X a matrix with eight columns, which correspond to the function input.
#' @param z1,z2 vectors with the spatial coordinates.
#'
#' @return a three-dimensional array with the output maps. The first two dimensions correspond to the maps dimension. The third dimension corresponds to the number of maps (or simulations).
#'
#'
#' @examples
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
#' @export
Campbell2D <- function(X,z1,z2){
  Zgrid <- expand.grid(z1=z1, z2=z2)

  n<-nrow(X) # number of design points
  nz1 <- length(z1); nz2 <- length(z2) # dimensions of the spatial output (maps, images... )

  Y <- lapply(1:n, function(i){
    Xi <- X[i,] # i-th design
    y1 <- Xi[1]*exp(-((0.8*Zgrid$z1+0.2*Zgrid$z2-10*Xi[2])**2)/(60*Xi[1]**2))
    y2 <- (Xi[2]+Xi[4])*exp(((0.5*Zgrid$z1+0.5*Zgrid$z2)*Xi[1])/500)
    y3 <- (Xi[5]*(Xi[3]-2))*exp(-(0.4*Zgrid$z1+0.6*Zgrid$z2-20*Xi[6])**2/(40*Xi[5]**2))
    y4 <- (Xi[6]+Xi[8])*exp(((0.3*Zgrid$z1+0.7*Zgrid$z2)*Xi[7])/250)
    return(matrix(y1+y2+y3+y4,nrow=nz1,ncol=nz2))
  })


  Ymaps<- array(unlist(Y),dim=c(nz1,nz2,n))

  return(Ymaps)
}


