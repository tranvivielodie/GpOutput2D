#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title Wavelet coefficients in matrix form from a \code{dwt.2d} object
#'
#' @description organizes the wavelet coefficients from a 2D DWT into a single matrix.
#' @author Tran Vi-vi Elodie PERRIN
#'
#' @param x an object of class \code{\link{dwt.2d}}.
#' @param nrows,ncols the number of rows and columns.
#' @param ... additional arguments to be passed to or from methods.
#'
#' @return an object of class \code{matrix.dwt.2d} corresponding to a matrix with wavelet coefficients.
#'
#' @method as.matrix dwt.2d
#'
#' @keywords internal
#' @export
as.matrix.dwt.2d <- function(x,nrows,ncols,...){
  d<-c(nrows,ncols) # matrix dimensions
  MatW <- matrix(nrow=d[1],ncol=d[2])  # matrix with wavelet coefficients
  NamesJ <- matrix(rep("",d[1]*d[2]),ncol=d[2]) # function names (tensor product)

  J <- attr(x,"J")# depth of the decomposition
  for(j in 1:J){
    dj <-d/(2**j) # map dimsensions at level j

    # empty matrix for level j
    d0 <- dj*2 # dimension
    MatWj <-matrix(nrow=d0[1],ncol=d0[2])
    NamesWj <- matrix(rep("",d0[1]*d0[2]),ncol=d0[2])  # function names (tensor product)

    # HL
    NamesW <-  paste("HL",j,sep="") # wavelet function names
    wHL <-x[[NamesW]] # wavelet function
    idx1 <- 1:(d0[1]-dj[1]) # rows
    idx2<-(d0[2]-dj[2]+1):d0[2] # columns
    MatWj[idx1,idx2]<-wHL
    NamesWj[idx1,idx2] <- NamesW

    # HH
    NamesW <-  paste("HH",j,sep="") # wavelet function names
    wHH <-x[[NamesW]] # wavelet function
    idx1 <-(d0[1]-dj[1]+1):d0[1]# rows
    idx2<-(d0[2]-dj[2]+1):d0[2] # columns
    MatWj[idx1,idx2]<-wHH
    NamesWj[idx1,idx2] <- NamesW

    # LH
    NamesW <-  paste("LH",j,sep="") # wavelet function names
    wLH <-x[[NamesW]] # wavelet function
    idx1 <-(d0[1]-dj[1]+1):d0[1]# rows
    idx2<-1:(d0[2]-dj[2]) # columns
    MatWj[idx1,idx2]<-wLH
    NamesWj[idx1,idx2] <- NamesW


    # Output matrix
    if(j==1){
      MatW<-MatWj
      NamesJ <- NamesWj
    }else{
      idx1<-1:d0[1]
      idx2<-1:d0[2]
      MatW[idx1,idx2]<-MatWj
      NamesJ[idx1,idx2] <- NamesWj
    }# end ifelse

  }# end for j

  # scale function
  NameScale <- paste("LL",J,sep="")
  wScale <-x[[NameScale]] # scale function
  dJ <-dim(wScale)# dim scale function
  MatW[1:dJ[1],1:dJ[2]]<-wScale
  NamesJ[1:dJ[1],1:dJ[2]]<-NameScale

  attr(MatW,"type")<-NamesJ
  attr(MatW,"J")<-J # depth of the wavelet decomposition
  attr(MatW,"wf")<-attr(x,"wavelet") # wavelet type
  attr(MatW,"boundary")<-attr(x,"boundary")
  attr(MatW,"dim")<-c(nrows,ncols) # maps dimension
  class(MatW)<- "matrix.dwt.2d "
  return(MatW)
} # end as.matrix.dwt.2d

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title \code{dwt.2d} object by using a matrix of wavelet coefficients
#'
#' @description creates an object of class \code{\link{dwt.2d}} by using
#'              a matrix of wavelet coefficients.
#'
#' @author Tran Vi-vi Elodie PERRIN
#'
#' @param x a object of class "matrix.dwt.2d ", which is a matrix which contains wavelet coefficients.
#'
#' @importFrom waveslim dwt.2d
#'
#' @return an object of class \code{\link{dwt.2d}}.
#'
#' @keywords internal
#' @export
as.dwt.2d <- function(x){
  # matrix with the associate wavelet tensorisation
  typ <- attr(x,"type")
  # wavelet type
  wf<-attr(x,"wf")
  # depth of the wavlet decomposition
  J<-attr(x,'J')
  # Boundary method
  boundary<-attr(x,'boundary')

  # "empty" dwt.2d object
  d <- attr(x,"dim") # map dimension
  zero <- matrix(rep(0,d[1]*d[2]),nrow=d[1],ncol=d[2]) # zero matrix
  res <- dwt.2d(zero,wf=wf,J=J,boundary=boundary)

 # asign wavelet coefficients on the "empty" dwt.2d object
 NamesW <- names(res) # names of the wavelet tensor

 for(nm in NamesW){
    wnm <-res[[nm]]
    dnm <- dim(wnm)
    res[[nm]]<- matrix(x[which(typ==nm)],nrow=dnm[1],ncol=dnm[2])
 }# end for nm

  return(res)
}# end as.dwt.2d


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Multiple two dimensional wavelet transform
#'
#' @description performs two dimensional wavelet transforms on all images/maps of a data set.
#' @author Tran Vi-vi Elodie PERRIN
#'
#' @param x a matrix (for one map) or a three dimensional array (for several maps).
#'          Two first dimensions correspond to maps dimensions.
#'          The third one is the number of maps.
#' @param wf name of the wavelet filter to use in the decomposition.
#' @param J depth of the decomposition, must be a number less than or equal to log(min{M,N},2), with  M and N are respectively the number of rows and columns of each map.
#' @param boundary a character string which specified the method used for side effect. Only "periodic" is currently implemented.
#'
#' @importFrom waveslim dwt.2d
#'
#' @return an object of class "Wavelet2D", which corresponds to
#'          a matrix (for one map) or an array (for several maps)
#'          of same dimension as \code{x}, with the wavelet coefficients.
#'
#' @seealso \code{\link{Inverse2D.Wavelet2D}}  \code{\link{Inverse2D}}
#'
#' @examples
#' # inputs of the Campbell2D function
#' x1<-rep(-1,8); x2<-rep(5,8); x3<-c(5,3,1,-1,5,3,1,-1)
#' X <- rbind(x1,x2,x3)
#'
#' # spatial domain of the Campbell2D output
#' nz<-64
#' z<-seq(-90,90,length=nz)
#'
#' # Campbell2D function
#' Y = Campbell2D(X,z,z)
#'
#' # Wavelet transform
#' w <-Wavelet2D(Y,wf="d4",J=2)
#'
#' # Inverse wavelet transform
#' hatY <- Inverse2D(w)
#'
#' @export
Wavelet2D <-function(x, wf ,J=1, boundary ="periodic"){
  d <- dim(x) # maps dimension
  n<- d[3] # d[3] = number of maps

  if(length(d)==2){# if only one map
    w<-dwt.2d(x,wf=wf,J=J,boundary=boundary)# wavelet transform
    yi<-as.matrix(w,d[1],d[2])
    copy<-yi # asign wavelet coefficients
  }else{
    copy <- x # copy of x characteristics
    for(i in 1:n){
      w<-dwt.2d(x[,,i],wf=wf,J=J,boundary=boundary) # wavelet transform
      yi <-as.matrix(w,d[1],d[2])
      copy[,,i]<-yi # asign wavelet coefficients
    }# end for i
  }# end ifelse


  res <-copy
  attr(res,"type")<-attr(yi,"type")
  attr(res,"J")<-J  # depth of wavelet decomposition
  attr(res,"wf")<-wf #wavelet type
  attr(res,"boundary")<-boundary
  attr(res,"dim")<-d

  class(res)<-"Wavelet2D"

 return(res)
} # end Wavelet2D

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Multiple two-dimensional inverse wavelet transform
#'
#' @description performs inverse wavelet transform of each set of wavelet coefficients from a data set.
#' @author Tran Vi-vi Elodie PERRIN
#'
#' @param object an object of class \code{\link{Wavelet2D}}.
#' @param ... other arguments
#'
#' @return a three dimensional array.
#' The first two dimensions correspond to maps dimensions.
#' The third dimension corresponds to the size of the data set
#'
#' @importFrom waveslim idwt.2d
#' @method Inverse2D Wavelet2D
#'
#' @seealso \code{\link{Wavelet2D}}
#'
#' @export
Inverse2D.Wavelet2D <- function(object,...){
  # matrix with the associate wavelet tensorisation
  typ <- attr(object,"type")
  # wavelet type
  wf<-attr(object,"wf")
  # depth of the wavlet decomposition
  J<-attr(object,'J')
  # Boundary method
  boundary<-attr(object,'boundary')

  # dimensions
  d <- attr(object,"dim") # map dimension
  if(length(d)==2){n=1}else{n=d[3]}# number of images/maps


  # map approximations
  if(length(d)==2){
    res <-array(dim=c(d,n))
    }else{
    res <- array(dim=d)
    }# end ifelse
  for(i in 1:n){
    # wavelet coeffcients
    wi <- object[,,i]

    attr(wi,"type")<-attr(object,"type")
    attr(wi,"J")<-attr(object,"J") # depth of the wavelet decomposition
    attr(wi,"wf")<-attr(object,"wf") # wavelet type
    attr(wi,"boundary")<-attr(object,"boundary")
    attr(wi,"dim")<-d[1:2] # maps dimension

    class(wi)<-"matrix.dwt.2d"

    # inverse transform
    w <- as.dwt.2d(wi) # dwt.2d object
    res[,,i]<-idwt.2d(w)
  }# end i

  return(res)
}# end Inverse2D.Wavelet2D
