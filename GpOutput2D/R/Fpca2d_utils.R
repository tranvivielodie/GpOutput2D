#' @title Plot for Fpca2d object
#'
#' @description plots characteristic of Fpca2d decomposition : eigenfunctions,
#'              coefficients total mean proportion of energy, coefficients truncation etc.
#'
#' @param x an object of class "Fpca2d".
#' @param type a character vector, which the characteristics of \code{x}
#'        are plotted. The default is c("inertia","energy","MeanPoe","SelectedCoeff").
#'        It means that the \code{barplot} of proportion of inertia of each principal
#'        component (\code{type}="inertia"), the one of the number of coefficients
#'        estimated with PCA according to the total mean proportion of energy (\code{type}="energy"),
#'        the image of the coefficients mean proportion of energy (\code{type}="MeanPoe"),
#'        and the image which indicates which coefficients are estimated
#'        by PCA or by empirical mean (\code{type}="SelectedCoeff"). These graphics can be separetely plotted.
#'        The eigenfunctions can be plotted be specifying (\code{type}="eigenfunctions").
#'
#' @param PC a vector with the numbers of which principal components are plotted. The default is c(1,2).
#' @param p a vector with the proportion of energy which are plotted. The default is seq(0.5,1,by=0.05).
#' @param z1,z2 spatial coordinates.
#' @param ... plot parameters (see \code{\link{plot}}, \code{\link{image}} or \code{\link{image.plot}}).
#'
#' @import graphics
#' @importFrom fields image.plot
#' @method plot Fpca2d
#'
#' @export
plot.Fpca2d<-function(x,type=c("inertia","energy","MeanPoe"),
                      PC=1:2,
                      p=seq(0.5,1,by=0.05),
                      z1=NULL,z2=NULL,...){


  #########################################################
  ### proportion of inertia of each principal component ###
  #########################################################
  if("inertia"%in%type){
    pi <- round(x$sdev[PC]**2/sum(x$sdev**2),3) # proportion of inertia
    names(pi)<-PC

    # barplot
    bp<- barplot(pi, ylim=c(0,1),main="Proportion of inertia",col="skyblue",
              xlab="principal components", ylab="proportion of inertia")
    abline(h=0)
    text(bp, pi+0.05, labels = paste(pi*100,"%",sep=""),cex=1,col="blue")
  }# end if

  #########################################################
  ### proportion of energy / number of coefficients ###
  #########################################################
  if("energy"%in%type){
    poe <- as.vector(attr(x,"mean_poe")) # proportion of energy
    idx <- order(poe,decreasing=TRUE)# coefficients order
    CumPoe <- cumsum(poe[idx])

    # number of coefficients
    ncoeff <- sapply(p, function(pi){
       # if p==1 (<=> total energy)
        if(pi==1){
          # choose the minimum number of coefficients such as CumPoe=1
          ncoeffi <- min(which(CumPoe==1))
        }else{
          # choose the number of coefficients such as CumPoe<=p
          ncoeffi <- length(which(CumPoe<=pi))
        }# end ifelse (p==1)
      return(ncoeffi)
    })# end ncoeff

    names(ncoeff)<-p

    # barplot
    bp<-barplot(ncoeff,main="Number of coefficients", col="orange",xaxt="n",
            ylim=c(0,(K<-length(poe))), # K is the total number of coefficient
            xlab="Total mean proportion of energy", ylab="number of coefficients")

    abline(h=0)
    axis(1, at=bp, labels = p,cex=1)

    # percentage of kept coefficients
    p.coeff <- round(ncoeff/K,3)*100

    text(bp, ncoeff/2, labels = paste(p.coeff,"%",sep=""),cex=1)
  }# end if

  ##############################################
  ### configuration error (due to S3 method) ###
  ##############################################
  if(("eigenfunctions"%in%type) & ("MeanPoe"%in%type)){
    stop("eigenfunctions and MeanPoe cannot be plot with one run.
         Please use separate plot.")
  } # end if
  if(("eigenfunctions"%in%type) & ("SelectedCoeff"%in%type)){
    stop("eigenfunctions, and SelectedCoeff cannot be plot with one run.
         Please use separate plot.")
  }# end if
  if((!("MeanPoe"%in%type)) & ("SelectedCoeff"%in%type)){
    stop("The image of SelectedCoeff cannot be plot without MeanPoe.
         On the other hand, MeanPoe can be plotted without SelectedCoeff.")
  }# end if

  ######################
  ### eigenfunctions ###
  ######################
  if("eigenfunctions"%in%type){
    pi <- round(x$sdev[PC]**2/sum(x$sdev**2),3) # proportion of inertia

      if(is.null(z1)|(is.null(z1))){
        for(i in PC){
          image.plot(x$EigFct[,,i],...)
          title(paste("Principal component number",i,"\n",
                      paste(pi[i]*100,"%",sep="")))
        } # end for i
      }else{
        for(i in PC){
          image.plot(z1,z2,x$EigFct[,,i],...)
          title(paste("Principal component number",i,"\n",
                      paste(pi[i]*100,"%",sep="")))
        } # end for i
      }# end ifelse z1 and z2
  }# end if eigenfunctions

  ###############
  ### MeanPoe ###
  ###############
  if("MeanPoe"%in%type){
    poe <- attr(x,"mean_poe")
      plot(poe,...)
  }# end if MeanPoe

} # end plot.Fpca2d

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################


#' @title Get approximation with scores of Fpca2d
#'
#' @param scores a matrix where columns correpond to scores of principal components.
#' @param fpca an object of class \code{Fpca2d}.
#'
#' @return approximation of scores.
#'
#'
#' @keywords internal
inverse_Fpca2d <-function(scores,fpca){

  # number of principal components
  nPC <- ifelse(is.matrix(scores), # condition
                ncol(scores), # TRUE
                length(scores)) # FALSE

  # number of scores
  n <- ifelse(is.matrix(scores), # condition
              nrow(scores), # TRUE
              1) # FALSE

  # coefficients of wavelet or B-splines basis
  coeff <- attr(fpca,"coeff")
  d<-dim(coeff)[1:2]
  K <- d[1]*d[2] # size of the functional basis

  # number of coefficients on PCA
  ncoeff <-attr(fpca,"ncoeff")

  # poe
  poe <-attr(fpca,"mean_poe")

  # coefficients on PCA
  idx_order <- order(poe, decreasing=TRUE)
  idx_pca <-idx_order[1:ncoeff]


  mu <-as.vector(apply(coeff, c(1,2), mean)) # mean coefficients

  #"---------------------------------"
  # estimation of coefficients on PCA
  #"---------------------------------"
  # rotation matrix
  rot <- attr(fpca,"pca")$rotation
  coeff_pca <- rot%*%t(scores)

  # estimation of all coefficients
  if(isTRUE(fpca$center)){
    # add mean values
    res_coeff <- sapply(1:n, function(i){
      ci <- rep(0,K)
      ci[idx_pca] <- coeff_pca[,i]
      return(ci+mu)
    })# end res_coeff
  }else{
    res_coeff <- sapply(1:n, function(i){
      ci <- rep(0,K)
      # coefficients on PCA
      ci[idx_pca] <- coeff_pca[,i]
      # NOT on PCA
      ci[-idx_pca] <- mu[-idx_pca]
      return(ci)
    })# end res_coeff
  }# end ifelse center

  # coefficients
  if(attr(fpca,"method")=="Bsplines"){
    res_coeff_mat <- res_coeff # keep matrix for B-splines
    res_coeff <- array(res_coeff_mat,dim=c(d[1],d[2],n))
    attr(res_coeff,"matrix")<-res_coeff_mat
    rm(res_coeff_mat) # memory deallocation
  }else{
    res_coeff <- array(res_coeff,dim=c(d[1],d[2],n))
  }# end ifelse



  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #  if wavelet basis
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(attr(fpca,"method")=="wavelet"){

    attr(res_coeff,"type")<-attr(fpca,"type")
    attr(res_coeff,"J")<-attr(fpca,"J")  # depth of wavelet decomposition
    attr(res_coeff,"wf")<-attr(fpca,"wf") #wavelet type
    attr(res_coeff,"boundary")<-attr(fpca,"boundary")
    attr(res_coeff,"dim")<-c(d,n)

    class(res_coeff)<-"Wavelet2D"

    # approximation
    return(Inverse2D(res_coeff))

  }# end if wavelet

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #  if B-splines basis
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(attr(fpca,"method")=="Bsplines"){
    OPhi<-attr(fpca,"SplinesBasis")
    return(Inverse2D(OPhi,res_coeff))
  }# end if Bsplines

} # end inverse_Fpca2d
