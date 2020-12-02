#' @title Prediction accuracy
#'
#' @description Prediction accuracy is measured by the spatial RMSE and Q2.
#' \code{error.predict} also computes RMSE and Q2 of the FPCA scores.
#'
#' @param y an array with the real two dimensional functional outputs.
#' The two first dimensions correspond to data dimensions,
#'  which are denoted M and N. The third one is the size of the data set.
#' @param py predictions which are given by \code{\link{predict.km_Fpca2d}} or \code{\link{predict.gp_Fpca2d}}.
#' @param fpca an object of class \code{\link{Fpca2d}}.
#' @param sd.epsilon a value with \code{y} standard deviation threshold. Q2 are NA, where
#' spatial standard deviations are less than or equal to \code{sd.epsilon}.
#' The default is 0. This is to avoid infinite values of Q2.
#' @param scores.epsilon a value with  \code{y}  score standard deviation threshold. Q2 of scores are NA, where
#' score standard deviations are less than or equal to \code{scores.epsilon}.
#' \code{scores.epsilon} are NA. The default is 0. This is to avoid infinite values of Q2.
#' @param rtx.scores a logical value. If TRUE, rmse and Q2 of the predicted scores are returned.
#' The default is FALSE
#'
#' @importFrom stats sd predict
#'
#' @return
#' \itemize{
#'   \item if rtx.scores=FALSE, a list with the following items is returned :
#'         \itemize{
#'         \item \code{rmse} : the spatial rmse of \code{y} prediction.
#'         \item \code{Q2} : the spatial Q2 of \code{y} prediction.
#'         }
#'   \item if rtx.scores=TRUE, a list with the following items is returned :
#'         \itemize{
#'         \item \code{scores} :
#'            \itemize{
#'                  \item \code{rmse} : rmse of \code{y} score prediction.
#'                   \item \code{Q2} : Q2 of \code{y} score prediction.
#'             }
#'             \item \code{y} :the same list as for the case rtx.scores=FALSE
#'         }
#' }
#'
#'
#' @export
error.predict <- function(y,py,fpca,sd.epsilon=0,scores.epsilon=0,rtx.scores=FALSE){

  ###################
  ### images/maps ###
  ###################
  res.2D <- list()# RMSE and Q2 of predicted images/maps

  # y variance
  sd.2D <-apply(y, c(1,2), sd) # standard variance
  var.2D<-sd.2D**2

  # maps MSE
  mse2D <- apply((y - py$mean)**2,c(1,2),mean)
  res.2D$rmse <- sqrt(mse2D)

  #maps Q2
  Q2 <- 1-mse2D/var.2D
  Q2[which(sd.2D<=sd.epsilon)]<-NA
  res.2D$Q2 <- Q2

  ###################
  ###   Scores   ###
  ###################

  if(isTRUE(rtx.scores)){
    n<-dim(y)[3] # number of maps

    # wavelet or B-splines coefficients
    if(attr(fpca,"method")=="wavelet"){
      # wavelet coefficients of y
      wy <- matrix(Wavelet2D(y,
                             wf=attr(fpca,"wf"),
                             J=attr(fpca,"J"),
                             boundary = attr(fpca,"boundary")),ncol=n)
    }

    # B-splines coefficients
    if(attr(fpca,"method")=="Bsplines"){
      # basis
      OPhi <-attr(fpca,"SplinesBasis")
      wy <- matrix(coef(OPhi,y),ncol=n)
    }

    #---------------------
    # coefficients on PCA
    #---------------------
    # poe
    poe <-attr(fpca,"mean_poe")
    ncoeff <-attr(fpca,"ncoeff")

    # coefficients on PCA
    idx_order <- order(poe, decreasing=TRUE)
    idx_pca <-idx_order[1:ncoeff]
    #--------------------

    # scores of y on fpca$pca
    pca <-attr(fpca,"pca")
    xy <- predict(pca,newdata=t(wy[idx_pca,]))

    # Score variance
    sd.scores <-apply(xy, 2, sd)
    var.scores<- sd.scores**2

    # predicted scores
    nPC <- ncol(xy) # number of principal components
    HatScores <- sapply(1:nPC, function(i){
      py$scores_predict[[i]]$mean
    }) # end HatScores

    # score errors
    res.scores <- list()

    # rmse.scores
    mse.scores <- apply((xy-HatScores)**2, 2, mean)
    res.scores$rmse <- sqrt(mse.scores)

    # Q2.scores
    Q2.scores <- 1-mse.scores/var.scores
    Q2.scores[which(sd.scores <= scores.epsilon)]=NA
    res.scores$Q2<-Q2.scores

    # function return
    res <-list()
    res$scores <- res.scores
    res$y<-res.2D
    class(res)<-"error.predict"
    return(res)

  }else{
    res <-res.2D
    class(res)<-"error.predict"
    return(res)
    }# end ifelse


}# end rmse
