#' @title Predict values and confidence intervals at newdata for a km_Fpca2d object
#'
#' @param object an object of class \code{\link{km_Fpca2d}}.
#' @param newdata a vector, matrix or data frame containing the points where to perform predictions.
#' @param type a character string corresponding to the kriging family, to be chosen between simple kriging ("SK"), or universal kriging ("UK").
#' @param compute an optional boolean. If FALSE, only the kriging mean is computed. If TRUE, the kriging variance (here, the standard deviation is returned) and confidence intervals are computed too.
#' @param ... see \code{\link{predict.km}}.
#'
#' @importFrom stats predict
#' @importFrom graphics image
#' @method predict km_Fpca2d
#'
#' @return a list with the following items :
#' \itemize{
#'   \item \code{scores_predict} :  prediction of scores from two-dimensional FPCA
#'   (see  \code{\link{Fpca2d}} and \code{\link{predict.km}}).
#'   \item \code{mean} : an array with three dimensions which contains image (or map) kriging means.
#'   The two first dimensions corresponds to the output data dimensions. The third correspond to the number of predictions.
#'  \item \code{sd} : it is return only if compute = TRUE.
#'  If TRUE, an array which contains image (or map) prediction standard deviation.
#'  \item \code{lower95,upper95} :  bounds of the 95% confidence interval computed at \code{newdata}.
#'   Not computed if compute=FALSE.
#' }
#'
#' @seealso \code{\link{predict.km}}
#'
#'
#' @export
predict.km_Fpca2d <-function(object,newdata,type,compute = TRUE,...){

  # number of principal component
  nPC <- length(object)


  # number of prediction
  n<-nrow(newdata)

  # fpca object
  fpca <-attr(object,"fpca")

  # return
  res <-list()

  # prediction on each principal component
  p_PC <-lapply(1:nPC, function(i){
    pi <-predict(object=object[[i]],newdata=newdata,type=type,
                 compute = compute,...)
    return(pi)
  })# end p.PC

  res$scores_predict <-p_PC

  #%%%%%%%%%%%%%
  # prediction
  #%%%%%%%%%%%%%
  p_mean <-sapply(1:nPC, function(i){p_PC[[i]]$mean}) # scores
  hatY <- inverse_Fpca2d(p_mean,fpca)# maps
  res$mean <-hatY

  d <- dim(fpca$EigFct)# eigenfunction dimensions

  # scores sd
  p_sd <-sapply(1:nPC, function(i){p_PC[[i]]$sd}) # scores standard deviation
  if(isTRUE(compute)){ # if true return map prediction sd**2
    p_sd2 <- p_sd**2

    #%%%%%%%%%%%%%%%%%%%%%%%%%
    # map prediction variance
    #%%%%%%%%%%%%%%%%%%%%%%%%%
    sd2<-matrix(nrow=d[1]*d[2],ncol=n)
    xi <- matrix(fpca$EigFct,ncol=nPC)# eigenfunctions on vectors


    if(nPC==1){
      for(i in 1:n){
        XiZ<-p_sd2[i]*(xi%*%t(xi))
        sd2[,i] <- diag(XiZ)
      }# end for i
      sd2 <-array(sd2,dim = c(d[1:2],n))
    }else{
      for(i in 1:n){
        XiZ<-xi%*%diag(p_sd2[i,])%*%t(xi)
        sd2[,i] <- diag(XiZ)
      }# end for i
      sd2 <-array(sd2,dim = c(d[1:2],n))
    } # end ifelse
    res$sd <- sqrt(sd2) # standard deviation

    # 95% confidence interval
    res$lower95 <- hatY - 1.96*res$sd
    res$upper95 <- hatY + 1.96*res$sd
  }# end if

  return(res)

}# end predict.km_Fpca2d

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################

#' @title Predict values and confidence intervals at newdata for a km_Fpca2d object
#'
#' @param object an object of class \code{\link{gp_Fpca2d}}.
#' @param newdata a vector, matrix or data frame containing the points where to perform predictions.
#' @param type a character string corresponding to the kriging family, to be chosen between simple kriging ("SK"), or universal kriging ("UK").
#' @param compute an optional boolean. If FALSE, only the kriging mean is computed. If TRUE, the kriging variance (here, the standard deviation is returned) and confidence intervals are computed too.
#' @param ... see \code{\link{predict.km}}.
#'
#' @importFrom stats predict
#' @importFrom graphics image
#' @method predict gp_Fpca2d
#'
#' @return a list with the following items :
#' \itemize{
#'   \item \code{scores_predict} :  prediction of scores from two-dimensional FPCA
#'   (see  \code{\link{Fpca2d}} and \code{\link{predict.km}}).
#'   \item \code{mean} : an array with three dimensions which contains image (or map) kriging means.
#'   The two first dimensions corresponds to the output data dimensions. The third correspond to the number of predictions.
#'  \item \code{sd} : it is return only if compute = TRUE.
#'  If TRUE, an array which contains image (or map) prediction standard deviation.
#'  \item \code{lower95,upper95} :  bounds of the 95% confidence interval computed at \code{newdata}.
#'   Not computed if compute=FALSE.
#' }
#'
#' @seealso \code{\link{predict.km}}
#'
#' @export
predict.gp_Fpca2d <-function(object,newdata,type,compute = TRUE,...){

  # number of principal component
  nPC <- length(object)

  # number of prediction
  n<-nrow(newdata)

  # fpca object
  fpca <-attr(object,"fpca")

  # return
  res <-list()

  # prediction on each principal component
  p_PC <-lapply(1:nPC, function(i){
    pi <-predict(object=object[[i]],newdata=newdata,type=type,
                 compute = compute,...)
    return(pi)
  })# end p.PC

  res$scores_predict <-p_PC

  #%%%%%%%%%%%%%
  # prediction
  #%%%%%%%%%%%%%
  p_mean <-sapply(1:nPC, function(i){p_PC[[i]]$mean}) # scores
  hatY <- inverse_Fpca2d(p_mean,fpca)# maps
  res$mean <-hatY

  d <- dim(fpca$EigFct)# eigenfunction dimensions

  # scores sd
  p_sd <-sapply(1:nPC, function(i){p_PC[[i]]$sd}) # scores standard deviation
  if(isTRUE(compute)){ # if true return map prediction sd**2
    p_sd2 <- p_sd**2

    #%%%%%%%%%%%%%%%%%%%%%%%%%
    # map prediction variance
    #%%%%%%%%%%%%%%%%%%%%%%%%%
    sd2<-matrix(nrow=d[1]*d[2],ncol=n)
    xi <- matrix(fpca$EigFct,ncol=nPC)# eigenfunctions on vectors

    if(nPC==1){
      for(i in 1:n){
        XiZ<-p_sd2[i]*(xi%*%t(xi))
        sd2[,i] <- diag(XiZ)
      }# end for i
      sd2 <-array(sd2,dim = c(d[1:2],n))
    }else{
      for(i in 1:n){
        XiZ<-xi%*%diag(p_sd2[i,])%*%t(xi)
        sd2[,i] <- diag(XiZ)
      }# end for i
      sd2 <-array(sd2,dim = c(d[1:2],n))
    } # end ifelse


    res$sd <- sqrt(sd2) # standard deviation

    # 95% confidence interval
    res$lower95 <- hatY - 1.96*res$sd
    res$upper95 <- hatY + 1.96*res$sd
  }# end if

  return(res)

}# end predict.km_Fpca2d
