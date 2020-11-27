#' @title Gaussian Process Model on principal components of \code{Fpca2d},
#' by using \code{kergp} package
#'
#' @description  the function \code{gp} of the \code{\link{kergp}} package is use to
#' fit kriging models on each principal component modeled on \code{\link{Fpca2d}} (for more details, see \code{\link{gp}} ).
#'
#'
#' @param formula an object of class "formula" (or a list of "formula" which the length is equal to the number of modeled principal component)
#' specifying the linear trend of the kriging model (see \code{\link{lm}}) on each principal component.
#' This formula should concern only the input variables (\code{design}), and not the output (\code{response}).
#'  The default is ~1, which defines a constant trend on each principal component.
#' @param design a data frame representing the design of experiments.
#' The ith row contains the values of the d input variables corresponding
#' to the ith evaluation.
#' @param response n object of class \code{Fpca2d} which contains eigen
#' decomposition of the model/function ouput.
#' @param cov a covariance kernel object or call
#' (or a list of covariance kernel objects or call )
#'
#' @param estim Logical. If TRUE, the model parameters are estimated
#' by Maximum Likelihood. The initial values can then be specified using
#' the parCovIni and varNoiseIni arguments of mle,covAll-method passed though dots. If FALSE, a simple Generalized Least Squares estimation will be used, see gls,covAll-method. Then the value of varNoise must be given and passed through dots in case noise is TRUE.
#' @param ... other inputs of \code{\link{gp}}.
#'
#' @seealso \code{\link{gp}} \code{\link{kergp}}
#'
#' @importFrom kergp gp
#' @importFrom stats as.formula
#'
#' @return a list of object of class \code{\link{gp}} for each modeled principal component.
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
#' # change X on data.frame
#' colnames(X)<-paste("x",1:8,sep="")
#'
#' ############
#' ### FPCA ###
#' ############
#'
#' ### by using wavelet basis ###
#' fpca_w<- Fpca2d(Y,method="Wavelets",
#'                 wf="d4", J=1, # wavelet parameters
#'                 ncoeff=1200, rank.=2) # FPCA configuration
#'
#' #####################
#' ###               ###
#' ### Kriging model ###
#' ###               ###
#' #####################
#'
#' #------------------------------------#
#' #------------------------------------#
#' #   Example by using wavelet basis   #
#' #------------------------------------#
#' #------------------------------------#
#'
#' #--------------------------------------------#
#' #  Same kernel for all principal components  #
#' #--------------------------------------------#
#'
#' ## kernel
#' myCov <- covTS(inputs = colnames(X),
#' kernel = "k1Matern5_2",
#' dep = c(range = "input"),
#' value = c(range = 0.4))
#'
#' myGp<- gp_Fpca2d(design=X, response=fpca_w, cov=myCov,estim=FALSE)
#'
#' #-------------------------------------------------------------#
#' #  Different kernel and formula for each principal component  #
#' #-------------------------------------------------------------#
#'
#' \dontrun{
#' ## kernel of firt principal component
#' myCov1<-myCov
#'
#' ## kernel of second principal component
#' myCov2 <- covTS(inputs = colnames(X),
#' kernel = "k1Matern3_2",
#' dep = c(range = "input"),
#' value = c(range = 0.4))
#'
#' ## List of both kernels
#' myCovList <- list(myCov1,myCov2)
#'
#' ## Gp model
#' myGp2<- gp_Fpca2d(formula=list(~1,~x1+x2+x3+x4+x5+x6+x7+x8),
#'            design=X, response=fpca_w, cov=myCovList,estim=FALSE)
#' }
#' ##################
#' ### Prediction ###
#' ##################
#'
#' NewX<-matrix(runif(5*8,min=-1,max=5),ncol=8) # newdata
#' RealY <- Campbell2D(NewX,z,z)# real maps
#'
#' # change NewX on data.frame
#' colnames(NewX)<-colnames(X)
#'
#' #------------------------------------#
#' #------------------------------------#
#' #   Example by using wavelet basis   #
#' #------------------------------------#
#' #------------------------------------#
#' pw.UK <- predict(myGp,NewX,"UK")
#'
#' ###############################
#' ### Prediction RMSE and Q2  ###
#' ###############################
#'
#' #------------------------------------#
#' #------------------------------------#
#' #   Example by using wavelet basis   #
#' #------------------------------------#
#' #------------------------------------#
#' err.pw.UK <-error.predict(RealY,pw.UK,fpca_w,rtx.scores=TRUE)
#'
#' ### scores ###
#' print(err.pw.UK$scores$rmse)
#' print(err.pw.UK$scores$Q2)
#'
#' ### images/maps ###
#' library(fields)
#' image.plot(err.pw.UK$y$rmse, main="RMSE")
#' image.plot(err.pw.UK$y$Q2, main="Q2")
#'
#' @export
gp_Fpca2d <-function(formula=~1, design, response, cov,estim=TRUE,...){
  # scores of fpca
  y<-response$x; nPC <-ncol(y)

  #________________
  # formula to list
  #________________
  if(class(formula)=="formula"){
    res_i <-formula
    formulaPC <-(foreach(i=1:nPC,.combine = list,.multicombine = TRUE)%do%
                 res_i) # end foreach
  }else{
    formulaPC = formula
  }# end ifelse

  #---------------------------------------------
  # write formulas in the right form for each PC
  #---------------------------------------------
  formula <- lapply(1:nPC, function(i){
    output <- "yi" # output name of gp

    # as.character(formula)
    formula_character <- as.character(formulaPC[[i]])
    formule <- foreach(i=1:length(formula_character),.combine = paste,.multicombine = TRUE)%do%{
      formula_character[i]
    } # end formule

    past_formula <- paste(output,formule,sep="")

    return(as.formula(past_formula))
  })# end formula

  #________________

  #________________
  # function to list
  #________________
  if((!is.list(cov))){
    res_i <-cov
    cov <-(foreach(i=1:nPC,.combine = list,.multicombine = TRUE)%do%
                res_i)# end foreach
  }# end if
  #________________
  #%%%%%%%%%%%%%
  ### Models ###
  #%%%%%%%%%%%%%
  designNames<-colnames(design) # input names

  m <-lapply(1:nPC,function(i){
    yi<-as.numeric(y[,i])
    data.input <- data.frame(design,yi=yi)

    mi<-gp(formula=formula[[i]], data=data.input, inputs = designNames,cov=cov[[i]],
           estim=estim,...)
    return(mi)
  }) # end m

  class(m)<-"gp_Fpca2d"
  attr(m,"fpca")<-response

  return(m)

}# end gp_Fpca2d
