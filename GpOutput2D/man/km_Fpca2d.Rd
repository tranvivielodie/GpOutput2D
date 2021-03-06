% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/km_Fpca2d.R
\name{km_Fpca2d}
\alias{km_Fpca2d}
\title{Gaussian Process Model on principal components of \code{Fpca2d},
by using \code{DiceKriging} package.}
\usage{
km_Fpca2d(
  formula = ~1,
  design,
  response,
  parallel = FALSE,
  covtype = "matern5_2",
  coef.trend = NULL,
  coef.cov = NULL,
  coef.var = NULL,
  nugget = NULL,
  noise.var = NULL,
  lower = NULL,
  upper = NULL,
  parinit = NULL,
  multistart = 1,
  kernel = NULL,
  ...
)
}
\arguments{
\item{formula}{an object of class "formula"
(or a list of "formula" which the length is equal to the number of modeled principal components)
specifying the linear trend of the kriging model (see \code{\link{lm}}) on each principal component.
 This formula should concern only the input variables (\code{design}), and not the output (\code{response}).
 The default is ~1, which defines a constant trend on each principal component.}

\item{design}{a data frame representing the design of experiments.
The ith row contains the values of the d input variables corresponding
to the ith evaluation.}

\item{response}{an object of class \code{Fpca2d} which contains eigen
decomposition of the model/function ouput.}

\item{parallel}{a logical value specifying if parallelization is done on principal components.
The default is FALSE.}

\item{covtype}{optional character string or vector of character strings
specifying the covariance structure to be used on each modeled principal component
(see \code{\link{km}} for possible inputs of \code{covtype}).
If a vector, the length should be equal to the number of modeled principal components.}

\item{coef.trend, coef.cov, coef.var}{optional vectors or matrices containing
the values for the trend, covariance and variance parameters.
If matrices, the number of rows should be equal to the number of modeled principal components.
For details, see \code{\link{km}}).}

\item{nugget}{an optional variance value or vector standing for the homogeneous nugget effect.
If vector, the length should be equal to the number of modeled principal components.}

\item{noise.var}{an optional vector or matrix containing the noise variance
at each observation on each modeled principal component.}

\item{lower, upper}{optional vectors or matrices containing the bounds of the correlation parameters
of each principal component for optimization. For details, see \code{\link{km}}).
If matrices, the number of rows should be equal to the number of modeled principal components.}

\item{parinit}{an optional vector or matrix containing the initial values for the variables to be optimized over.
For details, see \code{\link{km}}).}

\item{multistart}{an optional integer indicating the number of initial points from which running the BFGS optimizer.
(see \code{\link{km}}).}

\item{kernel}{an optional function or list of functions containing a new covariance structure
for each principal component. At this stage, the parameters must be provided as well, and are not estimated.}

\item{...}{other parameters of \code{\link{km}} function from \code{DiceKriging}.
(see \code{\link{DiceKriging}})}
}
\value{
a list of object of class \code{\linkS4class{km}} for each modeled principal component.
}
\description{
the function \code{km} of the \code{\link{DiceKriging}} package is use to
fit kriging models on each principal component modeled on \code{\link{Fpca2d}} (for more details, see \code{\link{km}} ).
}
\examples{

################################
###   Learning Sample        ###
################################
n<-200 # size of the learning sample
nz<-64; z <- seq(-90,90,length=nz) # spatial domain

### inputs of Campbell2D ###
library(lhs)
library(DiceDesign)

x <- maximinLHS(n=n,k=8)
X <-maximinSA_LHS(x)$design
X<-X*6 -1

# Campbell2D
Y <- Campbell2D(X,z,z)

# change X on data.frame
colnames(X)<-paste("x",1:8,sep="")

############################
###   Test Sample        ###
############################

NewX<-matrix(runif(5*8,min=-1,max=5),ncol=8) # newdata
RealY <- Campbell2D(NewX,z,z)# real maps

# change NewX on data.frame
colnames(NewX)<-colnames(X)


######################################
#
#   Example by using wavelets
#
######################################

############
### FPCA ###
############

fpca_w<- Fpca2d(Y,method="Wavelets",
                wf="d4", J=1, # wavelet parameters
                ncoeff=1200, rank.=2) # FPCA configuration

#####################
### Kriging model ###
#####################

mw <- km_Fpca2d(design=X,response=fpca_w,control=list(trace=FALSE))

#--------------------------------------------------------
# (same example) To fix different kernel and formula
#                for each principal component
#--------------------------------------------------------
\dontrun{
mw <- km_Fpca2d(formula=list(~.,~1)),
                design=X,response=fpca_w,
                covtype = c("matern5_2","matern3_2"),
                control=list(trace=FALSE))
}

#--------------------------------------------------------
# (same example) how to use the multistart argument of km
#--------------------------------------------------------
\dontrun{
nCores <- 2
require(doParallel)
cl <-  makeCluster(nCores)
registerDoParallel(cl)
mw <- km_Fpca2d(design=X,response=fpca_w,multistart=4,control=list(trace=FALSE))
stopCluster(cl)
}

######################
###   Prediction   ###
######################

pw.UK <- predict(mw,NewX,"UK")

####################
### RMSE and Q2  ###
####################

err.pw.UK <-error.predict(RealY,pw.UK,fpca_w,rtx.scores=TRUE)

### scores ###
print(err.pw.UK$scores$rmse)
print(err.pw.UK$scores$Q2)

### images/maps ###
library(fields)
image.plot(err.pw.UK$y$rmse, main="RMSE")
image.plot(err.pw.UK$y$Q2, main="Q2")


######################################
#
#   Example by using B-splines
#
######################################

\dontrun{
############
### FPCA ###
############

### using B-splines basis ###

# knots for B-splines basis
K<-35
z.knots <- seq(-90,90,length=K)

fpca_Bs<- Fpca2d(Y,method="Bsplines",
                z1=z,z2=z,z1.knots=z.knots,z2.knots=z.knots, ortho="GS",
                expand_knots=TRUE,# B-splines parameters
                ncoeff=1225, rank.=5) # FPCA configuration


#####################
### Kriging model ###
#####################

mB <- km_Fpca2d(design=X,response=fpca_Bs,control=list(trace=FALSE))

##################
### Prediction ###
##################

pB.UK <- predict(mB,NewX,"UK")

########################
###   RMSE and Q2    ###
########################

err.pB.UK <-error.predict(RealY,pB.UK,fpca_Bs,rtx.scores=TRUE)

### scores ###
print(err.pB.UK$scores$rmse)
print(err.pB.UK$scores$Q2)

### images/maps ###
library(fields)
image.plot(err.pB.UK$y$rmse, main="RMSE")
image.plot(err.pB.UK$y$Q2, main="Q2")
}
}
\seealso{
\code{\link{km}}  \code{\link{DiceKriging}}
}
\author{
Tran Vi-vi Elodie PERRIN
}
