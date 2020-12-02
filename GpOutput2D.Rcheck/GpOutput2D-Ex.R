pkgname <- "GpOutput2D"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "GpOutput2D-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('GpOutput2D')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Campbell2D")
### * Campbell2D

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Campbell2D
### Title: Example of function called Campbell2D
### Aliases: Campbell2D

### ** Examples

# inputs of the Campbell2D function
x1<-rep(-1,8); x2<-rep(5,8); x3<-c(5,3,1,-1,5,3,1,-1)
X <- rbind(x1,x2,x3)

# spatial domain of the Campbell2D output
nz<-64 # root of the size of the spatial domain
z<-seq(-90,90,length=nz)

# Campbell2D function
Y = Campbell2D(X,z,z)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Campbell2D", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Fpca2d")
### * Fpca2d

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Fpca2d
### Title: Two-dimensional Functional Principal Component Analysis (FPCA).
### Aliases: Fpca2d

### ** Examples


################################
### two-dimensional data set ###
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

############
### FPCA ###
############


### using wavelet basis ###
fpca_w<- Fpca2d(Y,method="Wavelets",
                wf="d4", J=1, # wavelet parameters
                ncoeff=1200, rank.=5) # FPCA configuration

plot(fpca_w,type="eigenfunctions")
plot(fpca_w)

## Not run: 
##D ### using B-splines basis ###
##D 
##D # knots for B-splines basis
##D K<-35
##D z.knots <- seq(-90,90,length=K)
##D 
##D fpca_Bs<- Fpca2d(Y,method="Bsplines",
##D                 z1=z,z2=z,z1.knots=z.knots,z2.knots=z.knots, # wavelet parameters
##D                 ncoeff=1225, rank.=2) # FPCA configuration
##D 
##D plot(fpca_Bs,type="eigenfunctions")
##D plot(fpca_Bs)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Fpca2d", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("OrthoNormalBsplines2D")
### * OrthoNormalBsplines2D

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: OrthoNormalBsplines2D
### Title: Creating a two-dimensional othonormal B-splines basis
### Aliases: OrthoNormalBsplines2D

### ** Examples


#####################################
#
#  Create two-dimensional data set
#
####################################

# inputs of the Campbell2D function
x1<-rep(-1,8); x2<-rep(5,8); x3<-c(5,3,1,-1,5,3,1,-1)
X <- rbind(x1,x2,x3)

# spatial domain of the Campbell2D output
nz<-64 # root of the size of the spatial domain
z<-seq(-90,90,length=nz)

# Campbell2D function
Y = Campbell2D(X,z,z)

###########################################
#
#  Create two-dimensional B-splines basis
#
###########################################

# knots for B-splines basis
K<-35
z.knots <- seq(-90,90,length=K)

# Generating a two-dimensional othonormal B-splines basis
OPhi <-OrthoNormalBsplines2D(z,z,z.knots,z.knots,ortho="GS")

########################
#
#  Get control points
#
########################

coeff_Y<-coef(OPhi,Y)

#######################
#
#  Approximation of Y
#
########################

hatY <-Inverse2D(OPhi,coeff_Y)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("OrthoNormalBsplines2D", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Wavelet2D")
### * Wavelet2D

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Wavelet2D
### Title: Multiple two dimensional wavelet transform
### Aliases: Wavelet2D

### ** Examples

# inputs of the Campbell2D function
x1<-rep(-1,8); x2<-rep(5,8); x3<-c(5,3,1,-1,5,3,1,-1)
X <- rbind(x1,x2,x3)

# spatial domain of the Campbell2D output
nz<-64
z<-seq(-90,90,length=nz)

# Campbell2D function
Y = Campbell2D(X,z,z)

# Wavelet transform
w <-Wavelet2D(Y,wf="d4",J=2)

# Inverse wavelet transform
hatY <- Inverse2D(w)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Wavelet2D", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gp_Fpca2d")
### * gp_Fpca2d

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gp_Fpca2d
### Title: Gaussian Process Model on principal components of 'Fpca2d', by
###   using 'kergp' package
### Aliases: gp_Fpca2d

### ** Examples


################################
### two-dimensional data set ###
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

############
### FPCA ###
############

### by using wavelet basis ###
fpca_w<- Fpca2d(Y,method="Wavelets",
                wf="d4", J=1, # wavelet parameters
                ncoeff=1200, rank.=2) # FPCA configuration

#####################
###               ###
### Kriging model ###
###               ###
#####################

#------------------------------------#
#------------------------------------#
#   Example by using wavelet basis   #
#------------------------------------#
#------------------------------------#

#--------------------------------------------#
#  Same kernel for all principal components  #
#--------------------------------------------#

## kernel
myCov <- covTS(inputs = colnames(X),
kernel = "k1Matern5_2",
dep = c(range = "input"),
value = c(range = 0.4))

myGp<- gp_Fpca2d(design=X, response=fpca_w, cov=myCov,estim=FALSE)

#-------------------------------------------------------------#
#  Different kernel and formula for each principal component  #
#-------------------------------------------------------------#

## Not run: 
##D ## kernel of firt principal component
##D myCov1<-myCov
##D 
##D ## kernel of second principal component
##D myCov2 <- covTS(inputs = colnames(X),
##D kernel = "k1Matern3_2",
##D dep = c(range = "input"),
##D value = c(range = 0.4))
##D 
##D ## List of both kernels
##D myCovList <- list(myCov1,myCov2)
##D 
##D ## Gp model
##D myGp2<- gp_Fpca2d(formula=list(~1,~x1+x2+x3+x4+x5+x6+x7+x8),
##D            design=X, response=fpca_w, cov=myCovList,estim=FALSE)
## End(Not run)
##################
### Prediction ###
##################

NewX<-matrix(runif(5*8,min=-1,max=5),ncol=8) # newdata
RealY <- Campbell2D(NewX,z,z)# real maps

# change NewX on data.frame
colnames(NewX)<-colnames(X)

#------------------------------------#
#------------------------------------#
#   Example by using wavelet basis   #
#------------------------------------#
#------------------------------------#
pw.UK <- predict(myGp,NewX,"UK")

###############################
### Prediction RMSE and Q2  ###
###############################

#------------------------------------#
#------------------------------------#
#   Example by using wavelet basis   #
#------------------------------------#
#------------------------------------#
err.pw.UK <-error.predict(RealY,pw.UK,fpca_w,rtx.scores=TRUE)

### scores ###
print(err.pw.UK$scores$rmse)
print(err.pw.UK$scores$Q2)

### images/maps ###
library(fields)
image.plot(err.pw.UK$y$rmse, main="RMSE")
image.plot(err.pw.UK$y$Q2, main="Q2")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gp_Fpca2d", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("km_Fpca2d")
### * km_Fpca2d

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: km_Fpca2d
### Title: Gaussian Process Model on principal components of 'Fpca2d', by
###   using 'DiceKriging' package.
### Aliases: km_Fpca2d

### ** Examples


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
## Not run: 
##D mw <- km_Fpca2d(formula=list(~.,~1)),
##D                 design=X,response=fpca_w,
##D                 covtype = c("matern5_2","matern3_2"),
##D                 control=list(trace=FALSE))
## End(Not run)

#--------------------------------------------------------
# (same example) how to use the multistart argument of km
#--------------------------------------------------------
## Not run: 
##D nCores <- 2
##D require(doParallel)
##D cl <-  makeCluster(nCores)
##D registerDoParallel(cl)
##D mw <- km_Fpca2d(design=X,response=fpca_w,multistart=4,control=list(trace=FALSE))
##D stopCluster(cl)
## End(Not run)

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

## Not run: 
##D ############
##D ### FPCA ###
##D ############
##D 
##D ### using B-splines basis ###
##D 
##D # knots for B-splines basis
##D K<-35
##D z.knots <- seq(-90,90,length=K)
##D 
##D fpca_Bs<- Fpca2d(Y,method="Bsplines",
##D                 z1=z,z2=z,z1.knots=z.knots,z2.knots=z.knots, ortho="GS",
##D                 expand_knots=TRUE,# B-splines parameters
##D                 ncoeff=1225, rank.=5) # FPCA configuration
##D 
##D 
##D #####################
##D ### Kriging model ###
##D #####################
##D 
##D mB <- km_Fpca2d(design=X,response=fpca_Bs,control=list(trace=FALSE))
##D 
##D ##################
##D ### Prediction ###
##D ##################
##D 
##D pB.UK <- predict(mB,NewX,"UK")
##D 
##D ########################
##D ###   RMSE and Q2    ###
##D ########################
##D 
##D err.pB.UK <-error.predict(RealY,pB.UK,fpca_Bs,rtx.scores=TRUE)
##D 
##D ### scores ###
##D print(err.pB.UK$scores$rmse)
##D print(err.pB.UK$scores$Q2)
##D 
##D ### images/maps ###
##D library(fields)
##D image.plot(err.pB.UK$y$rmse, main="RMSE")
##D image.plot(err.pB.UK$y$Q2, main="Q2")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("km_Fpca2d", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
