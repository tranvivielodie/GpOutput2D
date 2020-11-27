#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Display a color image of the wavelet coefficients mean proportion of energy.
#'
#' @param x an object \code{\link{MeanPoe.Wavelet2D}}.
#' @param legend.col a logical value. Colors legend is given if \code{TRUE} (default value).
#' @param main an overall title for the plot: see \code{\link{title}}.
#' @param ... additional graphical parameters for \code{\link{plot}}.
#'
#' @author Tran Vi-vi Elodie PERRIN
#'
#' @importFrom fields image.plot
#' @import graphics
#' @method plot MeanPoe.Wavelet2D
#'
#'
#' @keywords internal
#' @export
plot.MeanPoe.Wavelet2D<-function(x,legend.col=TRUE,
                                 main="Mean proportion of energy",...){
  J <- attr(x,"J")# depth of the wavelet decomposition
  d <- dim(x) # maps dimension

  # plot "coordinates"
  z1<-1:d[1]; z2 <- 1:d[2]

  # image according if color legend is required
  if(legend.col==TRUE){
    image.plot(z1,z2,x,xaxt="n",yaxt="n",xlab="",ylab="",...)
  }else{
    image(z1,z2,x,xaxt="n",yaxt="n",xlab="",ylab="",...)
    }


  ######################################
  ### square to separate resolutions ###
  ######################################
  # plot limits
  al <-par("usr")
  names(al)<-c("xmin","xmax","ymin","ymax")

  # # first resolution
  # hv<-(al[c("xmax","ymax")]-al[c("xmin","ymin")])/2

  # first resolution squares
  lines(x=rep((al["xmin"]+d[1]/2),2),
        y=c(al["ymin"],al["ymax"])) # vertical
  lines(x=c(al["xmin"],al["xmax"]),
        y=rep((al["ymin"]+d[2]/2),2)) # horizontal

  # label location of wavelet information (tensor)
  at.lab <- al - (al/4)
  # at1 <- d[1]-(d[1]/4); at2 <-d[2]-(d[2]/4)

  ## H1
  axis(1,at=at.lab["xmax"],labels = "H1")
  axis(2,at=at.lab["ymax"],labels = "H1")

  ## L1
  if(J==1){
    at.lab <- al/4# at1 <- (d[1]/4); at2 <-(d[2]/4)
    axis(1,at=at.lab["xmax"],labels = "L1")
    axis(2,at=at.lab["ymax"],labels = "L1")
  }else{ # for J>1
    al0 <-al # graphic limits
    for(j in 1:(J-1)){
      alj<-al/2
      dj<-d[1:2]/(2**j)
      d12j<-dj/2
      # hvj <-(alj[c("xmax","ymax")]-alj[c("xmin","ymin")])/2# d12<-dj/2

      # j_th resolution squares
      xminj <- al0["xmin"]; yminj <- al0["ymin"]

      # vertical
      x.vj<-rep(xminj+d12j[1],2)
      y.vj<-c(yminj,yminj+dj[2])
      lines(x=x.vj,y=y.vj)

      # horizontal
      x.hj<-c(xminj,xminj+dj[1])
      y.hj<-rep(yminj+d12j[2],2)
      lines(x=x.hj,y=y.hj)


      # xj1 <- rep(d12[1],2); yj1 <- c(-1,dj[2]) # line 1
      # xj2 <- c(-1,dj[1]); yj2 <- rep(d12[2],2)# line 2
      #
      # lines(xj1+0.5,yj1)
      # lines(xj2,yj2+0.5)

      # label location of wavelet information (tensor)
      at.lab <- alj - (alj/4)
      # at1 <- dj[1]-(dj[1]/4); at2 <-dj[2]-(dj[2]/4)

      ## Hj
      hj <- paste("H",j+1,sep="")
      axis(1,at=at.lab["xmax"],labels = hj)
      axis(2,at=at.lab["ymax"],labels = hj)

      # axis(1,at=at1+0.5,labels = hj)
      # axis(2,at=at2+0.5,labels = hj)

      ## Lj
      if(J==j+1){
        at.lab <-alj/4
        # at1 <- (dj[1]/4); at2 <-(dj[2]/4)
        lj <- paste("L",j+1,sep="")
        axis(1,at=at.lab["xmax"],labels = lj)
        axis(2,at=at.lab["ymax"],labels = lj)
        # axis(1,at=at1+0.5,labels = lj )
        # axis(2,at=at2+0.5,labels = lj )
      }# end if


      # values update
      al <-alj
      # hv<-hvj
    }# end for j
  }# end if

  title(main)
}# end plot

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Display a color image of the wavelet coefficients mean proportion of energy.
#'
#' @param x an object \code{\link{MeanPoe.coef_OrthoNormalBsplines2D}}.
#' @param legend.col a logical value. Colors legend is given if \code{TRUE} (default value).
#' @param main an overall title for the plot: see \code{\link{title}}.
#' @param ... additional graphical parameters for \code{\link{plot}}.
#'
#' @author Tran Vi-vi Elodie PERRIN
#'
#' @importFrom fields image.plot
#' @import graphics
#' @method plot MeanPoe.coef_OrthoNormalBsplines2D
#'
#'
#' @keywords internal
#' @export
plot.MeanPoe.coef_OrthoNormalBsplines2D<-function(x,legend.col=TRUE,
                                 main="Mean proportion of energy",...){

  d<-dim(x)

  if(legend.col==TRUE){
    image.plot(x,xaxt="n",yaxt="n",...)
  }else{
    image(x,xaxt="n",yaxt="n",...)
  }
  title(main)
}# end plot
