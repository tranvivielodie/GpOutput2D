% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OrthoNormalBsplines2D_utils.R
\name{coef.OrthoNormalBsplines2D}
\alias{coef.OrthoNormalBsplines2D}
\title{Get control points of two-dimensional B-splines basis}
\usage{
\method{coef}{OrthoNormalBsplines2D}(object, y, ...)
}
\arguments{
\item{object}{an object of class \code{\link{OrthoNormalBsplines2D}}
with functions of the B-splines basis.}

\item{y}{two-dimensional data to approximate on B-splines basis.}

\item{...}{other arguments.}
}
\value{
a matrix with control points on two-dimensional B-splines basis.
}
\description{
Get control points of two-dimensional B-splines basis
}
\seealso{
\code{\link{OrthoNormalBsplines2D}} \code{\link{Inverse2D.OrthoNormalBsplines2D}}
\code{\link{Inverse2D}}
}
