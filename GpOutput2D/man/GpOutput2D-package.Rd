\name{GpOutput2D-package}
\Rdversion{1.1}
\alias{GpOutput2D-package}
\alias{GpOutput2D}
\docType{package}

\title{Memodelling for model with two-dimensional functional output by using Functional Principal Component Analysis and kriging methods}

\description{Package for two-dimensional functional data analysis (for instance : data set of images or maps). Functions for performing Functional Principal Component Analysis (FPCA) on such data have been developped. Furthermore, the package contains methods for Computer Experiments by using kriging methods : metamodelling of models with two-dimensional functional output.}

\details{
  \tabular{ll}{
    Package : \tab GpOutput2D \cr
    Type : \tab Package \cr
    Version : \tab 0.4\cr
    Date : \tab 2020-11-15\cr
    License : \tab GPL (=2.0)\cr
    LazyLoad : \tab yes  \cr
  }
}

\author{Tran Vi-vi Élodie PERRIN (tran-vivi-elodie.perrin@outlook.fr)}

\note{ Important functions or methods :
  \itemize{
    \item \code{Fpca2d} : Functional Principal Component Analysis (FPCA) for two-dimensional data (images, maps,...).

   \item  \code{km_Fpca2d} : Gaussian Process Model on principal components of \code{Fpca2d}, by using \code{DiceKriging} package.

  \item   \code{gp_Fpca2d} : Gaussian Process Model on principal components of \code{Fpca2d}, by using \code{kergp} package.

  \item   \code{predict} :	Prediction of the objective function at new points using  \code{km_Fpca2d} or \code{gp_Fpca2d} model (Simple and Universal Kriging).
  }
}

\references{test}



