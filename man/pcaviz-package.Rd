\name{PCAviz-package}
\alias{PCAviz-package}
\alias{pcaviz-package}
\alias{PCAviz}
\docType{package}
\title{Visualizing Principal Component Analysis}
\description{The \pkg{PCAviz} package provides a simple interface for
  quickly creating visually compelling plots from Principal Component
  Analysis (PCA) and accompanying data. See the vignettes included in
  this package for examples of plots generated with the \pkg{PCAviz}
  package. The plotting functions are also suitable for other
  dimensionality reduction techniques such as Multidimensional Scaling
  (MDS) and Local Linear Embedding (LLE).
}
\section{Introduction}{
  The \dQuote{iris}, \dQuote{popres} and \dQuote{regmap} vignettes give
  a quick overview of the package's key features and functions.

  PCAviz is specifically designed for visualizing three types of
  information jointly: (1) PCA results; (2) continuous covariates such
  as geographical co-ordinates; and (3) categorical data such as group
  labels. The \code{\link{pcaviz}} S3 class combines all these data
  stored into a single object. See \link{pcaviz-plots} for details on
  the PCAviz plotting interface.

  The main functions operate on \code{\link{pcaviz}} objects.  Most
  functions in the PCAviz package start with \code{pcaviz_}, so a list
  of available PCAviz functions can be quickly scanned using tab
  completion. For a complete list of help pages, use \code{help(package
  = "PCAviz")}.
}
\author{John Novembre, Richard Williams, Hossein Pourreza, Yuexi Wang and
    Peter Carbonetto}
\seealso{
  \url{http://github.com/NovembreLab/PCAviz}
}
\examples{
\dontrun{
# Demonstrations of the PCAviz package.
vignette("iris")
vignette("popres")
vignette("regmap")
}
}
