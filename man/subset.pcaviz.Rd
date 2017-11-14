\name{subset.pcaviz}
\alias{subset.pcaviz}
\title{Subsetting of pcaviz objects}
\description{Return subset of observations in pcaviz object.}
\usage{
  \method{subset}{pcaviz}(x, subset, \dots)
}
\arguments{
  \item{x}{An object of class \code{"pcaviz"}.}

  \item{subset}{Logical expression indicating elements or rows to keep;
    missing values are taken as false.}

  \item{\dots}{Additional arguments to generic subset function. These
    additional arguments are not used here.}
}
\value{
  A pcaviz object in containing data and PCA results for the selected
  elements only.

  Note that unused levels are dropped from all factors (i.e.,
  categorical variables).
}
\seealso{
  \code{link{pcaviz}},
  \code{\link{subset}},
  \code{\link{droplevels}}
}
\examples{
  # Create the pcaviz object.
  data(iris)
  a <- pcaviz(prcomp(iris[1:4]),dat = iris)

  # Retain the subset of Iris samples in which the first PC is
  # positive. 
  b <- subset(a,PC1 > 0)
}
