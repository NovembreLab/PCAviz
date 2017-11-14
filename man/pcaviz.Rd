\name{pcaviz}
\alias{pcaviz}
\alias{summary.pcaviz}
\alias{print.pcaviz}
\alias{print.summary.pcaviz}
\title{Create and summarize pcaviz objects}
\description{Creates a "pcaviz" data structure from PCA results (e.g.,
  \code{\link[stats]{prcomp}} output, and accompanying data. The data
  structure is specifically designed to facilitate plotting PCA results
  and relationships between estimated principal components and related
  data (e.g., geographic co-ordinates).}
\usage{
pcaviz(out.pca, x = NULL, sdev = NULL, var = NULL, rotation = NULL,
       dat = NULL, pc.cols)

\method{summary}{pcaviz}(object, \dots)

\method{print}{summary.pcaviz}(x, n = 4, \dots)

\method{print}{pcaviz}(x, n = 4, \dots)
}
\arguments{
  \item{out.pca}{Output of \code{prcomp} or \code{princomp} (\pkg{stats}
    package), or output of \code{rpca} (\pkg{rsvd} package), providing
    results of principal components analysis. Alternatively, PCA results
    (or results of applying another dimensionality reduction technique)
    may be provided via inputs \code{x}, \code{rotation} and \code{dat}
    (see below).}
  
  \item{x}{For function "pcaviz", this is a matrix or data frame
    containing the data samples projected onto the principal components
    (also sometimes referred to as "rotated data" or "scores"). This is
    equivalent to the \code{x} component returned by functions
    \code{\link[stats]{prcomp}} and \code{\link[rsvd]{rpca}}, and the
    \code{scores} component returned by function
    \code{\link[stats]{princomp}}. For function \code{print.pcaviz},
    this should be an object of class "pcaviz"; for function
    \code{print.summary.pcaviz}, this should be the output of function
    \code{summary.pcaviz}.}

  \item{sdev}{Optional vector containing the standard deviations of the
    PCs or, equivalently, the square root of the eigenvalues. This is
    the same as the \code{sdev} component returned by functions
    \code{\link[stats]{prcomp}}, \code{\link[stats]{princomp}} and
    \code{\link[rsvd]{rpca}}.}

  \item{var}{Optional value specifying the total variance. This should
    be equal to the sum of the variances of the individual PCs. It can
    be computed as \code{sum(apply(X,2,var))}, in which \code{X} is the
    target matrix supplied as input to the PCA. This is the same as
    \code{sum(sdev^2)} when \code{sdev} is returned by
    \code{\link[stats]{prcomp}} or \code{\link[stats]{princomp}}, and is
    the same as the \code{var} component returned by
    \code{\link[rsvd]{rpca}}.}
  
  \item{rotation}{Optional matrix or data frame containing the
    eigenvectors (alternatively called "variable loadings" or
    "rotation matrix"). This is can be chosen as the \code{rotation}
    component returned by functions \code{\link[stats]{prcomp}} and
    \code{\link[rsvd]{rpca}}, or the \code{loadings} component returned
    by function \code{\link[stats]{princomp}}. Additional data columns
    beyond the eigenvectors may also be supplied so long as all PC
    loadings are included.}

  \item{dat}{Optional data frame supplying additional information about
    the data samples. It should be a data frame with one row for each
    data sample, and should have as many rows as the rotated data.
    Optionally, this data frame may also supply the rotated data, as
    indicated by input \code{pc.cols}. This is useful if the PCA results
    and accompanying data are already combined into one table.}
  
  \item{pc.cols}{Indicates which columns of input \code{dat} contain the
    principal components (rotated data). The names of the columns (not
    numeric indices) must be supplied. By default, when \code{pc.cols}
    is not specified and the rotated data are not already supplied by
    \code{x} or \code{out.pca}, the columns of \code{dat} with prefix
    "PC" are interpreted as principal components.}

  \item{object}{An object of class \code{"pcaviz"}.}
  
  \item{n}{Maximum number of variables of each category to display in
    the summary.}
  
  \item{\dots}{Additional print and summary arguments. These additional
    arguments are not used here.}
}
\value{
  
  \code{pcaviz} returns a list with class \code{"pcaviz"} containing the
  following components:

  \item{data}{Data frame containing the rotated variables as well as
    accompanying data.}

  \item{data.coltypes}{Names of the columns in \code{data} that contain
    the rotated variables.}

  \item{sdev}{If not \code{NULL}, this is a vector containing the
    standard deviations of the PCs; equivalently, the square roots of
    the eigenvalues.}
  
  \item{rotation}{If not \code{NULL}, this is a data frame containing
    the eigenvectors ("loadings") and, optionally, other data columns.}
  
  \item{basis}{A matrix with one column per PC that is used to keep
    track of the principal axes after transformation.}

  \item{transformed.pcs}{A logical value indicating whether at least one
    transformation (e.g., a rotation) has been applied to the PCs.
    Initially set to FALSE.}

  It is advised to avoid directly manipulating these components as much
  as possible. Instead, it is preferrable to use package functions
  designed for manipulating \code{pcaviz} objects, e.g.,
  \code{\link[=subset.pcaviz]{subset}}, \code{\link{pcaviz_rotate}},
  \code{\link{pcaviz_abbreviate_var}}.
}
\seealso{
  \code{\link{summary.pcaviz}}, \code{\link{pcaviz-plots}},
  \code{\link[stats]{prcomp}}, \code{\link[stats]{princomp}},
  \code{\link[rsvd]{rpca}}, \link{pcaviz-transformations}
}
\examples{

# Create pcaviz object from prcomp output.
data(iris)
a1 <- pcaviz(prcomp(iris[1:4]),dat = iris)

# Create pcaviz object from princomp output.
a2 <- pcaviz(princomp(iris[1:4]),dat = iris)

# An alternative way to create pcaviz object from prcomp output.
out <- prcomp(iris[1:4])
a3  <- pcaviz(x = out$x,sdev = out$sdev,var = sum(out$sdev^2),
              rotation = out$rotation,dat = iris)

\dontrun{

# Print summaries of the pcaviz objects.
summary(a1)
summary(a2)
summary(a3)

# These vignettes demonstrate various ways pcaviz objects
# can be created.
vignette("iris")
vignette("popres")
vignette("regmap")
}
}
