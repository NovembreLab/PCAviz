\name{pcaviz-transformations}
\title{Transformations of principal components}
\alias{pcaviz-transformations}
\alias{pcaviz-transformation}
\alias{matrix-transformations}
\alias{matrix-transformation}
\alias{pcaviz_rotate}
\alias{pcaviz_reflect}
\alias{pcaviz_scale}
\alias{pcaviz_translate}
\alias{pcaviz_transform2d}
\alias{matrix_rotate}
\alias{matrix_reflect}
\alias{matrix_scale}
\alias{matrix_translate}
\alias{matrix_transform2d}
\description{
  These functions define basic transformations of points in Euclidean
  space. This includes principal components (PCs) as a special
  case. Transformations include rotations, reflections, scaling and
  translation. For plotting PCs, these functions will typically be used
  to apply transformations along two dimensions (\emph{i.e.}, the
  Cartesian plane), but more general transformations on n-dimensional
  space are also possible. In all cases, points are stored as row
  vectors in a matrix, or in a data frame within a pcaviz object.
}
\usage{
pcaviz_rotate(x, angle, dims, units = c("degrees","radians"))

pcaviz_reflect(x, dims)

pcaviz_scale(x, scale, dims)

pcaviz_translate(x, a, dims)

pcaviz_transform2d(x, dims, angle = 0,reflect.x = FALSE,
                   reflect.y = FALSE,scale = c(1,1),a = c(0,0),
                   units = c("degrees","radians"))

matrix_rotate(X, angle, dims = 1:2, units = c("degrees","radians"))

matrix_reflect(X, dims)

matrix_scale(X, scale, dims)

matrix_translate(X, a, dims)

matrix_transform2d(X, dims = 1:2, angle = 0,reflect.x = FALSE,
                   reflect.y = FALSE,scale = c(1,1),a = c(0,0),
                   units = c("degrees","radians"))
}
\arguments{
  \item{x}{An object of class "pcaviz".}

  \item{X}{A matrix in which the points are stored as rows.}
  
  \item{dims}{The transformation is applied to these dimensions
    (\emph{i.e.}, column indices of \code{X}.) The dimensions must be
    specified as \code{character} for the \code{pcaviz_} functions, and
    \code{numeric} or \code{character} for the \code{matrix_}
    functions. For \code{pcaviz_rotate}, \code{pcaviz_transform2d},
    \code{matrix_rotate} and \code{matrix_transform2d}, this must be a
    vector of length 2. By default, \code{pcaviz_rotate} and
    \code{pcaviz_transform2d} apply the transformation to the first two
    PCs. See \code{\link{Extract}} for more details on selecting matrix
    columns with numeric and character indices.}
  
  \item{angle}{A scalar specifying the angle of rotation. Positive
    angles specify counter-clockwise rotations.}
  
  \item{units}{The angle of rotation is specified in these units.}
  
  \item{scale}{A numeric vector specifying the scaling factor for each
    selected dimension; must be a vector of the same length as
    \code{dims}.}
  
  \item{a}{A numeric vector specifying the translation of the selected
  dimensions; must be a vector of the same length as \code{dims}.}
  
  \item{reflect.x}{If TRUE, reflect first dimension, \code{dims[1]}.}
  
  \item{reflect.y}{If TRUE, reflect second dimension, \code{dims[2]}.}
}
\value{
  A pcaviz object with the transformed PCs, or a matrix with the same
  dimensions as \code{X} giving the transformed points.
}
\details{

\code{matrix_rotate} applies a counter-clockwise plane rotation to row
vectors in \code{X} when \code{angle} is positive. \code{pcaviz_rotate}
applies a counter-clockwise plane rotation to selected PCs in \code{x}.

\code{matrix_reflect} reflects points about the origin, individually
along each selected dimension. \code{pcaviz_reflect} reflects selected
PCs.

\code{matrix_scale} scales the selected dimensions by \code{scale}, and
\code{pcaviz_scale} scales the selected PCs by \code{scale}. Note that
it is also possible to use a negative scaling factor, which is
equivalent to a positive scaling combined with a reflection about the
origin.

\code{matrix_translate} translates row vectors in \code{X} by \code{a}
within the selected dimensions. \code{pcaviz_translate} translates
selected PCs by {a}.

\code{matrix_transform2d} and \code{pcaviz_transform2d} are convenience
functions that that apply a sequence of 2-d transformations in the
following order: (1) counter-clockwise rotation, (2) reflection, (3)
scaling and (4) translation. Note that exactly two dimensions (PCs) must
be selected.
}
\examples{

# A matrix storing 4 points in 2-d Euclidean space (i.e., the
# Cartesian plane).
X <- rbind(c(0,0),
           c(1,0),
           c(0,1),
           c(1,1))

# Rotate the points 45 degrees counter-clockwise.
matrix_rotate(X,45)

# Reflect the points along the first dimension; i.e., X axis
# in the Cartesian plane about the X=0 line.
matrix_reflect(X,1)

# Translate all points by (-1,2).
matrix_translate(X,a = c(-1,2),1:2)

# Apply multiple transformations to the points.
matrix_transform2d(X,angle = -90,reflect.x = TRUE,reflect.y = TRUE,
                   scale = c(2,5),a = c(-1,2))

# Create pcaviz object from prcomp output.
data(iris)
a <- pcaviz(prcomp(iris[1:4]),dat = iris)

# Rescale the first two PCs.
a <- pcaviz_scale(a,scale = c(0.5,1.25),dims = c("PC1","PC2"))

# Rotate the first two PCs 15 degrees counter-clockwise.
a <- pcaviz_rotate(a,angle = 15,dims = c("PC1","PC2"))
}
