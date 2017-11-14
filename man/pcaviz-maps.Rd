\name{pcaviz-maps}
\title{PCAviz maps}
\alias{pcaviz-maps}
\alias{pcaviz-map}
\alias{overlay_map_world}
\alias{overlay_map_europe}
\alias{overlay_map_usa}
\description{"Overlays" for adding maps to existing ggplots.}
\usage{
  overlay_map_world(g)
  overlay_map_europe(g)
  overlay_map_usa(g)
}
\arguments{
  \item{g}{A "ggplot" object.}
}
\value{
  All these functions return a \code{ggplot} object.
}
\details{
  These functions are all examples of PCAviz plot "overlays"; that
  is, functions that can be assigned to input argument \code{overlay}
  in \code{plot(x)}, where \code{x} is a \code{pcaviz} object; see
  \link{pcaviz-plots} for more details.
}
\seealso{
  \link{pcaviz-maps},
  \code{\link[ggplot2]{map_data}},
  \code{\link[maps]{world}},
  \code{\link[maps]{usa}}
}
\examples{
\dontrun{

# See the POPRES vignette for some examples of how these functions are
# used to add maps to PCAviz plots.
vignette("popres")
}
}
