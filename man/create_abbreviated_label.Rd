\name{create_abbreviated_label}
\alias{create_abbreviated_label}
\title{Generate abbreviated strings}
\description{Generate abbreviated strings, allowing for both a
customized abbreviation (e.g., specific to country names) and a generic
abbreviation.}
\usage{
create_abbreviated_label(x, abbrv.generic = abbreviate, abbrv.custom = NULL)
}
\arguments{
  \item{x}{A factor or character vector to be abbreviated.}
  
  \item{abbrv.generic}{Function accepting a single argument that returns
  a character vector containing abbreviations for the character strings
  in its first (and only) argument.}

  \item{abbrv.custom}{Takes the same form as the \code{abbrv.genetic}
    function, but can also return \code{NA} for any character strings
    that are unrecognized.}
}
\value{
  \code{create_abbreviated_label} returns either a character vector with
  abbreviated strings, or a factor in which the levels are abbreviated.
}
\details{
  \code{create_abbreviated_label} first applies the custom abbreviation,
  then applies the generic abbreviation to any entries that are
  unrecognized by the custom abbreviation (as indicated by \code{NA}).
}
\seealso{
  \code{\link[base]{abbreviate}}
}
\examples{

  # Define custom abbreviations for some Canadian provinces and
  # territories.
  canada.post.abbrv <- function (x)
    unlist(list(Quebec  = "QC",Ontario = "ON",Manitoba = "MB",
                Alberta = "AB",Saskatchewan = "SK",Yukon = "YT"))[x]

  # Function create_abbreviated_label can be applied to a factor or a
  # character vector.
  x <- c("Quebec","Alberta","Ontario","Alberta","Yukon")
  create_abbreviated_label(x,abbrv.custom = canada.post.abbrv)
  create_abbreviated_label(factor(x),abbrv.custom = canada.post.abbrv)
}
