\name{regmap}
\docType{data}
\alias{regmap}
\title{Principal components analysis of A. thaliana RegMap data}
\description{
  Data on 1,307 \emph{Arabidopsis thaliana} accessions collected around
  the world. Data include geographic co-ordinates and results of
  applying principal components analysis (PCA) to genotypes
  at 214,051 single nucleotide polymorphisms (SNPs).
}
\usage{data(regmap)}
\format{
  \code{regmap} is a data frame with 1,307 rows and 21 columns, in which
  each row represents an \emph{A. thaliana} accession included
  in the RegMap panel. The last 11 columns (labeled \code{"PC1"} through
  \code{"PC10"}) store the projection of the samples onto the first
  10 PCs.
}
\source{
  \emph{Arabidopsis thaliana} SNP genotypes were obtained from the
  website of Joy Bergelson's lab at the University of Chicago,
  \url{http://bergelson.uchicago.edu/?page_id=790}.
}
\details{
  PCs were computed using \code{rpca} from the \pkg{rsvd} package after
  first converting the genotypes to a 1,307 x 214,051 binary matrix.
}
\references{
  Horton, M. W., Hancock, A. M., Huang, Y. S., Toomajian, C., Atwell,
  S., Auton, A., \emph{et al} (2012). Genome-wide patterns of genetic
  variation in worldwide Arabidopsis thaliana accessions from the RegMap
  panel. \emph{Nature Genetics} \bold{44}, 212–216.
}
\examples{

# View the first few data samples.
data(regmap)
head(regmap)

\dontrun{
# See the vignette for more examples.
vignette("regmap")
}
}
\keyword{datasets}
