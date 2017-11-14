\name{popres}
\docType{data}
\alias{popres}
\title{Principal components analysis of genetic variation in POPRES sample}
\description{
  Data on 1,387 European subjects from the POPRES panel (Nelson et al,
  2008), and results from principal components analysis (PCA) of
  genotypes at 197,146 single nucleotide polymorphisms (SNPs).
}
\usage{data(popres)}
\format{\code{popres$x} is a data frame with 1,387 rows and 10 columns,
  in which each row corresponds to an individual:
  
  \describe{
    \item{id}{sample id}
    \item{country}{country of origin}
    \item{abbrv}{abbreviated country label}
    \item{color}{mapping from countries to colors}
    \item{longitude}{geographic longitude, in degrees}
    \item{latitude}{geographic latitude, in degrees}
    \item{PC1--4}{projection of sample onto first 4 PCs}
  }

  \code{popres$rotation} is a data frame containing the eigenvectors
  ("loadings") computed from the genotype data at a subset of the SNPs.
  The first two columns give the chromosome number (1--22) and base-pair
  position of each SNP.
}
\details{
  Note that some of the country designations do not reflect current
  political boundaries.

  Also note that "country" and "abbrv" columns do not correspond
  exactly; for example, abbreviations "GB" (Great Britain) and "Sct"
  (Scotland) are both assigned to samples from the United Kingdom.

  Geographic coordinates are the same for each country, so
  all samples with the same country of origin are assigned the same
  latitude and longitude.
}
\source{Data were retrieved and compiled from 
\url{http://github.com/NovembreLab/Novembre_etal_2008_misc}. To reduce
file size, only the top 4 PCs were retained.
}
\references{
  Nelson, M. R., Bryc, K., King, K. S., Indap, A., Boyko, A. R.,
  Novembre, J., \emph{et al} (2008). The population reference sample,
  POPRES: a resource for population, disease, and pharmacological
  genetics research. \emph{American Journal of Human Genetics} \bold{83},
  347--358.
  
  Novembre, J., Johnson, T., Bryc, K., Kutalik, Z., Boyko, A. R., Auton,
  A., Indap, A., King, K. S., Bergmann, S., Nelson, M. R., Stephens, M.
  and Bustamante, C. D. (2008) Genes mirror geography within Europe.
  \emph{Nature} \bold{456}, 274-–274.
}
\examples{

# View the first few samples and genetic markers.
data(popres)
head(popres$x[c(1:3,7:8)])
head(popres$rotation)

\dontrun{
# See the vignette for more examples.
vignette("popres")
}
}
\keyword{datasets}
