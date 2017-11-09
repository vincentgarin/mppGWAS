\name{EUNAM_pheno}
\alias{EUNAM_pheno}
\docType{data}
\title{
  Reduced phenotype data from Maize EU-NAM population
}
\description{
  Reduced phenotype data from the Maize EU nested association mapping (NAM)
  population (Bauer et al., 2013).
}
\usage{data(EUNAM_pheno)}
\format{
  data.frame
}
\details{
  Dry matter content (DMC) scores for the 300 offspring genotypes of the marker
  matrix \code{\link{EUNAM_geno}}. The values represent the adjusted means
  computed from plot data using the same procedure as Lehermeier et al. (2014).
  The original plot data are available here
  http://www.genetics.org/content/198/1/3/suppl/DC1.

}
\source{
  http://www.genetics.org/content/198/1/3/suppl/DC1
}
\references{
  Bauer, E., Falque, M., Walter, H., Bauland, C., Camisan, C., Campo, L., ... &
    Altmann, T. (2013). Intraspecific variation of recombination rate in maize.
  Genome biology, 14(9), R103.

  Lehermeier, C., Krämer, N., Bauer, E., Bauland, C., Camisan, C., Campo, L., ...
  & Moreau, L. (2014). Usefulness of multiparental populations of maize
  (Zea mays L.) for genome-based prediction. Genetics, 198(1), 3-16.

}
\examples{
  data(EUNAM_pheno)
}

\seealso{
  \code{\link{EUNAM_geno}}, \code{\link{EUNAM_map}}
}

\keyword{datasets}
