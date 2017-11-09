\name{EUNAM_geno}
\alias{EUNAM_geno}
\docType{data}
\title{
  Reduced genotype data maize EU-NAM Dent population
}
\description{
  Selection of markers and genotypes from the maize EU nested association mapping
  (NAM) population (Bauer et al., 2013).
}
\usage{data(EUNAM_geno)}
\format{
  \code{data.frame}
}
\details{
  Sample of the marker matrix of the EU-NAM population. The selection correspond
  to 2182 markers coming from chromosomes 2 and 3 also present in
  \code{\link{EUNAM_map}}.The 300 genotypes come from six crosses:
    (F353 x B73), (F353 x D06), (F353 x D09), (F353 x EC169), (F353 x F252) and
  (F353 x F618). The selection contains 50 per cross. These genotypes correspond
  to the selected phenotypic values in \code{\link{EUNAM_pheno}}. The data are
  available on http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50558.
}
\source{
  http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50558
}
\references{
  Bauer, E., Falque, M., Walter, H., Bauland, C., Camisan, C., Campo, L., ... &
    Altmann, T. (2013). Intraspecific variation of recombination rate in maize.
  Genome biology, 14(9), R103.

}
\examples{
  data(EUNAM_geno)
}

\seealso{
  \code{\link{EUNAM_pheno}}, \code{\link{EUNAM_map}}
}

\keyword{datasets}
