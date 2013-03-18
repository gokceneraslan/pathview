\name{korg}
\alias{korg}
\alias{bods}
\docType{data}
\title{
  Mapping data on KEGG species code and corresponding Bioconductor gene
  annotation package
}
\description{
  Data on mapping between KEGG code, scientific name and common name of species or
  organisms included in KEGG databases, and gene annotation package names in Bioconductor
}
\usage{
data(korg)
data(bods)
}
\format{
  korg is a character matrix of 2299 rows and 3 columns on the mapping
  between KEGG code, scientific name and common name of species or
  organisms included in KEGG databases.
  bods is a character matrix of 19 rows and 3 columns on the mapping
  between gene annotation package names in Bioconductor, common name and
  KEGG code of most common research species.
}
\source{
  http://www.genome.jp/kegg-bin/get_htext?br08601.keg
  http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
}
\examples{
data(korg)
data(bods)
head(korg)
head(bods)
}
\keyword{datasets}
