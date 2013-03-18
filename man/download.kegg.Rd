\name{download.kegg}
\alias{download.kegg}

\title{
Download KEGG pathway graphs and associated KGML data
}
\description{
This is the downloader function for KEGG pathways, automatically
download graph images and associated KGML data.
}
\usage{
download.kegg(pathway.id = "00010", species = "hsa", kegg.dir = ".")
}

\arguments{
  \item{pathway.id}{
character, 5-digit KEGG pathway IDs. Default pathway.id="00010".
}
  \item{species}{
character, either the KEGG code, scientific name or the common name of
the target species. When KEGG ortholog pathway is considered,
species="ko". Default species="hsa", it is equivalent to use either
"Homo sapiens" (scientific name) or "human" (common name). 
}
  \item{kegg.dir}{
character, the directory of KEGG pathway data file (.xml) and image file
(.png). Default kegg.dir="." (current working directory).
}
}
\details{
  Species can be specified as either kegg code, scientific name or the
  common name. Scientific name and the common name are always mapped to
  kegg code first.
  Length of species should be either 1 or the same as pathway.id, if
  not, the same set of pathway.id will be applied to all species.
}
\value{
a named character vector, either "succeed" or "failed", indicating the
download  status of corresponding pathways.
}
\references{
  Luo, W. et al. Pathview: a pathway based data integration and visualization
  tool. submitted, 2013.
}
\author{
  Weijun Luo <luo_weijun@yahoo.com>
}


\seealso{
  \code{\link{pathview}} the main function,
  \code{\link{node.info}} the parser,
}
\examples{
data(demo.paths)
sel.2paths=demo.paths$sel.paths[1:2]
download.kegg(pathway.id = sel.2paths, species = "hsa")
#pathway files should be downloaded into current working directory
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
