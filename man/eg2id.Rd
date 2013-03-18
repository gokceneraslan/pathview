\name{eg2id}
\alias{eg2id}
\alias{id2eg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mapping between Entrez Gene and other IDs
}
\description{
These auxillary gene ID mappers connect Entrez Gene ID to external gene,
transcript or protein IDs or vise versa. 
}
\usage{
eg2id(eg, category = gene.idtype.list[1:2], org = "Hs", pkg.name = NULL)
id2eg(ids, category = gene.idtype.list[1], org = "Hs", pkg.name = NULL)
}
\arguments{
  \item{eg}{
character, input Entrez Gene IDs.
}
  \item{ids}{
character, input gene/transcript/protein IDs to be converted to Entrez
Gene IDs.
}
  \item{category}{
character, for \code{eg2id} the output ID types to map from Entrez Gene,
d to be c("SYMBOL", "GENENAME"); for \code{id2eg}, the input ID type  to
be mapped to Entrez Gene, default to be "SYMBOL".
}
  \item{org}{
character, the two-letter abbreviation of organism name, used to
determine the gene annotation package. Default org="Hs". Only effective when
pkg.name is not NULL.
}
  \item{pkg.name}{
character, name of the gene annotation package. This package should be
one of the standard annotation packages from Bioconductor, such as
"org.Hs.eg.db". Check \code{data(bods); bods} for a full list of standard
annotation packages. You may also use your custom annotation package
built with AnnotationDbi, the Bioconductor Annotation Database
Interface. Default pkg.name=NULL, hence argument \code{org} should be
specified.
}
}
\details{
KEGG uses Entrez Gene ID as its standard gene ID. Therefore, all gene
data need to be mapped to Entrez Genes when working with KEGG
pathways. Function \code{id2eg} does this mapping. On the other hand, we
frequently want to check or show gene symbols or full names instead of
the less informative Entrez Gene ID when working with KEGG gene nodes,
Function \code{eg2id} does this reverse mapping.
These functions are
written as part of the Pathview mapper module, they are equally useful
for other gene ID or data mapping tasks.
The use of these functions depends on gene annotation packages like
"org.Hs.eg.db", which are Bioconductor standard. IFf no such packages not available for
your interesting organisms, you may build one with Bioconductor
AnnotationDbi package.
}
\value{
a 2- or multi-column character matrix recording the mapping between input IDs to
the target ID type(s).
}
\references{
  Luo, W. et al. Pathview: a pathway based data integration and visualization
  tool. submitted, 2013.
}
\author{
  Weijun Luo <luo_weijun@yahoo.com>
}
\seealso{
  \code{\link{cpd2kegg}} etc the auxillary compound ID mappers,
  \code{\link{mol.sum}} the auxillary molecular data mapper,
  \code{\link{node.map}} the node data mapper function.
}
\examples{
data(gene.idtype.list)
#generate simulated gene data named with non-KEGG/Entrez gene IDs
gene.ensprot <- sim.mol.data(mol.type = "gene", id.type = gene.idtype.list[4], 
    nmol = 50000)
#construct map between non-KEGG ID and KEGG ID (Entrez gene)
id.map.ensprot <- id2eg(ids = names(gene.ensprot), 
    category = gene.idtype.list[4], org = "Hs")
#Map molecular data onto Entrez Gene IDs
gene.entrez <- mol.sum(mol.data = gene.ensprot, id.map = id.map.ensprot)
#check the results
head(gene.ensprot)
head(id.map.ensprot)
head(gene.entrez)

#map Entrez Gene to Gene Symbol and Name
eg.symbname=eg2id(eg=id.map.ensprot[,2])
#entries with more than 1 Entrez Genes are not mapped
head(eg.symbname)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
