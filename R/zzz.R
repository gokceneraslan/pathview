.onLoad <- function(libname, pkgname) {
  pnames=rownames(installed.packages())
  if("pathview" %in% pnames){
    data(gene.idtype.list, package ="pathview")
    data(cpd.simtypes, package ="pathview")
  }
disclaimer="##############################################################################\nPathview is an open source software package distributed under GNU General Public License version 3 (GPLv3). Details of GPLv3 is available at http://www.gnu.org/licenses/gpl-3.0.html.\n\nThe pathview downloads and uses KEGG data. Academic users may freely use the KEGG website at http://www.kegg.jp/ or its mirror site at GenomeNet http://www.genome.jp/kegg/. Academic users may also freely link to the KEGG website. Non-academic users may use the KEGG website as end users for non-commercial purposes, but any other use requires a license agreement (details at http://www.kegg.jp/kegg/legal.html).\n##############################################################################\n\n"
packageStartupMessage(wordwrap(disclaimer, 80))
}
