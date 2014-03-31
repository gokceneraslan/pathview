eg2id <-
function(eg, category=gene.idtype.list[1:2], org="Hs", pkg.name=NULL){
category=toupper(category)
  if(is.null(pkg.name)) pkg.name=paste("org", org, "eg.db", sep=".")
  pkg.on=require(pkg.name, character.only = TRUE)
    if(!pkg.on) {
          source("http://bioconductor.org/biocLite.R")
          biocLite(pkg.name, suppressUpdates =TRUE)
          pkg.on=require(pkg.name, character.only = TRUE)
          if(!pkg.on) stop(paste("Fail to install/load gene annotation package ", pkg.name, "!",  sep=""))
        }
  pkg.name=gsub("[.]db", "", pkg.name)
    data(gene.idtype.list)
  idi=category %in% gene.idtype.list
  if(sum(idi)==0) stop("No correct target ID type!")
  category=category[idi]
  
  EntrezGene=eg
  annot=cbind(EntrezGene)
  for(ci in category){
    ids=character(length(eg))
    bimap=eval(as.name(paste(pkg.name, ci, sep="")))
    mapped=eg %in% mappedkeys(bimap)
    if(sum(mapped)>0) {
      map.res=AnnotationDbi::mget(eg[mapped], bimap)
#      mr.len=sapply(map.res, length)
#      print(c(ci,mean(mr.len)))
      ids[mapped]=sapply(map.res, paste, sep="", collapse="; ")
    }
    annot=cbind(annot, ids)
  }
  colnames(annot)[-1]=category
  return(annot)
}

