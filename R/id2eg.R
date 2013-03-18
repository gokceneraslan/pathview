id2eg <-
function(ids, category=gene.idtype.list[1], org="Hs", pkg.name=NULL){
  category=toupper(category)
  if(is.null(pkg.name)) pkg.name=paste("org", org, "eg.db", sep=".")
  pkg.on=require(pkg.name, character.only = TRUE)
    if(!pkg.on) {
          source("http://bioconductor.org/biocLite.R")
          biocLite(pkg.name)
          pkg.on=require(pkg.name, character.only = TRUE)
          if(!pkg.on) stop(paste("Fail to install/load gene annotation package ", pkg.name, "!",  sep=""))
        }
  pkg.name=gsub("[.]db", "", pkg.name)
  data(gene.idtype.list)  
  idi=category %in% gene.idtype.list
  if(!idi) stop("Invalid source ID type!")
  
  annot=cbind(ids)
  egs=character(length(ids))
  rev.map=paste(pkg.name, category, 2, "EG", sep="")
  if(exists(rev.map))  {
    bimap=eval(as.name(rev.map))
  } else {
    for.map=paste(pkg.name, category, sep="")
    for.map=eval(as.name(for.map))
    if(class(for.map)!="AnnDbBimap") for.map=as(for.map,"AnnDbBimap")
    bimap=revmap(for.map)
  }
  
  mapped=ids %in% mappedkeys(bimap)
  if(sum(mapped)>0) {
    map.res=AnnotationDbi::mget(ids[mapped], bimap)
    mr.len=sapply(map.res, length)
#    print(c(category,mean(mr.len)))
    if(any(mr.len>1)){
    map.res[mr.len>1]=lapply(map.res[mr.len>1], unique)
    mr.len=sapply(map.res, length)
  }
#    print(c(category,mean(mr.len)))
    egs[mapped]=sapply(map.res, paste, sep="", collapse="; ")
  }

  annot=cbind(annot, EntrezGene=egs)

  colnames(annot)[1]=category
  return(annot)
}

