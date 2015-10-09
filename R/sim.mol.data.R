sim.mol.data=function(mol.type=c("gene","gene.ko","cpd")[1], id.type=NULL, species="hsa", discrete=FALSE, nmol=1000, nexp=1, rand.seed=100)
{
  msg.fmt="\"%s\" is not a good \"%s\" \"%s\" ID type for simulation!"
  msg.fmt2="\"%s\" has only %i unique IDs!"
  set.seed(rand.seed)

  if(species!="ko"){
    species.data=kegg.species.code(species, na.rm=T, code.only=FALSE)
    species=species.data["kegg.code"]
  } else if(mol.type=="gene") mol.type="gene.ko"
  
  
  if(mol.type=="gene"){
    if(is.null(id.type)) id.type="KEGG"
    id.type=toupper(id.type)

    data(bods)
    data(gene.idtype.bods)
    org19=bods[,"kegg code"]
    
    if(!species %in% c(org19, "ko")){
      if(!id.type %in% c("ENTREZ","KEGG")){
        msg=sprintf(msg.fmt, id.type, species, mol.type)
        stop(msg)
      }
      if(is.na(species.data["ncbi.geneid"])){
        if(!is.na(species.data["kegg.geneid"])){
          msg.fmt3="Only native KEGG gene ID is supported for species \"%s\"!"
          msg=sprintf(msg.fmt3, species)
          message("Note: ", msg)
        } else{
          msg.fmt3="Simulation is not supported for species \"%s\"!"
          msg=sprintf(msg.fmt3, species)
          stop(msg)
        }
      }
      gid.map=keggConv("ncbi-geneid",species)
      if(id.type=="KEGG") {
        all.mn=gsub(paste(species, ":", sep=""), "", names(gid.map))
      } else all.mn=gsub("ncbi-geneid:", "", gid.map)

    } else if(species %in% org19){
      if(id.type=="ENTREZ") id.type="ENTREZID"
      if(id.type=="KEGG") {
        gid.map=keggConv("ncbi-geneid",species)        
        all.mn=gsub(paste(species, ":", sep=""), "", names(gid.map))
      } else if(id.type %in% gene.idtype.bods[[species]]){
        idx=which(bods[,3]==species)
        pkg.name=bods[idx,1]
        pkg.on=requireNamespace(pkg.name)
        if(!pkg.on) {
          source("http://bioconductor.org/biocLite.R")
          biocLite(pkg.name, suppressUpdates =TRUE)
          pkg.on=requireNamespace(pkg.name)
          if(!pkg.on) stop(paste("Fail to install/load gene annotation package ", pkg.name, "!",  sep=""))
        }
        
        db.obj <- eval(parse(text=paste0(pkg.name, "::", pkg.name)))
        all.mn <-keys(db.obj, keytype=id.type)
      } else stop("Wrong gene ID type!")
    }            
  } else if(mol.type=="cpd"){
    data(cpd.accs)
    data(cpd.simtypes)
    data(rn.list)
    accn=cpd.accs$ACCESSION_NUMBER
    if(is.null(id.type)) id.type="KEGG COMPOUND accession"
    if(!id.type %in% cpd.simtypes){
      msg=sprintf(msg.fmt, id.type, mol.type)
      stop(msg)
    }
    all.mn=unique(as.character(accn[rn.list[[id.type]]]))
  } else if(mol.type=="gene.ko"){
    data(ko.ids)
    all.mn=ko.ids
  } else stop("Invalid mol.type!") 

  nuids=length(all.mn)
  if(nmol>nuids){
    msg=sprintf(msg.fmt2, id.type, nuids)
    message("Note: ", msg)
    nmol=nuids
  }
  sel.mn=sample(all.mn, nmol)
  if(discrete) return(sel.mn)
  sel.mn.data=matrix(rnorm(nmol*nexp), ncol=nexp)
  rownames(sel.mn.data)=sel.mn
  colnames(sel.mn.data)=paste("exp", 1:nexp, sep="")
  return(sel.mn.data[, 1:nexp])
}

