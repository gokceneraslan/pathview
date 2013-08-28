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

      if(!species %in% c("hsa", "ko")){
      if(!id.type %in% c("ENTREZ","KEGG")){
        msg=sprintf(msg.fmt, id.type, species, mol.type)
        stop(msg)
      }
      if(is.na(species.data["ncbi.geneid"])){
        if(!is.na(species.data["kegg.geneid"])){
          msg.fmt3="Only native KEGG gene ID is supported for species \"%s\"!"
          msg=sprintf(msg.fmt3, species)
          message(msg)
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

    } else if(species=="hsa"){
      if(id.type=="KEGG") id.type="ENTREZ"
      data(gene.idtype.list)

      if(!id.type %in% c("ENTREZ",gene.idtype.list)){
        msg=sprintf(msg.fmt, id.type, species, mol.type)
        stop(msg)
      }
      pkg.name="org.Hs.eg.db"
      require(pkg.name, character.only = TRUE)
      pkg.name = gsub("[.]db", "", pkg.name)
      if(id.type=="ENTREZ"){
        bi.map = paste(pkg.name, "SYMBOL", sep = "")
        bimap = eval(as.name(bi.map))
      } else {
        rev.map = paste(pkg.name, id.type, 2, "EG", sep = "")
        if (exists(rev.map)) {
          bimap = eval(as.name(rev.map))
        }
        else {
          for.map = paste(pkg.name, id.type, sep = "")
          for.map = eval(as.name(for.map))
          if (class(for.map) != "AnnDbBimap")
            for.map = as(for.map, "AnnDbBimap")
          bimap = revmap(for.map)
        }
      }
      all.mn=keys(bimap)
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
    message(msg)
    nmol=nuids
  }
  sel.mn=sample(all.mn, nmol)
  if(discrete) return(sel.mn)
  sel.mn.data=matrix(rnorm(nmol*nexp), ncol=nexp)
  rownames(sel.mn.data)=sel.mn
  colnames(sel.mn.data)=paste("exp", 1:nexp, sep="")
  return(sel.mn.data[, 1:nexp])
}

