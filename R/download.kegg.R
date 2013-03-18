download.kegg <-
function (pathway.id = "00010", species = "hsa", kegg.dir = ".")
  {
      npath=length(pathway.id)
      if(species!="ko") species=kegg.species.code(species, na.rm=T)
      nspec=length(species)

      if(npath!=1 | nspec!=1) {
        species=rep(species, npath)
        pathway.id=rep(pathway.id, each=nspec)
      }
      pathway.id <- paste(species, pathway.id, sep = "")
      uidx=!duplicated(pathway.id)
      pathway.id=pathway.id[uidx]
      species=species[uidx]
      npath=length(pathway.id)
      
      xml.fnames=paste(pathway.id, ".xml", sep="")
      png.fnames=paste(pathway.id, ".png", sep="")
      xml.fmt="http://www.genome.jp/kegg-bin/download?entry=%s&format=kgml"
      png.fmt="http://www.genome.jp/kegg/pathway/%s/%s"
      success=rep(NA, npath)
      names(success)=pathway.id
      warn.fmt="Download of %s xml and png files failed!"
      
      for (i in 1:npath) {
      msg=sprintf("Downloading xml files for %s, %d/%d pathways..", pathway.id[i], i, length(pathway.id))
      print(msg)
      xml.url=sprintf(xml.fmt,  pathway.id[i])
      xml.target=sprintf("%s/%s", kegg.dir, xml.fnames[i])
      xml.status=try(download.file(xml.url, xml.target, quiet=T), silent=T)

      msg=sprintf("Downloading png files for %s, %d/%d pathways..", pathway.id[i], i, length(pathway.id))
      print(msg)
      png.url=sprintf(png.fmt, species[i], png.fnames[i])
      png.target=sprintf("%s/%s", kegg.dir, png.fnames[i])
      png.status=try(download.file(png.url, png.target, quiet=T, mode="wb"), silent=T)

      success[i]=ifelse(png.status==0, "succeed", "failed")
      if(class(png.status)=="try-error"){
        warn.msg=sprintf(warn.fmt, pathway.id[i])
        message(warn.msg)
      }
    }
      return(success)
  }

