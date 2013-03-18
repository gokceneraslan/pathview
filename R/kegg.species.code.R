kegg.species.code <-
function(species="hsa", na.rm=FALSE){
      nspec=length(species)
      data(korg)
      ridx=match(species, korg) %% nrow(korg)
      nai=is.na(ridx)
      if(sum(nai)>0) {
        na.msg=sprintf("Unknown species '%s'!", paste(species[nai], sep="", collapse="', '"))
        message(na.msg)
        }
      if(sum(nai)==nspec) {
        stop.msg="All species are invalid!"
        stop(stop.msg)
        }
      if(na.rm) ridx=ridx[!nai]
      species=korg[ridx,1]
      return(species)
    }

