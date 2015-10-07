kegg.species.code <-
function(species="hsa", na.rm=FALSE, code.only=TRUE){
      nspec=length(species)
      if(!exists("korg")) data(korg)

      ridx=match(species, korg[,1:3]) %% nrow(korg)
      nai=is.na(ridx)
      if(sum(nai)>0) {
        na.msg=sprintf("Unknown species '%s'!", paste(species[nai], sep="", collapse="', '"))
        message("Note: ", na.msg)
        }
      if(sum(nai)==nspec) {
        stop.msg="All species are invalid!"
        stop(stop.msg)
        }
      if(any(ridx[!nai]==0)) ridx[!nai & ridx==0]=nrow(korg)
      if(na.rm) ridx=ridx[!nai]
      if(code.only) coln=1 else coln=c(1,4:6)
      species.info=korg[ridx,coln]
      return(species.info)
    }

