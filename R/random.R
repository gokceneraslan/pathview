random <-
function(x, na.rm=TRUE){
if(na.rm) x=x[!is.na(x)]
sample(x, 1)
}

