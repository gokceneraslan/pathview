render.kegg.node <-
function(plot.data, cols.ts, img, same.layer=TRUE, type=c("gene","compound")[1], text.col="black", cex=0.25){
  width=ncol(img)
  height=nrow(img)
if(type=="gene"){
  if(same.layer!=T){
    rect(plot.data$x-plot.data$width/2+0.5, height-plot.data$y-plot.data$height/2+0.25,
         plot.data$x+plot.data$width/2-0.5, height-plot.data$y+plot.data$height/2-0.25,
         col=cols.ts, border=NA)
    text(plot.data$x, height-plot.data$y, labels = as.character(plot.data$labels),
         cex = cex, col = text.col)
    return(invisible(1))
  } else{
    col.rgb=col2rgb(cols.ts)/255
    img2=img
    pidx=cbind(ceiling(plot.data$x-plot.data$width/2)+1,
      floor(plot.data$x+plot.data$width/2)+1,
      ceiling(plot.data$y-plot.data$height/2)+1,
      floor(plot.data$y+plot.data$height/2)+1)
    nn=nrow(plot.data)
    for(i in 1:nn){
      node.rgb=img2[pidx[i,3]:pidx[i,4],pidx[i,1]:pidx[i,2], 1:3]
      node.rgb.sum=apply(node.rgb,c(1,2), sum)
      blk.ind=which(node.rgb.sum==0,arr.ind=T)
      node.rgb=array(col.rgb[,i],dim(node.rgb)[3:1])
      node.rgb=aperm(node.rgb, 3:1)
      for(j in 1:3) node.rgb[cbind(blk.ind,j)]=0
      img2[pidx[i,3]:pidx[i,4],pidx[i,1]:pidx[i,2], 1:3]=node.rgb
    }
    return(img2)
  }
} else if(type=="compound"){
  if(same.layer!=T){
    circles(plot.data$x, height-plot.data$y, plot.data$width[1], cols=cols.ts, border=NULL, n=16)
    return(invisible(1))
  } else{
    col.rgb=col2rgb(cols.ts)/255
    blk=c(0,0,0)
    img2=img
    w=ncol(img)
    h=nrow(img)
    cidx=rep(1:w, each=h)
    ridx=rep(1:h, w)
    nn=nrow(plot.data)
    pidx=lapply(1:nn, function(i){
      ii=which((cidx-plot.data$x[i])^2+(ridx-plot.data$y[i])^2<(plot.data$width[i])^2)
      imat=cbind(cbind(ridx, cidx)[rep(ii,each=3),],1:3)
      imat[,1:2]=imat[,1:2]+1
      ib=which(abs((cidx-plot.data$x[i])^2+(ridx-plot.data$y[i])^2-(plot.data$width[i])^2)<=8)
      ibmat=cbind(cbind(ridx, cidx)[rep(ib,each=3),],1:3)
      ibmat[,1:2]=ibmat[,1:2]+1
      return(list(fill=imat,border=ibmat))
    })
    for(i in 1:nn){
      img2[pidx[[i]]$fill]=col.rgb[,i]
      img2[pidx[[i]]$border]=blk
    }
    return(img2)
  }
} else stop("unrecognized node type!")
}

