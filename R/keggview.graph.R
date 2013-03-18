keggview.graph <-function(
                           plot.data.gene=NULL,
                           plot.data.cpd=NULL,
                           cols.ts.gene=NULL,
                           cols.ts.cpd=NULL,
                           node.data,
                           path.graph,
                           pathway.name,
                           out.suffix="pathview",
                           pdf.size=c(7,7),


                          same.layer=TRUE,
                          match.data=TRUE,
                           rankdir=c("LR","TB")[1],
                           is.signal=TRUE,
                           split.group=F,
                           afactor=1,

                           text.width=15, #k
                           cex=0.5,
                           map.cpdname=FALSE, #k
                           cpd.lab.offset=1.0,


                          discrete=list(gene=FALSE, cpd=FALSE),
                          limit=list(gene=1, cpd=1),
                         bins=list(gene=10, cpd=10),
                         both.dirs=list(gene=T, cpd=T),
                          low = list(gene = "green", cpd = "blue"),
                         mid = list(gene = "gray", cpd = "gray"),
                         high = list(gene = "red", cpd = "yellow"),
                           na.col="transparent",
         
                         new.signature=TRUE,
                         plot.col.key=TRUE,
                         key.align="x",
                         key.pos="topright",
                         sign.pos="bottomright",#g
                           ...){
  
  gR1=path.graph

                                        #group nodes mapping and merge
  grp.idx=node.data$size>1
  if(sum(grp.idx)>0 & !split.group){
    sub2grp=cbind(unlist(node.data$component[grp.idx], use.names=F), rep(names(grp.idx)[grp.idx], node.data$size[grp.idx]))
    du.idx=duplicated(sub2grp[,1])
    if(sum(du.idx)>0){
      du.rn=sub2grp[,1] %in% sub2grp[du.idx,1]
      message("warning: reconcile groups sharing member nodes!")
      print(sub2grp[du.rn,])
      du.grps=sub2grp[du.idx,]
      rn=which(du.idx)
      for(r in rn){
        comps=node.data$component[[sub2grp[r,2]]]
        comps=comps[comps!=sub2grp[r,1]]
        node.data$component[[sub2grp[r,2]]]=comps
        node.data$size[sub2grp[r,2]]=node.data$size[sub2grp[r,2]]-1
      }
      sub2grp=sub2grp[!du.idx,]
    }
    rownames(sub2grp)=sub2grp[,1]
  } else sub2grp=NULL
  
  if(sum(grp.idx)>0 & !split.group){
    for(gn in names(grp.idx)[grp.idx]){
      gR1=combineKEGGnodes(node.data$component[[gn]], gR1, gn)
    }
  } else if(split.group){
    gR1=subGraph(nodes(gR1)[node.data$size==1], gR1)
  }
  nNames=nodes(gR1)
  nSizes=node.data$size[nNames]

                                        #unconnected nodes processing
  deg=degree(gR1)
  deg=deg$inDegree+deg$outDegree
  if(is.signal & sum(deg<1)>0){
    gR2=subKEGGgraph(nNames[deg>0], gR1)
    nNames=nNames[deg>0]
    nSizes=nSizes[deg>0]
    if(!is.null(sub2grp)){
                                        #  sub.idx=!sub2grp[,1] %in% names(deg[deg<1])
      sub.idx=sub2grp[,1] %in% nNames |sub2grp[,2] %in% nNames
    } else sub.idx=0
  } else {
    gR2=gR1
    if(!is.null(sub2grp)){
      sub.idx=rep(T, nrow(sub2grp))
    } else sub.idx=0
  }

  if(length(nNames)<2){
        msg=sprintf("%s not rendered, 0 or 1 connected nodes!\nTry \"kegg.native=T\" instead!", pathway.name)
            message(msg)
            return(list())
      }
  
                                        #give up the KEGG positions, use graphviz layout


                                        #general attributes
  attrs=list()
  attrs$graph$rankdir="LR"
  attrs$node <- list(fixedsize=FALSE)

                                        #node attributes
  ntype=node.data$type[nNames]
  cpd.idx=ntype=="compound"
  map.idx=ntype=="map"
  rect.idx=!(cpd.idx|map.idx)
  nAttr=list()
  nAttr$label=rep('', length(nNames))
                                        #nAttr$shape=rep('rectangle', length(nNames))
  shapes=node.data$shape[nNames]
  if(any(cpd.idx)) shapes[cpd.idx]="ellipse"
  if(any(map.idx)) shapes[map.idx]="plaintext"
                                        #if(any(map.idx)) shapes[map.idx]="rectangle"
  nAttr$shape=shapes
  nAttr$height=.75*17/46*nSizes*afactor
  nAttr$width=rep(.75, length(nNames))*afactor
  if(any(cpd.idx)){
    nAttr$height[cpd.idx]=nAttr$height[cpd.idx]*1.5
    nAttr$width[cpd.idx]=nAttr$width[cpd.idx]*1.5
  }
  if(any(map.idx)){
    nAttr$height[map.idx]=nAttr$height[map.idx]*1.5
    nAttr$width[map.idx]=nAttr$width[map.idx]*2
  }
  nAttr<- lapply(nAttr, function(x) {names(x) <- nNames
                                     x})

  na.col=colorpanel2(1, low=na.col, high=na.col)
  fillcol=rep(na.col, length(nNames))
  names(fillcol)=nNames

  
                                        #edge attributes
  subdisplay <- subtypeDisplay(gR2)
  if(length(subdisplay)<1) eAttrs=list() else{
    na.rn=apply(subdisplay, 2, function(x) sum(is.na(x))==7)
  if(sum(na.rn)>0) subdisplay[,na.rn]=KEGGEdgeSubtype[KEGGEdgeSubtype[,1]=="others",rownames(subdisplay)]
  eLabel <- subdisplay["label", ]
  eCol <- subdisplay["color", ]
  eTextCol <- subdisplay["fontcolor", ]
  eLty <- subdisplay["style", ]
  eArrowhead <- subdisplay["arrowhead", ]
  if (ncol(subdisplay) == 1) {
    tmp <- colnames(subdisplay)[1]
    names(eLabel) <- names(eCol) <- names(eTextCol) <- tmp
    names(eLty) <- names(eArrowhead) <- tmp
  }
  eAttrs <- list(lty = eLty, col = eCol, textCol = eTextCol, 
                 label = eLabel, arrowhead = eArrowhead)
  }
  
  gR2.layout=gR2
  edgeRenderInfo(gR2.layout)=eAttrs
  layoutType=ifelse(is.signal, "dot", "neato")
  gR2.layout <- layoutGraph(gR2.layout, attrs = attrs, nodeAttrs=nAttr, layoutType=layoutType)
  nri=nodeRenderInfo(gR2.layout)
  loc=list(x=nri$nodeX, y=nri$nodeY)
  if(sum(rect.idx)>0){
  w.unit=min(nri$lWidth[rect.idx])
  h.unit=min(nri$height[rect.idx])
}
  cni=nSizes>1
  if(sum(cni)>0){
    xloc=rep(loc[[1]][cni], nSizes[cni])
    sn.y=unlist(sapply(nSizes[cni], function(x) seq(-(x-1)/2, (x-1)/2,1)),use.names =F)
    yloc=rep(loc[[2]][cni], nSizes[cni])+h.unit*sn.y
  } else xloc=yloc=NULL
  xloc.nd=c(xloc,loc[[1]][nSizes==1 & rect.idx])
  yloc.nd=c(yloc,loc[[2]][nSizes==1 & rect.idx])
                                        #labs=gsub("-", "- ", node.data$labels)
                                        #labs=sapply(labs,wordwrap,len=text.width)
                                        #labs=gsub("- ", "-", labs)
  labs=node.data$labels
  labs[nNames[map.idx]]=sapply(labs[nNames[map.idx]],wordwrap,width=text.width, break.word=F)
  labs[nNames[cpd.idx]]=sapply(labs[nNames[cpd.idx]],wordwrap,width=text.width, break.word=T)

  if(sum(cpd.idx)>0){
    ell.col=fillcol[cpd.idx]
    ell.col[ell.col==na.col]=NA
    w.e=min(nri$lWidth[cpd.idx])
    h.e=min(nri$height[cpd.idx])
    xloc.e=loc[[1]][cpd.idx]
    yloc.e=loc[[2]][cpd.idx]
  }


  cols.ts.gene=cbind(cols.ts.gene)
  cols.ts.cpd=cbind(cols.ts.cpd)
  nc.gene=max(ncol(cols.ts.gene),0)
  nc.cpd=max(ncol(cols.ts.cpd),0)#@
  nplots=max(nc.gene,nc.cpd)
  pn.suffix=colnames(cols.ts.gene)
  if(length(pn.suffix)<nplots)  pn.suffix=colnames(cols.ts.cpd)
  if(length(pn.suffix)<nplots)  pn.suffix=1:nplot
  if(length(pn.suffix)==1) {
    pn.suffix=out.suffix
  } else pn.suffix=paste(out.suffix, pn.suffix, sep=".")

  if(match.data & nc.gene!=nc.cpd){
  if(nc.gene>nc.cpd) cols.ts.cpd= cols.ts.cpd[, rep(1:nc.cpd, nplots)[1:nplots]]
  if(nc.gene<nc.cpd) cols.ts.gene= cols.ts.gene[, rep(1:nc.gene, nplots)[1:nplots]]
  nc.gene=nc.cpd=nplots
  }

  if(!is.null(cols.ts.gene)){
    nidx.gene=which(nNames %in% rownames(cols.ts.gene))
    cidx.gene=match(nNames[nidx.gene], rownames(cols.ts.gene))
  }
    if(!is.null(cols.ts.cpd)){
    nidx.cpd=which(nNames %in% rownames(cols.ts.cpd))
    cidx.cpd=match(nNames[nidx.cpd], rownames(cols.ts.cpd))
  }

  out.fmt="Working in directory %s"
  wdir=getwd()
  out.msg=sprintf(out.fmt, wdir)
  message(out.msg)
  out.fmt="Writing image file %s"
  
for(np in 1:nplots){
  if(!is.null(cols.ts.gene) & nc.gene>=np){
    fillcol[nidx.gene]=cols.ts.gene[cidx.gene,np] #need to be added in pathview.graph
#  vector indexing gets NA instead of out-of-bound error
#   if(is.null(sub.grp))  get NULL or character(0)
    cn.col=cols.ts.gene[,np][sub2grp[sub.idx,1]]
  } else cn.col[sci.gene]=fillcol[sub2grp[sub.idx,1]]
  rect.col=c(cn.col,fillcol[nSizes==1 & rect.idx])
  rect.col[rect.col==na.col]=NA

  if(!is.null(cols.ts.cpd) & nc.cpd>=np){
    fillcol[nidx.cpd]=cols.ts.cpd[cidx.cpd,np]
  if(sum(cpd.idx)>0){
    ell.col=fillcol[cpd.idx]
    ell.col[ell.col==na.col]=NA
  }
  }
  gfile=paste(pathway.name, pn.suffix[np],"pdf", sep=".")
  out.msg=sprintf(out.fmt, gfile)
  message(out.msg)

  pdf.width=ifelse(same.layer,1.5,1)* pdf.size[1]
    pdf(gfile, width=pdf.width,height=pdf.size[2])
    op <- par(no.readonly = TRUE) # save default, for resetting...
    if(same.layer) nf <- layout(cbind(1,2), c(2,1))#, TRUE)
    rg=renderGraph(gR2.layout)
    gri=graphRenderInfo(rg)
    par(mai=gri$mai, usr=gri$usr)
  
                                        #text(loc, label = labs[nNames], cex = cex)
  if(sum(rect.idx)>0){
    rect(xloc.nd-w.unit,yloc.nd-h.unit/2, xloc.nd+w.unit,yloc.nd+h.unit/2, col=rect.col)
  rect.col[]=NA
  }
  if(sum(cpd.idx)>0) {
    ellipses(xloc.e, yloc.e, w.e, h.e/2, cols=ell.col)
    ell.col[]=NA
  }
  if(sum(cni)>0){
    if(sum(sub.idx)>0) text(xloc, yloc, label = labs[sub2grp[sub.idx,1]], cex = cex)
  }
                                        #text(loc, label = labs[nNames], cex = cex)
  if(sum(!cpd.idx)>0) text(loc[[1]][!cpd.idx], loc[[2]][!cpd.idx], label = labs[nNames[!cpd.idx]], cex = cex)
  if(sum(cpd.idx)>0) {
    if(map.cpdname & !is.null(cols.ts.cpd))  yloc.et=yloc.e+h.e*cpd.lab.offset else  yloc.et=yloc.e
    text(xloc.e, yloc.et, label = labs[nNames[cpd.idx]], cex = cex)
  }

  
  pv.pars=list()
  pv.pars$gsizes=c(gri$bbox[2,1],gri$bbox[2,2])
  if(sum(rect.idx)>0){
    pv.pars$nsizes=c(w.unit,h.unit)
  } else pv.pars$nsizes=c(w.e,h.e)
  pv.pars$op=op
  pv.pars$key.cex=cex*1.5
  pv.pars$key.lwd=1
  pv.pars$sign.cex=1.2*cex

  off.sets=c(x=0,y=0)
  align="n"

 ucol.gene=unique(as.vector(cols.ts.gene))
 na.col.gene=ucol.gene %in% c(na.col, NA)

  if(plot.col.key & !is.null(cols.ts.gene) & !all(na.col.gene))  {
    off.sets=col.key(limit=limit$gene, bins=bins$gene, both.dirs=both.dirs$gene, discrete=discrete$gene, graph.size=pv.pars$gsizes,
      node.size=pv.pars$nsizes, key.pos=key.pos, cex=pv.pars$key.cex, lwd=pv.pars$key.lwd, low=low$gene, mid=mid$gene, high=high$gene, align="n")
    align=key.align
    
  }
  
 ucol.cpd=unique(as.vector(cols.ts.cpd))
 na.col.cpd=ucol.cpd %in% c(na.col, NA)
  if(plot.col.key & !is.null(cols.ts.cpd) & !all(na.col.cpd)) {
    off.sets=col.key(limit=limit$cpd, bins=bins$cpd, both.dirs=both.dirs$cpd, discrete=discrete$cpd, graph.size=pv.pars$gsizes, node.size=pv.pars$nsizes, key.pos=key.pos, off.sets=off.sets, cex=pv.pars$key.cex, lwd=pv.pars$key.lwd, low=low$cpd, mid=mid$cpd, high=high$cpd, align=align)
  }
  
  if(new.signature) pathview.stamp(position=sign.pos, graph.sizes=pv.pars$gsizes, on.kegg=F, cex = pv.pars$sign.cex)
  kegg.legend(edges.only=same.layer)
  par(pv.pars$op)
  dev.off()
}
  
  return(pv.pars)
}

