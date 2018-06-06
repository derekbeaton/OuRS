

require(ggplot2)
# require(gridExtra)


mcd.outliers.plot <- function(
  target.data,
  dists, # dists$rob.md output from cont.mcd
  corrmax.res, # output from cont.corrmax (verify)
  boot.res, # output from cont.boot.sup.u
  rob.mb.lower.cut=0.75, # (lower) cutoff for in-betweeners; if same as upper.cut, in.betweeners will be an empty set
  rob.mb.upper.cut=0.9, # cutoff for in-betweeners (upper) and extreme outliers (lower)
  corrmax.threshold=(1/ncol(corrmax.res))*2, # corrmax threshold, i.e. what minimum percent contrib to include in list
  nvars=3
){

  # copied from mcd.outliers.list()

  rob.mds <- dists$rob.md

  vec.boot.md <- c(boot.res)
  lower.cut <- sort(vec.boot.md)[(length(vec.boot.md) * rob.mb.lower.cut)]
  upper.cut <- sort(vec.boot.md)[(length(vec.boot.md) * rob.mb.upper.cut)]

  # proportion of all MDs that this participant possesses (all MDs sum to 1)
  MD_prop <- rob.mds/sum(rob.mds)

  outliers <- names(which(rob.mds >= upper.cut))
  ## if lower.cut and upper.cut are the same then in.betweeniers will be empty.
  in.betweeniers <- names(which(rob.mds >= lower.cut & rob.mds < upper.cut ))



  # identifying top contributing variables for each participant
  contrib_prop <- t(apply(corrmax.res,1,function(i){i/sum(i)}))
  idcontrib <- apply(contrib_prop,1,function(i){sort(i[which(i >= corrmax.threshold)],decreasing=T)})

  topcontribs <- t(sapply(idcontrib,function(i){names(i)[1:3]})) # limited to 3 here...
  colnames(topcontribs) <- paste0("VAR",1:3)

  mapdat <- merge(as.data.frame(MD_prop),data.frame(topcontribs),by=0)
  colnames(mapdat)[1] <- "ID"
  mapdat$FLAG <- ifelse(mapdat$ID %in% outliers,"O",ifelse(mapdat$ID %in% in.betweeniers,"B","I"))

  # For global plots
  distsdat <- merge(merge(as.data.frame(dists$md),as.data.frame(dists$rob.md),by=0),as.data.frame(dists$rob.chid),by.x=1,by.y=0)
  colnames(distsdat) <- c("ID","MD","ROB.MD","ROB.CHI2")

  mapdat <- merge(mapdat,distsdat)


  # Plots
  listplots <- list()

# Individual outlier plots


  for(i in outliers){ # or outliers and inbetweeniers

    tmpvars <- c(as.matrix(mapdat[which(mapdat$ID == i),grep("VAR",colnames(mapdat))]))
    tmpvars <- tmpvars[which(!is.na(tmpvars))]
    tmpdat <- data.frame(target.data[,tmpvars])
    colnames(tmpdat) <- paste0("VAR",1:length(tmpvars))
    
    tmpdat <- merge(tmpdat,mapdat[,c("ID","FLAG")],by.x=0,by.y=1)
    tmpdat$PART_d <- ifelse(tmpdat$Row.names == i,"FLAG","")
    tmpdat$PART_c <- ifelse(tmpdat$Row.names == i,1,0)
    
    plotvars <- ifelse(nvars > length(tmpvars), length(tmpvars),nvars) # the number of variables desired/able to plot

    if(plotvars == 1){
      
      listplots[[i]] <- ggplot(tmpdat,aes(x=tmpvars,y=VAR1,col=FLAG,shape=PART_d,size=PART_c)) + geom_jitter() +
        labs(title=i,x="",y=tmpvars[1],color=element_blank()) + scale_shape_discrete(guide=F) + scale_size(guide=F,range=c(3,5)) +
        scale_color_manual(breaks=c("I","B","O"),labels=c("Inlier","In-betweenlier","Outlier"),values=c("blue","grey","red")) +
        theme_classic()
      
    }else if(plotvars == 2){

      listplots[[i]] <- ggplot(tmpdat,aes(x=VAR1,y=VAR2,col=FLAG,shape=PART_d,size=PART_c)) + geom_point() +
        labs(title=i,x=tmpvars[1],y=tmpvars[2],color=element_blank()) + scale_shape_discrete(guide=F) + scale_size(guide=F,range=c(3,5)) +
        scale_color_manual(breaks=c("I","B","O"),labels=c("Inlier","In-betweenlier","Outlier"),values=c("blue","grey","red")) +
        theme_classic()
      
    }else if(plotvars == 3){

      listplots[[i]] <- ggplot(tmpdat,aes(x=VAR1,y=VAR2,col=VAR3,shape=PART_d,size=PART_c)) + geom_point() +
        labs(title=i,x=tmpvars[1],y=tmpvars[2],color=tmpvars[3]) + scale_shape_discrete(guide=F) + scale_size(guide=F,range=c(3,5)) +
        theme_classic()
      
    }else{
      cat("Problem: Either not set up or unexpected result"); stop()
    }
  }

return(listplots)

}
