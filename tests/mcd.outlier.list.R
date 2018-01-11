
  ## does not work. Will need to revisit this.
mcd.outliers.list <- function(
  rob.mds,
  corrmax.res,
  boot.res,
  rob.mb.lower.cut=.9,
  rob.mb.upper.cut=.9,
  corrmax.threshold=1/ncol(corrmax.res),
  longver = T,
  widever = T,
  corrver = T
){

  vec.boot.md <- c(boot.res)
  lower.cut <- sort(vec.boot.md)[(length(vec.boot.md) * rob.mb.lower.cut)]
  upper.cut <- sort(vec.boot.md)[(length(vec.boot.md) * rob.mb.upper.cut)]

  outliers <- which(rob.mds >= upper.cut)
  inliers <- which(rob.mds < lower.cut)
      ## if lower.cut and upper.cut are the same then in.betweeniers will be empty.
  in.betweeniers <- which(rob.mds >= lower.cut & rob.mds < upper.cut )

  md.prop <- rob.mds/sum(rob.mds)

    allvar <- merge(md.prop,corrmax.res,by=0)
    rownames(allvar) <- rownames(corrmax)
    allvar.outliers <- allvar[outliers,-1]


    idcontrib <- apply(allvar.outliers[,-1],1,function(i){sort(i[which(i >= corrmax.threshold)],decreasing=T)})


    longlist <- lapply(idcontrib,function(i){data.frame(VAR=names(i),CONTRIB=i)})
    longlist1 <- do.call("rbind",longlist)
    longlist1$ID <- substr(rownames(longlist1),1,14)
    longlist2 <- merge(MD_prop,longlist1,by.x=0,by.y="ID")
    colnames(longlist2)[1] <- "ID"
    longlist.final <- longlist2[order(longlist2$MD_prop,longlist2$CONTRIB,decreasing = T),]
    wide.ncols <- max(sapply(idcontrib,function(i){length(i)}))


    widelist <-matrix(NA,nrow=length(idcontrib)*2,ncol=wide.ncols+1)
    colnames(widelist) <- c("ID",paste0("VAR",1:wide.ncols))
    widelist[,"ID"] <- rep(names(idcontrib),each=2)

      ## is loop necessary?
    for(j in names(idcontrib)){
      tmp <- idcontrib[[j]]
      widelist[which(widelist[,"ID"] == j),][1,2:(length(tmp)+1)] <- names(tmp)
      widelist[which(widelist[,"ID"] == j),][2,2:(length(tmp)+1)] <- tmp
    }

    # merge back MD and order
    widelist1 <- merge(MD_prop,widelist,by.x=0,by.y="ID")
    colnames(widelist1)[1] <- "ID"
    widelist.final <- widelist1[order(widelist1$MD_prop,decreasing = T),]




}
