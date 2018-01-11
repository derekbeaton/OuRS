
##### Making the lists of outliers #####

outliers.list <- function(
  dists,
  corrmax.res,
  boot.res,
  in.betweeners = T,
  outfile,
  corrmax.threshold,
  longver = T,
  widever = T,
  corrver = T
){

  # Filenames
  mydate <- toupper(format(Sys.Date(),"%Y%b%d"))
  longfile <- sprintf("%s//Outliers_LongList_%s.csv",outfile,mydate)
  widefile <- sprintf("%s//Outliers_WideList_%s.csv",outfile,mydate)
  corrfile <- sprintf("%s//Outliers_CorrmaxList_%s.csv",outfile,mydate)

  # Identifying Outliers

  vec.boot.md <- c(boot.res)
  lower.cut <- sort(vec.boot.md)[(length(vec.boot.md) * .75)]
  upper.cut <- sort(vec.boot.md)[(length(vec.boot.md) * .9)]

  # extreme.outliers <- names(which(dists >= upper.cut))
  outliers <- names(which(dists >= upper.cut))

  if(in.betweeners){
    # between.outliers <- names(which(dists < upper.cut & dists > lower.cut))
    outliers <- names(which(dists >= lower.cut))
  }

      ## clindat3 does not exist
  corrmax.sum <- matrix(NA, nrow=nrow(corrmax.res),ncol=ncol(clindat3),dimnames=list(rownames(corrmax.res),colnames(clindat3))) ### CLINDAT3 NEEDS TO BE CALLED INTO THE FUNCTION

  for(i in colnames(clindat3)){
    # cat(sprintf("%s: %s\n",i,paste(grep(paste0("^",i),colnames(corrmax.res),val=T),collapse=", ")))
    corrmax.sum[,i] <- rowSums(corrmax.res[,grep(paste0("^",i),colnames(corrmax.res))])
      ## this may not work if similar enough stems exist.
  }

      ## this should be !=
  if(all(rowSums(corrmax.sum) == 100)){cat("PROBLEM WITH SUMMING!\n");browser()}

  MD_prop <- data.frame(MD_prop=dists/sum(dists))

  allvar <- merge(MD_prop,corrmax.sum,by=0)
  rownames(allvar) <- allvar$Row.names
  allvar.outliers <-allvar[outliers,-1]
  # allvar.outliers <- allvar.outliers[order(allvar.outliers$MD_prop,decreasing = T),]

  if(longver|widever){
    idcontrib <- apply(allvar.outliers[,-1],1,function(i){sort(i[which(i >= corrmax.threshold)],decreasing=T)})

    if(longver){
      longlist <- lapply(idcontrib,function(i){data.frame(VAR=names(i),CONTRIB=i)})
      longlist1 <- do.call("rbind",longlist)
      longlist1$ID <- substr(rownames(longlist1),1,14)
      longlist2 <- merge(MD_prop,longlist1,by.x=0,by.y="ID")
      colnames(longlist2)[1] <- "ID"
      longlist.final <- longlist2[order(longlist2$MD_prop,longlist2$CONTRIB,decreasing = T),]
      write.csv(longlist.final,file=longfile,na="",row.names=F)
    } # end of longver

    if(widever){
      wide.ncols <- max(sapply(idcontrib,function(i){length(i)}))

      widelist <-matrix(NA,nrow=length(idcontrib)*2,ncol=wide.ncols+1)
      colnames(widelist) <- c("ID",paste0("VAR",1:wide.ncols))
      widelist[,"ID"] <- rep(names(idcontrib),each=2)

      for(j in names(idcontrib)){
        tmp <- idcontrib[[j]]
        widelist[which(widelist[,"ID"] == j),][1,2:(length(tmp)+1)] <- names(tmp)
        widelist[which(widelist[,"ID"] == j),][2,2:(length(tmp)+1)] <- tmp
      }

      # merge back MD and order
      widelist1 <- merge(MD_prop,widelist,by.x=0,by.y="ID")
      colnames(widelist1)[1] <- "ID"
      widelist.final <- widelist1[order(widelist1$MD_prop,decreasing = T),]
      write.csv(widelist.final,file=widefile,na="",row.names=F)
    } # end of widever
  } # end of longver | widever

  if(corrver){
    select.contrib <- apply(allvar.outliers[,-1],2,function(i){any(i >= corrmax.threshold)})
    corrlist.final <- data.frame(ID=rownames(allvar.outliers),allvar.outliers[,c("MD_prop",sort(names(select.contrib)))])
    write.csv(corrlist.final,file=corrfile,na="",row.names=F)
  } # end of corrver

} # end of function
