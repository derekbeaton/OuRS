

####### This one is for continuous data ###########


## should be unpacked. Will need to revisit this.
mcd.outliers.list <- function(
  rob.mds, # dists$rob.md output from cont.mcd
  corrmax.res, # output from cont.corrmax (verify)
  boot.res, # output from cont.boot.sup.u
  rob.mb.lower.cut=0.75, # (lower) cutoff for in-betweeners; if same as upper.cut, in.betweeners will be an empty set
  rob.mb.upper.cut=0.9, # cutoff for in-betweeners (upper) and extreme outliers (lower)
  corrmax.threshold=(1/ncol(corrmax.res))*2, # corrmax threshold, i.e. what minimum percent contrib to include in list
  list.type = "list" # options are "long" or "wide" atm
){

  # add check that list.type has "long" and "wide" options only


  vec.boot.md <- c(boot.res)
  lower.cut <- sort(vec.boot.md)[(length(vec.boot.md) * rob.mb.lower.cut)]
  upper.cut <- sort(vec.boot.md)[(length(vec.boot.md) * rob.mb.upper.cut)]

  outliers <- which(rob.mds >= upper.cut)
  inliers <- which(rob.mds < lower.cut)
  ## if lower.cut and upper.cut are the same then in.betweeniers will be empty.
  in.betweeniers <- which(rob.mds >= lower.cut & rob.mds < upper.cut )

  # proportion of all MDs that this participant possesses (all MDs sum to 1)
  MD_prop <- rob.mds/sum(rob.mds)

  # identifying top contributing variables for each participant
  contrib_prop <- t(apply(corrmax.res,1,function(i){i/sum(i)}))
  idcontrib <- apply(contrib_prop,1,function(i){sort(i[which(i >= corrmax.threshold)],decreasing=T)})

  if(list.type == "long"){
    list1 <- lapply(idcontrib,function(i){data.frame(VAR=names(i),CONTRIB_prop=i)})
    list2 <- do.call("rbind",list1)
    list2$ID <- substr(rownames(list2),1,14)

  }else if(list.type == "wide"){
    list1 <- sapply(idcontrib,function(i){paste(names(i),collapse="; ")})
    list2 <- data.frame(ID=names(list1),VARS=list1)

  }else{
    # Not an option yet!
    stop()
  }

  # merging contrib info with MD_prop
  list3 <- merge(as.data.frame(MD_prop),list2,by.x=0,by.y="ID")
  colnames(list3)[1] <- "ID"

  list.outliers <- list3[which(list3$ID %in% names(outliers)),]
  list.outliers <- list.outliers[order(list.outliers$MD_prop,decreasing = T),] # verify contrib order is still in long

  if(length(in.betweeniers) != 0){
    list.in.betweeniers <- list3[which(list3$ID %in% names(in.betweeniers)),]
    list.in.betweeniers <- list.in.betweeniers[order(list.in.betweeniers$MD_prop,decreasing = T),]

    return(list(extreme.outliers=list.outliers,in.betweeniers=list.in.betweeniers))

  }else{
    return(list(extreme.outliers=list.outliers))
  }


}
