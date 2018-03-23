

cat.mcd.outliers.list <- function(cat.mcd.res, corrmax.res, boot.res, lower.bound=0.75, upper.bound=0.9, corrmax.threshold.type = "level", variable.map, corrmax.threshold,  output.type = "list" ){


  ### trying to capture that variable-wise, not level-wise is wanted.

  if(corrmax.threshold.type=="variable"){
    if(length(variable.map)==ncol(corrmax.res)){ # valid variable.map
      # make some magic happen...
      corrmax.res <- corrmax.res %*% make.data.nominal(as.matrix(variable.map))
    }else{
      warning("Invalid variable map. Will default to level-wise list")
      corrmax.threshold.type <- "level"
    }
  }else{
    if(corrmax.threshold.type!="level"){
      warning("corrmax.threshold.type not recognized. Will default to level-wise list")
      corrmax.threshold.type <- "level"
    }
  }

  if(missing(corrmax.threshold) | is.nan(corrmax.threshold) | is.na(corrmax.threshold) | !is.numeric(corrmax.threshold) | is.infinite(corrmax.threshold) | is.null(corrmax.threshold)){
    corrmax.threshold <- (100/ncol(corrmax.res))*2 ## I could move this below...
  }

  ## should perform a corrmax.threshold check to make sure it's a valid value within the ranges of corrmax.res
    ## else change to default
  #findInterval(1:22,c(5,20),left.open = T)
  #corrmax.threshold quantile(c(corrmax.res),probs = )
  if( findInterval(corrmax.threshold,quantile(c(corrmax.res),probs = c(.1,.9)))!=1  ){
    warning("corrmax.threshold is below 10% or above 90%. Setting to default.")
    corrmax.threshold <- (100/ncol(corrmax.res))*2 ## I could move this below...
  }

  lower.cut <- sort(c(boot.res))[(length(boot.res) * lower.bound)]
  upper.cut <- sort(c(boot.res))[(length(boot.res) * upper.bound)]

  outliers <- which(cat.mcd.res$dists$rob.md >= upper.cut)
  inliers <- which(cat.mcd.res$dists$rob.md < lower.cut)
  betweenliers <- which(cat.mcd.res$dists$rob.md >= lower.cut & cat.mcd.res$dists$rob.md < upper.cut )

  proportional.obs.variance <- (cat.mcd.res$dists$rob.md/sum(cat.mcd.res$dists$rob.md))*100
  top.contributions.per.obs <- apply(corrmax.res,1,function(i){sort(i[which(i >= corrmax.threshold)],decreasing=T)})

  if(output.type == "long"){

    output.temp <- lapply(top.contributions.per.obs,function(i){data.frame(VAR=names(i),CONTRIBUTION.PERCENTAGE=i)})
    output.structure <- do.call("rbind",output.temp)
    output.structure$ID <- rep(names(output.temp),unlist(lapply(output.temp,nrow)))
    rownames(output.structure) <- NULL

  }else if(output.type == "list"){

    output.temp <- sapply(top.contributions.per.obs,function(i){paste(names(i),collapse="; ")})
    output.structure <- data.frame(ID=names(output.temp),VARS=output.temp)

  }else if(output.type == "wide"){

    output.temp <- apply(corrmax.res,2,function(i){ifelse(i >= corrmax.threshold,i,NA)})
    output.structure <- data.frame(ID=rownames(output.temp),output.temp)

  }

  output.structure <- merge(as.data.frame(proportional.obs.variance),output.structure,by.x=0,by.y="ID")
  colnames(output.structure)[1] <- "ID"

  list.outliers <- output.structure[which(output.structure$ID %in% names(outliers)),]
  list.outliers <- list.outliers[order(list.outliers$proportional.obs.variance,decreasing = T),] # verify contrib order is still in long

  list.betweenliers <- output.structure[which(output.structure$ID %in% names(betweenliers)),]
  list.betweenliers <- list.betweenliers[order(list.betweenliers$proportional.obs.variance,decreasing = T),]

  if(output.type=="wide"){
    if(any(apply(list.outliers,2,function(i){all(is.na(i))}))){
      list.outliers <- list.outliers[,-which(apply(list.outliers,2,function(i){all(is.na(i))}))]
    }
    if(any(apply(list.betweenliers,2,function(i){all(is.na(i))}))){
      list.betweenliers <- list.betweenliers[,-which(apply(list.betweenliers,2,function(i){all(is.na(i))}))]
    }
  }

  ##replace NAs with blanks here.
  return(list(extreme.outliers=list.outliers,betweenliers=list.betweenliers))

}
