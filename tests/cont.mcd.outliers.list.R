
# turn these into help file with roxygen2 or devtools.
#
# cont.mcd.res$dists$rob.md, # now the full results; can switch distances by the dist.type function.
# corrmax.res, # output from cont.corrmax (verify)
# boot.res, # output from cont.boot.sup.u
# lower.bound=0.75, # (lower) cutoff for in-betweeners; if same as upper.cut, in.betweeners will be an empty set
# upper.bound=0.9, # cutoff for in-betweeners (upper) and extreme outliers (lower)
# corrmax.threshold=(1/ncol(corrmax.res))*2, # corrmax threshold, i.e. what minimum percent contrib to include in list
# list.type = "list", # options are "long", "list", or "wide"
# dist.type = "mahal"


## someday soon we'll need to update the code so that output from any of the functions can pass through to the reporting function
  ## however for now we'll keep everything in "parallel" streams with the cat. and cont. prefixes.
cont.mcd.outliers.list <- function(cont.mcd.res, corrmax.res, boot.res, lower.bound=0.75, upper.bound=0.9, corrmax.threshold=((100/ncol(corrmax.res))*2),  output.type = "list" ){


    # can bring this back later.
  # if( !(dist.type %in% c("md","cd","od")) ){
  #   warning("Distance type not recognized. Setting to Mahalanobis distance (`md``).")
  #   dist.type <- "md"
  # }
  if(!(output.type %in% c("long","list","wide"))){
    warning("Output type not recognized. Setting to wide (`wide`).");
    list.type <- "wide"
  }

  ## should perform a corrmax.threshold check to make sure it's a valid value within the ranges of corrmax.res
    ## else change to default
  if( findInterval(corrmax.threshold,quantile(c(corrmax.res),probs = c(.1,.9)))!=1  ){
    warning("corrmax.threshold is below 10% or above 90%. Setting to default.")
    corrmax.threshold <- (100/ncol(corrmax.res))*2 ## I could move this below...
  }

  lower.cut <- sort(c(boot.res))[(length(boot.res) * lower.bound)]
  upper.cut <- sort(c(boot.res))[(length(boot.res) * upper.bound)]

  outliers <- which(cont.mcd.res$dists$rob.md >= upper.cut)
  inliers <- which(cont.mcd.res$dists$rob.md < lower.cut)
  betweenliers <- which(cont.mcd.res$dists$rob.md >= lower.cut & cont.mcd.res$dists$rob.md < upper.cut )


  proportional.obs.variance <- (cont.mcd.res$dists$rob.md/sum(cont.mcd.res$dists$rob.md))*100
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

  if(list.type=="wide"){
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
