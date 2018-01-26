

require(ggplot)


mcd.outliers.plot <- function(
  target.data,
  rob.mds, # dists$rob.md output from cont.mcd
  corrmax.res, # output from cont.corrmax (verify)
  boot.res, # output from cont.boot.sup.u
  rob.mb.lower.cut=0.75, # (lower) cutoff for in-betweeners; if same as upper.cut, in.betweeners will be an empty set
  rob.mb.upper.cut=0.9, # cutoff for in-betweeners (upper) and extreme outliers (lower)
  corrmax.threshold=(1/ncol(corrmax.res))*2 # corrmax threshold, i.e. what minimum percent contrib to include in list
){

  # copied from mcd.outliers.list()

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

  browser()

  topcontribs <- t(sapply(idcontrib,function(i){names(i)[1:2]}))
  colnames(topcontribs) <- c("VAR1","VAR2")

  mapdat <- merge(as.data.frame(MD_prop),data.frame(topcontribs),by=0)
  colnames(mapdat)[1] <- "ID"
  mapdat$FLAG <- ifelse(mapdat$ID %in% outliers,"O",ifelse(mapdat$ID %in% in.betweeniers,"B","I"))
  mapdat <- mapdat[order(mapdat$MD_prop,decreasing = T),]


  # eventually make a plot matrix with the "main" plot in the centre, but now just make a bunch of plots to start

  for(i in 1:length(outliers)){ # or outliers and inbetweeniers
    tmpvars <- c(as.matrix(mapdat[i,c("VAR1","VAR2")]))
    tmpdat <- target.data[,]
  }








}
