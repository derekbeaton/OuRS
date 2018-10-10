  ### like with cat.mcd, assumed only categorical but is actually not...
    ### we need generalizations of these.
#two.fold.repeated.ca <- function(DATA, make.data.disjunctive=F, iters=500,sh1.size=.5,k=0){
two.fold.repeated.ca <- function(DATA, iters=500,sh1.size=.5,k=0){

  if(sh1.size > .9 | sh1.size < .1){
    warning("sh1.size is greater than or equal to 90% or less than or equal to 10%. Setting sh1.size to .5")
    sh1.size <- .5
  }

    ## A weird choice I made but whatever!
  if(nrow(DATA) < 20){
    warning("DATA has less than 20 rows. sh1.size will automatically be set to 50%")
    sh1.size <- .5
  }

  ## do the intitial CA here anyways. We'll be able to use it for a lot of these things.
  ca.res <- ca(DATA,k=k,compact=F)
  max.rank <- length(ca.res$d.orig)

  sh1.orders <- matrix(NA,iters,ceiling(nrow(DATA)*sh1.size))
  sh2.orders <- matrix(NA,iters,nrow(DATA)-ncol(sh1.orders))
  sh.dets <- matrix(NA,iters,2)
  #fi.score.cors <- u.score.cors <-
  u.score.cors <- array(NA,dim=c(max.rank,max.rank,iters)) -> v.loadings.cors


    ### the pred.fi may not be any different.
  #pred.fi.array <- pred.p.array <- pred.u.array <- array(NA,dim=c(nrow(DATA),min(dim(DATA)),iters))
  pred.fi.array <- pred.u.array <- array(NA,dim=c(nrow(DATA),min(dim(DATA)),iters))
  rownames(pred.fi.array) <- rownames(pred.u.array) <- rownames(DATA)

  preproc.data <- ca.preproc(DATA)
  profiles <- sweep(preproc.data$Ox,1,preproc.data$m,"/")
  for(i in 1:iters){

    ## the sort is simply for aesthetics...
    sh1 <- sort(sample(nrow(DATA),ceiling(nrow(DATA)*sh1.size)))
    sh2 <- sort(setdiff(1:nrow(DATA),sh1))

    sh1.res <- gsvd(preproc.data$weightedZx[sh1, ], k=min(k,max.rank))
    sh2.res <- gsvd(preproc.data$weightedZx[sh2, ], k=min(k,max.rank))

    sh.dets[i,] <- c(geometric.mean(sh1.res$d^2),geometric.mean(sh2.res$d^2))
    sh1.orders[i,] <- sh1
    sh2.orders[i,] <- sh2



    sh1.projs <- cat.sup.fi.u(profiles[sh1, ], preproc.data$m[sh1], preproc.data$w, sh2.res$v, sh2.res$d)
      pred.fi.array[sh1,,] <- sh1.projs$sup.fi
      pred.u.array[sh1,,] <- sh1.projs$sup.u
    sh2.projs <- cat.sup.fi.u(profiles[sh2, ], preproc.data$m[sh2], preproc.data$w, sh1.res$v, sh1.res$d)
      pred.fi.array[sh2,,] <- sh2.projs$sup.fi
      pred.u.array[sh2,,] <- sh2.projs$sup.u


    v.loadings.cors[1:min(c(length(sh1.res$d),length(sh2.res$d))),1:min(c(length(sh1.res$d),length(sh2.res$d))),i] <- cor(
      sh1.res$v[,1:min(c(length(sh1.res$d),length(sh2.res$d)))],sh2.res$v[,1:min(c(length(sh1.res$d),length(sh2.res$d)))])

    u.score.cors[1:min(c(length(sh1.res$d),length(sh2.res$d))),1:min(c(length(sh1.res$d),length(sh2.res$d))),i] <- cor(
      rbind(sh1.res$u[,1:min(c(length(sh1.res$d),length(sh2.res$d)))],sh2.res$u[,1:min(c(length(sh1.res$d),length(sh2.res$d)))]),
      pred.u.array[c(sh1,sh2),1:min(c(length(sh1.res$d),length(sh2.res$d))),i])

  }

  return( list(pred.fi.array=pred.fi.array,pred.u.array=pred.u.array,sh1.orders=sh1.orders,sh2.orders=sh2.orders,sh.dets=sh.dets,v.loadings.cors=v.loadings.cors,u.score.cors=u.score.cors) )

}
