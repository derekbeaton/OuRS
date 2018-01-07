## via: http://users.mct.open.ac.uk/paul.garthwaite/Mahalanobis.R

  ##
rawdata <- read.table("http://users.mct.open.ac.uk/paul.garthwaite/SwissNotes.txt",header=TRUE)
IndivObs <- read.table("http://users.mct.open.ac.uk/paul.garthwaite/SingleItems.txt",header=TRUE)

attach(rawdata)
attach(IndivObs)
namesA <- names(rawdata)
datamat <- cbind(rawdata)
individuals <- cbind(IndivObs)

nrA=nrow(datamat)
nrB=nrow(individuals)
ncA=ncol(datamat)
ncB=ncol(individuals)

#  obtain the sample mean and sample variance
meanA<-colMeans(datamat)
SampVar <- var(datamat)
InvVarA <- solve(SampVar)

# Start that part of the Garthwaite-Koch partition that is the same
# for all items
InvSD <- rep(0,ncA)
for(i in 1:ncA){InvSD[i] <- 1/sqrt(SampVar[i,i])}
Dmat<-diag(InvSD)
DSD <- Dmat %*% SampVar %*% Dmat



# We next calculate the inverse of the symmetric square-root of DSD
eig <- eigen(DSD)
InvRootEig <- rep(0,ncA)
for(i in 1:ncA){InvRootEig[i] <- 1/sqrt(eig$val[i])}
InvDSDhalf <-  eig$vec %*% diag(InvRootEig) %*% t(eig$vec)





ItemVal<-as.numeric(colMeans(individuals))
# Take each individual item in turn
for(i2 in 1:nrB)
{
  for(j in 1:ncA){ItemVal[j] <- individuals[i2,j]}
  diff <- meanA -ItemVal
  # diff contains the difference between the sample mean and
  # individual item of current interest

  # Next calculate Mahalanobis distance
  Qform <- t(diff) %*% InvVarA %*% diff
  Tsquare <- Qform[1,1]
  W <- InvDSDhalf %*% Dmat %*% diff
  GKcontrib <- diag(W %*% t(W))
  Percentage <- round(100 * GKcontrib /sum(diag(GKcontrib)), digits=3)
  Correlation <- round(diag(solve(InvDSDhalf)),digits=4)
  Contribution <- as.numeric(round(GKcontrib, digits=3))
  GKframe <- data.frame(cbind(Contribution, Percentage, Correlation))
}


## KS CODE


#mcd1@cov

# Based on 2016 paper: Contributions to Quadratic Form
smpsd <- sqrt(diag(cov(datamat)))
d <- diag(1/smpsd)      ## DB NOTE: a fraction can make a diag() panic.
dsd <- d %*% cov(datamat) %*% d   ## DB Q: Is this converting cov to cor?

trans1 <- sweep(IndivObs,2, colMeans(datamat)) ## DB Q: Center data by robust mean?
#corrmax <- sqrtMat(dsd,-1/2) %*% d ## DB Q: inverse sqrt of cor times data?
corrmax <- dsd %^% (-1/2) %*% d
colnames(corrmax) <- rownames(corrmax) <- colnames(IndivObs)
contrib <- apply(trans1,1,function(i){corrmax %*% as.matrix(i)}) # DB Q: data times correlation?
w2 <- t(contrib**2) ## DB Q: squared "contributions"?


## MY CODE

db.sampvar <- var(datamat) ## "reference" or "population" data
db.sqrt.diag.sampvar <- sqrt(diag(db.sampvar))
db.invvar <- db.sampvar %^% (-1)
db.DSD <- sweep(sweep(db.sampvar,1,db.sqrt.diag.sampvar,"/"),2,db.sqrt.diag.sampvar,"/")
  ## so there is a way!
  #db.DSD2 <- tcrossprod(sweep(sweep(db.svd.res$v,2,db.svd.res$d,"*"),1,db.sqrt.diag.sampvar,"/"))  / (nrow(datamat)-1)
  #db.DSD2 / db.DSD
db.invDSDhalf <- db.DSD %^% (-1/2)
  #db.invDShalf2 <- (tcrossprod(sweep(sweep(db.svd.res$v,2,db.svd.res$d,"*"),1,db.sqrt.diag.sampvar,"/") %^% (-1/2)))   * sqrt((nrow(datamat)-1))
  #db.invDShalf2 / db.invDSDhalf
db.diff.dat <- expo.scale(IndivObs,center=colMeans(rawdata),scale=F) #signs are switched; I can just *-1
#db.Tsquares <- diag(db.diff.dat %*% db.invvar %*% t(db.diff.dat))
db.Ws <- t(sweep(db.invDSDhalf,2,db.sqrt.diag.sampvar,"/") %*% t(db.diff.dat))
  #sweep(db.invDSDhalf,2,db.sqrt.diag.sampvar,"/")
db.GK <- db.Ws * db.Ws
db.percentages <- sweep(db.GK,1,rowSums(db.GK),"/")*100

  ## don't know what this is really for...
db.corrs <- diag(db.invDSDhalf %^% (-1))


  ## I can get this even shorter, too, somehow...
db2.svd.res <- tolerance.svd(expo.scale(datamat,center=T,scale=F))
  #db2.sqrt.diag.sampvar <- sqrt(diag(tcrossprod(db2.svd.res$v %*% diag(db2.svd.res$d))) / (nrow(datamat)-1))
  #db2.invDShalf2 <- (tcrossprod(sweep(sweep(db2.svd.res$v,2,db2.svd.res$d,"*"),1,db2.sqrt.diag.sampvar,"/") %^% (-1/2)))   * sqrt((nrow(datamat)-1))
db2.sqrt.diag.sampvar <- sqrt(diag(tcrossprod(db2.svd.res$v %*% diag(db2.svd.res$d))))


sweep(db2.svd.res$v,2,db2.svd.res$d,"*") / (db2.svd.res$v %*% diag(db2.svd.res$d))

#sweep(db2.svd.res$v,2,db2.svd.res$d,"*") %*% t(sweep(db2.svd.res$v,2,db2.svd.res$d,"*"))

#sweep(db2.svd.res$v,2,db2.svd.res$d,"*") * sweep(db2.svd.res$v,2,db2.svd.res$d,"*")


db2.invDShalf <- (tcrossprod(sweep(sweep(db2.svd.res$v,2,db2.svd.res$d,"*"),1,db2.sqrt.diag.sampvar,"/") %^% (-1/2)))
db2.diff.dat <- expo.scale(IndivObs,center=colMeans(rawdata),scale=F) #signs are switched; I can just *-1
db2.Ws <- t(sweep(db2.invDShalf,2,db2.sqrt.diag.sampvar,"/") %*% t(db2.diff.dat))

db2.Ws / (db2.diff.dat %*% sweep(db2.invDShalf,1,db2.sqrt.diag.sampvar,"/"))

db2.contribs <- db2.Ws^2
#db2.percentages <- sweep(db2.contribs,1,rowSums(db2.contribs),"/")*100
db2.percentages <- sweep(db2.Ws*db2.Ws,1,rowSums(db2.contribs),"/")*100

  ## ok so we should return percentages and use those as they have meaning and they are just the contribs scaled.

  ## NOTE: I need to functionalize (for testing) the KS corr-max code
  ## and allow for input of arbitrary cov, center, data; I can get the vectors and SVs from robustbase's @cov and use that in my cont.corrmax call.
## ok what the hell
diff.svd <- tolerance.svd(db2.diff.dat)
w.svd <- tolerance.svd(db2.Ws)
indobs.svd <- tolerance.svd(expo.scale(IndivObs,scale=F))

  ## these are the same but the vectors are not...
rowSums(diff.svd$u^2)
rowSums(w.svd$u^2)
  ## so corrmax solution lives within the Mahal space of the cov by the target?
  ## OK so because this is true that means there *must* be a faster/easier solution through the SVD.


## what about this:
contribs.ca <- ca(db2.contribs)
percs.ca <- ca(db2.percentages)


## corrmax needs:
  # loadings, svs, sample size (can be compute dfrom data), data (duh)
    ## that should be it...

#our.corrmax <- corrmax(IndivObs,rawdata)


corrmax.temp <- function(target.data,rob.center=T,rob.scale=F,loadings,singular.values){

    ## if we want to do a correction for sample size in this case, all we need to do is scale up the sv/eigen values.

  diag.sampcov <- sqrt(diag(tcrossprod( sweep(loadings,2,singular.values,"*") )))
  inv.DSD.half <- (tcrossprod(sweep( sweep(loadings,2,singular.values,"*"),1,diag.sampcov,"/") %^% (-1/2)))
  target.data <- expo.scale(target.data,center=rob.center,scale=rob.scale) #signs are switched; I can just *-1
  W <- t(sweep(inv.DSD.half,2,diag.sampcov,"/") %*% t(target.data))
  W <- W*W

  return(sweep(W,1,rowSums(W),"/")*100)

}

datamat.cs <- expo.scale(datamat,T,F)
svd.res <- tolerance.svd(datamat.cs)

#cm.res <- corrmax.temp(IndivObs,rob.center=attributes(datamat.cs)$`scaled:center`,rob.scale = attributes(datamat.cs)$`scaled:scale`,svd.res$v,svd.res$d)
library(cellWise)
data("philips")
ours.res <- cont.mcd(philips,num.subsets = 500,alpha = .75)
cm.res <- cont.corrmax(philips,ours.res$best.center,ours.res$best.scale,ours.res$best.loadings,ours.res$best.svs)

philips2 <- expo.scale(philips,ours.res$best.center,F)
phil.svd <- tolerance.svd(philips2)
  ## somehow this is just the plain MD.


### categorical corrmax?
data("SNPS")
X <- make.data.nominal(SNPS)
preproc.res <- ca.preproc(X)


#tsvd.res <- tolerance.svd(crossprod(preproc.res$weightedZx))
wMah <- preproc.res$weightedZx %*% (crossprod(preproc.res$weightedZx) %^% (-1)) %*% t(preproc.res$weightedZx)
Mah <- preproc.res$Zx %*% (crossprod(preproc.res$Zx) %^% (-1)) %*% t(preproc.res$Zx)
diag(wMah) / diag(Mah)


  ## ok so the condtion of corrmax above is that the target data Mahals should end up the same as the W Mahals.
snps.res <- cat.mcd(SNPS,make.data.disjunctive = T,num.subsets = 500,alpha = .5)
  ## use this as reference

snps.svd <- tolerance.svd(preproc.res$weightedZx)

diag.sampcov <- sqrt(diag(tcrossprod( sweep(snps.res$best.loadings,2,snps.res$best.svs,"*") )))
inv.DSD.half <- (tcrossprod(sweep( sweep(snps.res$best.loadings,2,snps.res$best.svs,"*"),1,diag.sampcov,"/") %^% (-1/2)))

W <- t(sweep(inv.DSD.half,2,diag.sampcov,"/") %*% t(preproc.res$Zx))
W2 <- W*W

  ## well hot damn.
percs <- sweep(W2,1,rowSums(W2),"/")*100



snps.w.svd <- tolerance.svd(W)

  ## well I think we've solved it.
rowSums(snps.w.svd$u^2) / rowSums(snps.svd$u^2)
snps.w.svd$u / snps.svd$u


heatmap(percs * make.data.nominal(SNPS),Rowv=NA,Colv=NA)
