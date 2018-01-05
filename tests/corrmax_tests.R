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




db.sampvar <- var(datamat) ## "reference" or "population" data
db.sqrt.diag.sampvar <- sqrt(diag(db.sampvar))
db.invvar <- db.sampvar %^% (-1)

db.DSD <- sweep(sweep(db.sampvar,1,db.sqrt.diag.sampvar,"/"),2,db.sqrt.diag.sampvar,"/")
db.invDSDhalf <- db.DSD %^% (-1/2)
db.diff.dat <- expo.scale(IndivObs,center=colMeans(rawdata),scale=F) #signs are switched; I can just *-1
db.Tsquares <- diag(db.diff.dat %*% db.invvar %*% t(db.diff.dat))
db.Ws <- t(sweep(db.invDSDhalf,2,db.sqrt.diag.sampvar,"/") %*% t(db.diff.dat))
db.GK <- db.Ws * db.Ws
db.percentages <- sweep(db.GK,1,rowSums(db.GK),"/")*100

  ## don't know what this is really for...
db.corrs <- diag(db.invDSDhalf %^% (-1))

#our.corrmax <- corrmax(IndivObs,rawdata)
