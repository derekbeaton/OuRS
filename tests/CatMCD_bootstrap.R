data(SNPS)

snps.res <- cat.mcd(SNPS,make.data.disjunctive = T,num.subsets = 500,alpha = .5)
ca.preproc.res <- ca.preproc(make.data.nominal(SNPS))
profiles <- sweep(ca.preproc.res$Ox,1,ca.preproc.res$m,"/")



## BOOT!
boot.res <- cat.boot.sup.u(make.data.nominal(SNPS),snps.res$cov$loadings,snps.res$cov$singular.values,iters = 100)

## compare these two.
hist(c(boot.res))
hist(snps.res$dists$rob.md)


that.boot <- sample(nrow(SNPS),replace=T)

ca.preproc.data <- ca.preproc(make.data.nominal(SNPS)[that.boot,])
loadings <- snps.res$cov$loadings
singular.values <- snps.res$cov$singular.values

these.vecs <- sweep(loadings,1,sqrt(ca.preproc.data$w)/ca.preproc.data$w,"*")
these.vecs[is.nan(these.vecs)] <- 0

this.dumb.thing <- sweep(
  sweep(

    sweep(ca.preproc.data$Ox,1,ca.preproc.data$m,"/") %*% these.vecs

    ,2,singular.values,"/")
  ,1,sqrt(ca.preproc.data$m),"*")


#this.dumb.thing

rowSums(this.dumb.thing^2)

