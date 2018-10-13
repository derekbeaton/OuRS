rm(list=ls())
library(GSVD)

### some tests.

data("dfw.beer.survey")
nom.dat <- make.data.nominal(DFW.beer.survey)


tfca.res <- two.fold.repeated.ca(nom.dat,iters = 100)

ca.res <- ca(nom.dat)
max.rank <- length(ca.res$d.orig)
preproc.data <- ca.preproc(nom.dat)
profiles <- sweep(preproc.data$Ox,1,preproc.data$m,"/")


sh1 <- tfca.res$sh1.orders[100,]
sh2 <- sort(setdiff(1:nrow(nom.dat),sh1))

sh1.res <- tolerance.svd(preproc.data$weightedZx[sh1,])
sh2.res <- tolerance.svd(preproc.data$weightedZx[sh2,])

### (1) do I get the same results through the GSVD?
  ### note: this approach could be moved over to the cat.mcd (and its future versions) so that we can actually recover P & Q.
sh1.res_v2 <- gsvd(preproc.data$Zx[sh1,], (1/preproc.data$m)[sh1], 1/preproc.data$w)
sh2.res_v2 <- gsvd(preproc.data$Zx[sh2,], (1/preproc.data$m)[sh2], 1/preproc.data$w)



sh1.res_v2$u / sh1.res$u
sh1.res_v2$v / sh1.res$v

sh2.res_v2$u / sh2.res$u
sh2.res_v2$v / sh2.res$v

### (2) can I recover these same values through projections?
  fi1 <- profiles[sh1,] %*% sweep(sh1.res$v, 1, sqrt(preproc.data$w)/preproc.data$w,"*")
  p1 <- sweep(sweep(fi1, 2, sh1.res$d,"/"), 1, preproc.data$m[sh1],"*")
  u1 <- sweep(p1, 1, sqrt(preproc.data$m[sh1]),"/")

  fi2 <- profiles[sh2,] %*% sweep(sh2.res$v, 1, sqrt(preproc.data$w)/preproc.data$w,"*")
  p2 <- sweep(sweep(fi2, 2, sh2.res$d,"/"), 1, preproc.data$m[sh2],"*")
  u2 <- sweep(p2, 1, sqrt(preproc.data$m[sh2]),"/")


fi1 / sh1.res_v2$fi
p1 / sh1.res_v2$p
u1 / sh1.res_v2$u
#u1 / sh1.res$u

fi2 / sh2.res_v2$fi
p2 / sh2.res_v2$p
u2 / sh2.res_v2$u
#u2 / sh2.res$u

### fis obviously not equivalent but chids are
ca.res$fi[sh1,] / sh1.res_v2$fi

  ### all of these help clarify where the equivalencies exist as a bifactor analysis within the Chi2 space.
    ## some distances are equivalent and some are not.
    ## importantly I believe this inherently depends on the idea of subset CA... which side are we subsetting?

  ## equivalent
rowSums(ca.res$fi[sh1,]^2) / rowSums(sh1.res_v2$fi^2)
rowSums(ca.res$fi[sh2,]^2) / rowSums(sh2.res_v2$fi^2)
  ## not
rowSums(ca.res$fj^2) / rowSums(sh1.res_v2$fj^2)
rowSums(ca.res$fj^2) / rowSums(sh2.res_v2$fj^2)

  ## not
rowSums(ca.res$u[sh1,]^2) / rowSums(sh1.res_v2$u^2)
rowSums(ca.res$u[sh2,]^2) / rowSums(sh2.res_v2$u^2)
  ## equivalent
rowSums(ca.res$v^2) / rowSums(sh1.res_v2$v^2)
rowSums(ca.res$v^2) / rowSums(sh2.res_v2$v^2)

  ## not
rowSums(ca.res$p[sh1,]^2) / rowSums(sh1.res_v2$p^2)
rowSums(ca.res$p[sh2,]^2) / rowSums(sh2.res_v2$p^2)
  ## equivalent
rowSums(ca.res$q^2) / rowSums(sh1.res_v2$q^2)
rowSums(ca.res$q^2) / rowSums(sh2.res_v2$q^2)



### (3) now predict
  fi2.on1 <- profiles[sh2,] %*% sweep(sh1.res$v, 1, sqrt(preproc.data$w)/preproc.data$w,"*")
  p2.on1 <- sweep(sweep(fi2.on1, 2, sh1.res$d,"/"), 1, preproc.data$m[sh2],"*")
  u2.on1 <- sweep(p2.on1, 1, sqrt(preproc.data$m[sh2]),"/")

  fi1.on2 <- profiles[sh1,] %*% sweep(sh2.res$v, 1, sqrt(preproc.data$w)/preproc.data$w,"*")
  p1.on2 <- sweep(sweep(fi1.on2, 2, sh2.res$d,"/"), 1, preproc.data$m[sh1],"*")
  u1.on2 <- sweep(p1.on2, 1, sqrt(preproc.data$m[sh1]),"/")


  fi2.on.1_from.tfca <- tfca.res$pred.fi.array[sh2,,100]
  fi1.on.2_from.tfca <- tfca.res$pred.fi.array[sh1,,100]


  fi2.on1 / fi2.on.1_from.tfca
  fi1.on2 / fi1.on.2_from.tfca


  rowSums(fi1^2) / rowSums(ca.res$fi[sh1,]^2)
  rowSums(fi1.on2^2) / rowSums(ca.res$fi[sh1,]^2)

  rowSums(p1^2) / rowSums(ca.res$p[sh1,]^2)
  rowSums(p1.on2^2) / rowSums(ca.res$p[sh1,]^2)

  rowSums(u1^2) / rowSums(ca.res$u[sh1,]^2)
  rowSums(u1.on2^2) / rowSums(ca.res$u[sh1,]^2)



  rowSums(fi2^2) / rowSums(ca.res$fi[sh2,]^2)
  rowSums(fi2.on1^2) / rowSums(ca.res$fi[sh2,]^2)

  rowSums(p2^2) / rowSums(ca.res$p[sh2,]^2)
  rowSums(p2.on1^2) / rowSums(ca.res$p[sh2,]^2)

  rowSums(u2^2) / rowSums(ca.res$u[sh2,]^2)
  rowSums(u2.on1^2) / rowSums(ca.res$u[sh2,]^2)

  rowSums(u2^2) / rowSums(u2.on1^2)





  ### stupid test -- I forget the name of this technique...
  diag.median.r2s <- diag(apply(tfca.res$fi.score.cors^2,c(1,2),median))
  plot(diag.median.r2s)
  for(i in c(length(diag.median.r2s):2)){
    this.seq <- seq(1,i,1)
    abline(lm(diag.median.r2s[this.seq]~this.seq))
  }
