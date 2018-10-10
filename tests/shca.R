### some tests.

data("dfw.beer.survey")
nom.dat <- make.data.nominal(DFW.beer.survey)

ca.res <- ca(nom.dat)
preproc.data <- ca.preproc(nom.dat)
profiles <- sweep(preproc.data$Ox,1,preproc.data$m,"/")


sh1 <- sort(sample( 1:nrow(nom.dat), round(nrow(nom.dat)/2) ))
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
  u1 <- sweep(sweep(fi1, 2, sh1.res$d,"/"), 1, sqrt(preproc.data$m[sh1]),"*")

  fi2 <- profiles[sh2,] %*% sweep(sh2.res$v, 1, sqrt(preproc.data$w)/preproc.data$w,"*")
  u2 <- sweep(sweep(fi2, 2, sh2.res$d,"/"), 1, sqrt(preproc.data$m[sh2]),"*")


