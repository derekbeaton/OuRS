data(SNPS)

  ## will need these.
SNPs.nom <- make.data.nominal(SNPS)
ca.preproc.res <- ca.preproc(SNPs.nom)
profiles <- sweep(ca.preproc.res$Ox,1,ca.preproc.res$m,"/")

  ## do MCD
snps.res <- cat.mcd(SNPs.nom,make.data.disjunctive = F,num.subsets = 500,alpha = .75)

  ## do corrmax
snps.corrmax.res <- cat.corrmax(ca.preproc.res$weightedZx, loadings = snps.res$cov$loadings,singular.values = snps.res$cov$singular.values)

  ## do bootstrap
snps.boot.res <- cat.boot.sup.u(make.data.nominal(SNPS),snps.res$cov$loadings,snps.res$cov$singular.values,iters = 100)
