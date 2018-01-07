data(SNPS)

snps.res <- cat.mcd(SNPS,make.data.disjunctive = T,num.subsets = 500,alpha = .5)
ca.preproc.res <- ca.preproc(make.data.nominal(SNPS))

corrmax.res <- cat.corrmax(ca.preproc.res$weightedZx, loadings = snps.res$cov$loadings,singular.values = snps.res$cov$singular.values)
