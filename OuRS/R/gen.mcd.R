  ## data types will only allow for "cat" = categorical, "ord" for ordinal", "con" for continuous and "frq" for frequency

gen.mcd <- function(data, column.types=rep("cat",ncol(data)), alpha=.75, num.subsets=500, max.total.iters=num.subsets*20, top.sets.percent=.05, tol=.Machine$double.eps){
}
