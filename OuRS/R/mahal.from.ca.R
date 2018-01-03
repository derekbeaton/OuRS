mahal.from.ca <- function(data, row.weights, col.weights, loadings, singular.values){

  ## NOT EFFICIENT. MAKE MORE EFFICIENT

  p.scores <- diag(sqrt(row.weights)) %*% (data %*% (diag(sqrt(col.weights)*(1/col.weights)) %*% loadings) %*% diag(1/singular.values))
  ## can I make this easier? I know I can speed up computation with ordinary multiplications.
  return( rowSums(( p.scores )^2) )

}
