mahal.from.ca <- function(profiles, row.weights, col.weights, loadings, singular.values){

  ## NOT EFFICIENT. MAKE MORE EFFICIENT

  #fi.scores <- profiles %*% 
  u.scores <- diag(sqrt(row.weights)) %*% (profiles %*% (diag(sqrt(col.weights)*(1/col.weights)) %*% loadings) %*% diag(1/singular.values))
  ## can I make this easier? I know I can speed up computation with ordinary multiplications.
  return( rowSums(( u.scores )^2) )

}
