#escofier.coding
  # allow no norm, z norm, etc... through expo.scale

  ## this one is real easy.
escofier.coding <- function(DATA, center=T, scale="SS1"){

  DATA <- expo.scale(DATA,center=center,scale=scale)
  dat.col.names <- c(paste0(colnames(DATA),"-"),paste0(colnames(DATA),"+"))
  DATA <- cbind( (1-DATA)/2, (1+DATA)/2 )
  colnames(DATA) <- dat.col.names

  DATA <- as.matrix(DATA)
  attributes(DATA)$variable.map <- gsub("\\-","",gsub("\\+","",dat.col.names))
  
  return(DATA)
}
