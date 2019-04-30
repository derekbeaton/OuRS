thermometer.coding <- function(DATA){
  ## radically re-written and consequently simplified.
  ## the previous version was not correct.

  if(is.null(colnames(DATA))){
    warning("'colnames(DATA)' were NULL. Setting to 1:ncol(DATA).")
    colnames(DATA) <- as.character(1:ncol(DATA))
  }

  ## enforce proper thermometer: subtract the maxs and make this easy.
  DATA <- apply(DATA, 2, function(x){ max(x, na.rm = T)-x })
  mins <- apply(DATA, 2, min, na.rm=T)
  maxs <- apply(DATA, 2, max, na.rm=T)
  dat.col.names <- c(paste0(colnames(DATA),"+"),paste0(colnames(DATA),"-"))

  DATA <- as.matrix(
    cbind(
      sweep(sweep(DATA,2,maxs,"-")*-1,2,maxs,"/"),
      sweep(sweep(DATA,2,mins,"-"),2,maxs,"/")
      )
    )
  colnames(DATA) <- dat.col.names
  attributes(DATA)$variable.map <- gsub("\\-","",gsub("\\+","",dat.col.names))

  return(DATA)
}
