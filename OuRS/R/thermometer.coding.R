thermometer.coding <- function (DATA, mins, maxs)
{
  if (is.null(colnames(DATA))) {
    warning("'colnames(DATA)' were NULL. Setting to 1:ncol(DATA).")
    colnames(DATA) <- as.character(1:ncol(DATA))
  }

  dat.col.names <- c(paste0(colnames(DATA), "+"), paste0(colnames(DATA), "-"))

  ## these need a lot of checks:
    ### and these need to be cleaned up substantially...
  if(missing(maxs)){
    maxs <- apply(DATA, 2, max, na.rm = T)
  }
  if(length(maxs)!=ncol(DATA)){
    maxs <- apply(DATA, 2, max, na.rm = T)
  }
  if(any( maxs < apply(DATA, 2, max, na.rm = T))){
    maxs <- apply(DATA, 2, max, na.rm = T)
  }

  if(missing(mins)){
    mins <- apply(DATA, 2, min, na.rm = T)
  }
  if(length(mins)!=ncol(DATA)){
    mins <- apply(DATA, 2, min, na.rm = T)
  }
  if(any( mins > apply(DATA, 2, min, na.rm = T))){
    mins <- apply(DATA, 2, min, na.rm = T)
  }
  ## clean the above.

  DATAOUT <- sweep(cbind(sweep(DATA, 2, maxs, "-") * -1, sweep(DATA, 2, mins, "-")),2, (maxs-mins), "/")
  colnames(DATAOUT) <- dat.col.names
  attributes(DATAOUT)$variable.map <- gsub("\\-", "", gsub("\\+", "", dat.col.names))
  return(DATAOUT)
}
