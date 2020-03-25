### a data utils file
## put all data transform things here.

## a radically simplified version of "ExPosition::expo.scale" that ensures we have a scale and center returned from scale every time
#' @export

ours_scale <- function(DATA, center=TRUE, scale=TRUE, impute_NA_to_mean=T){

  column_names <- colnames(DATA)
  num_columns <- ncol(DATA)
  DATA <- scale(DATA, center = center, scale = scale)

  if(is.null(attributes(DATA)$`scaled:center`)){
    attributes(DATA)$`scaled:center` <- rep(0, num_columns)
    names(attributes(DATA)$`scaled:center`) <- column_names
  }

  if(is.null(attributes(DATA)$`scaled:scale`)){
    attributes(DATA)$`scaled:scale` <- rep(1, num_columns)
    names(attributes(DATA)$`scaled:scale`) <- column_names
  }

  ### impute mean here
  DATA

}



#' @export
#'

escofier.coding <- function(DATA, center=TRUE, scale=TRUE, impute_NA_to_mean=T){

  if(is.null(colnames(DATA))){
    warning("'colnames(DATA)' were NULL. Setting to 1:ncol(DATA).")
    colnames(DATA) <- as.character(1:ncol(DATA))
  }

    ## do *not* impute mean here.
  DATA <- ours_scale(DATA, center=center, scale=scale, impute_NA_to_mean = FALSE)
  dat.col.names <- c(paste0(colnames(DATA),"-"),paste0(colnames(DATA),"+"))
  DATA <- cbind( (1-DATA)/2, (1+DATA)/2 )
  colnames(DATA) <- dat.col.names

  DATA <- as.matrix(DATA)
  attributes(DATA)$variable.map <- gsub("\\-","",gsub("\\+","",dat.col.names))


  ### impute mean here
  DATA

}


#' @export
#'

thermometer_coding <- function (DATA, mins, maxs, impute_NA_to_mean=T)
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

  DATA <- sweep(cbind(sweep(DATA, 2, maxs, "-") * -1, sweep(DATA, 2, mins, "-")),2, (maxs-mins), "/")
  colnames(DATA) <- dat.col.names
  attributes(DATA)$variable.map <- gsub("\\-", "", gsub("\\+", "", dat.col.names))


  ### impute mean here
  DATA

}





### this one is a bit ugly but it works, so I'm leaving it alone for now.
#' @export
#'

disjunctive_coding <- function(datain, impute_NA_to_mean=T){

  if(is.null(colnames(datain))){
    warning("'colnames(datain)' were NULL. Setting to 1:ncol(datain).")
    colnames(datain) <- as.character(1:ncol(datain))
  }

  data_dims <- dim(datain)
  var_names <- colnames(datain)
  ind_names <- rownames(datain)

  new.col.count <- sum(apply(datain,2,function(x){ uniq <- unique(x); length(uniq) - sum(is.na(uniq))	}))
  dataout <- matrix(0,nrow(datain),new.col.count)
  beginner <- 0
  variable.map <- new_colnames <- matrix(0, 1, 0)
  for (i in 1:data_dims[2]) {
    unique_elements <- unique(datain[, i])
    unique_no_na <- unique_elements[!is.na(unique_elements)]
    mini.mat <- matrix(0, data_dims[1], length(unique_no_na))
    for (j in 1:ncol(mini.mat)) {
      mini.mat[which(datain[, i] == unique_no_na[j]), j] <- 1
      new_colnames <- cbind(new_colnames, paste(var_names[i], ".", unique_no_na[j], sep = ""))
      variable.map <- cbind(variable.map,var_names[i])
    }

    ## here I need to be able to allow for these to remain as NA.
    these.missing <- which(rowSums(mini.mat)==0)
    if(length(these.missing)>0){
      if(impute_NA_to_mean){
        replacement <- colSums(mini.mat)/sum(colSums(mini.mat))
      }else{
        replacement <- rep(NA,ncol(mini.mat))
      }
      mini.mat[these.missing,] <- matrix(replacement,length(these.missing),length(replacement),byrow=T)
    }

    ender <- beginner + length(unique_no_na)
    dataout[, (beginner + 1):ender] <- mini.mat
    beginner <- ender
  }
  colnames(dataout) <- c(new_colnames)
  rownames(dataout) <- ind_names
  dataout <- as.matrix(dataout)
  attributes(dataout)$variable.map <- c(variable.map)


  dataout
}



#' @export
#'

ca_preproc <- function(DATA){
  Ox <- DATA/sum(DATA)
  m <- rowSums(Ox)
  w <- colSums(Ox)
  Ex <- m %o% w
  Zx <- Ox - Ex
  weightedZx <- sweep(sweep(Zx,1,sqrt(m),"/"),2,sqrt(w),"/")
  return( list(m=m,w=w,Zx=Zx,Ox=Ox,Ex=Ex,weightedZx=weightedZx) )
}




#### need to revisit this one...
# categorical = "n"
# continuous= "z"
# ordinal = "o"
# counts/non-negatives = "f"; this effectively does nothing...
# counts/non-negatives to be normalized to 1 by column = "f1"
# do nothing - "x"

### DATA must be a data.frame for things to work generally, else it explodes.
### I generally discourage the use of this function---so much so that I may remove it entirely.
### each column must be correctly formed.

mixed_data_coding <- function(DATA, column.type = rep("x",ncol(DATA)), impute.NA.to.mean=F){

  if(class(DATA) != "data.frame"){
    stop("'DATA' must be a data.frame. To recode various types of data see make.data.nominal(), thermometer.coding(), and escofier.coding().")
  }

  if(is.null(colnames(DATA))){
    warning("'colnames(DATA)' were NULL. Setting to 1:ncol(DATA).")
    colnames(DATA) <- as.character(1:ncol(DATA))
  }
  # if(is.null(rownames(DATA))){
  #   warning("'rownames(DATA)' were NULL. Setting to 1:nrow(DATA).")
  #   rownames(DATA) <- as.character(1:nrow(DATA))
  # }

  # orig.colnames <- names(column.type) <- colnames(DATA)
  names(column.type) <- colnames(DATA)

  possible.types <- c("n","z","o","f","f1","x")
  if(any(!(column.type %in% possible.types))){
    warning("Unrecognized 'column.type'. Changing unrecognized types to 'x'")
    column.type[which(!(column.type %in% possible.types))] <-"x"
  }

  ## so I need to break these out into their types, do what needs to be done with the blocks, then somehow re-arrange...
  ## based on the variable maps produced, and the original column names, I should be able to establish an order...

  map.types <- lapply(possible.types,function(x){ column.type %in% x })
  ## may not be necessary.
  names(map.types) <- possible.types


  variable.map <- type.map <- c() #I'm lazy.
  DATA.out <- matrix(NA,nrow(DATA),0) #I'm so lazy.
  # I guess as I go through these I could use a flag, and then allocate the necessary amount of columns..?.

  if(any(map.types$n)){
    n.mat <- as.matrix(DATA[,which(map.types$n)])
    ## for safety...
    rownames(n.mat) <- rownames(DATA)
    colnames(n.mat) <- colnames(DATA[,which(map.types$n)])
    n.mat <- make.data.nominal(n.mat,impute.NA.to.mean=impute.NA.to.mean)
    type.map <- c(type.map,rep("n",ncol(n.mat)))
    variable.map <- c(variable.map,attributes(n.mat)$variable.map)
    DATA.out <- cbind(DATA.out,n.mat)
  }
  if(any(map.types$z)){
    z.mat <- as.matrix(DATA[,which(map.types$z)])
    ## for safety...
    rownames(z.mat) <- rownames(DATA)
    colnames(z.mat) <- colnames(DATA[,which(map.types$z)])
    z.mat <- escofier.coding(z.mat)
    type.map <- c(type.map,rep("z",ncol(z.mat)))
    variable.map <- c(variable.map,attributes(z.mat)$variable.map)
    DATA.out <- cbind(DATA.out,z.mat)
  }
  if(any(map.types$o)){
    o.mat <- as.matrix(DATA[,which(map.types$o)])
    ## for safety...
    rownames(o.mat) <- rownames(DATA)
    colnames(o.mat) <- colnames(DATA[,which(map.types$o)])
    o.mat <- thermometer.coding(o.mat)
    type.map <- c(type.map,rep("o",ncol(o.mat)))
    variable.map <- c(variable.map,attributes(o.mat)$variable.map)
    DATA.out <- cbind(DATA.out,o.mat)
  }
  if(any(map.types$f1)){
    f1.mat <- as.matrix(DATA[,which(map.types$f1)])
    ## for safety...
    rownames(f1.mat) <- rownames(DATA)
    colnames(f1.mat) <- colnames(DATA[,which(map.types$f1)])
    f1.mat <- sweep(f1.mat,2,colSums(f1.mat),"/")
    #attributes(f1.mat)$variable.map <- colnames(f1.mat)
    type.map <- c(type.map,rep("f1",ncol(f1.mat)))
    variable.map <- c(variable.map,colnames(f1.mat))
    DATA.out <- cbind(DATA.out,f1.mat)
  }
  if(any(map.types$f)){
    f.mat <- as.matrix(DATA[,which(map.types$f)])
    ## for safety...
    rownames(f.mat) <- rownames(DATA)
    colnames(f.mat) <- colnames(DATA[,which(map.types$f)])
    type.map <- c(type.map,rep("f",ncol(f.mat)))
    variable.map <- c(variable.map,colnames(f.mat))
    DATA.out <- cbind(DATA.out,f.mat)
  }
  if(any(map.types$x)){
    x.mat <- as.matrix(DATA[,which(map.types$x)])
    ## for safety...
    rownames(x.mat) <- rownames(DATA)
    colnames(x.mat) <- colnames(DATA[,which(map.types$x)])
    attributes(x.mat)$variable.map <- colnames(x.mat)
    type.map <- c(type.map,rep("x",ncol(x.mat)))
    variable.map <- c(variable.map,colnames(x.mat))
    DATA.out <- cbind(DATA.out,x.mat)
  }

  # for safety
  attributes(DATA.out)$variable.map <- variable.map
  attributes(DATA.out)$type.map <- type.map
  rownames(DATA.out) <- rownames(DATA)
  return(DATA.out)

}
