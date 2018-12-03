# categorical = "n"
# continuous= "z"
# ordinal = "o"
# counts/non-negatives = "f"; this effectively does nothing...
# counts/non-negatives to be normalized to 1 by column = "f1"
# do nothing - "x"

  ### DATA must be a data.frame for things to work generally, else it explodes.
  ### I generally discourage the use of this function---so much so that I may remove it entirely.
    ### each column must be correctly formed.
mixed.data.coding <- function(DATA, column.type = rep("x",ncol(DATA)), impute.NA.to.mean=F){

  if(class(DATA) != "data.frame"){
    stop("'DATA' must be a data.frame. To recode various types of data see make.data.nominal(), thermometer.coding(), and escofier.coding().")
  }

  if(is.null(colnames(DATA))){
    warning("'colnames(DATA)' were NULL. Setting to 1:ncol(DATA).")
    colnames(DATA) <- as.character(1:ncol(DATA))
  }
  if(is.null(rownames(DATA))){
    warning("'rownames(DATA)' were NULL. Setting to 1:nrow(DATA).")
    rownames(DATA) <- as.character(1:nrow(DATA))
  }

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
  DATA.out <- matrix(NULL,nrow(DATA),0) #I'm so lazy.
      # I guess as I go through these I could use a flag, and then allocate the necessary amount of columns..?.

  if(any(map.types$n)){
    n.mat <- as.matrix(DATA[,which(map.types$n)])
      ## for safety...
    colnames(n.mat) <- colnames(DATA[,which(map.types$n)])
    n.mat <- make.data.nominal(n.mat,impute.NA.to.mean=impute.NA.to.mean)
    type.map <- c(type.map,rep("n",ncol(n.mat)))
    variable.map <- c(variable.map,attributes(n.mat)$variable.map)
    DATA.out <- cbind(DATA.out,n.mat)
  }
  if(any(map.types$z)){
    z.mat <- as.matrix(DATA[,which(map.types$z)])
    ## for safety...
    colnames(z.mat) <- colnames(DATA[,which(map.types$z)])
    z.mat <- escofier.coding(z.mat)
    type.map <- c(type.map,rep("z",ncol(z.mat)))
    variable.map <- c(variable.map,attributes(z.mat)$variable.map)
    DATA.out <- cbind(DATA.out,z.mat)
  }
  if(any(map.types$o)){
    o.mat <- as.matrix(DATA[,which(map.types$o)])
    ## for safety...
    colnames(o.mat) <- colnames(DATA[,which(map.types$o)])
    o.mat <- thermometer.coding(o.mat)
    type.map <- c(type.map,rep("o",ncol(o.mat)))
    variable.map <- c(variable.map,attributes(o.mat)$variable.map)
    DATA.out <- cbind(DATA.out,o.mat)
  }
  if(any(map.types$f1)){
    f1.mat <- as.matrix(DATA[,which(map.types$f1)])
    ## for safety...
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
    colnames(f.mat) <- colnames(DATA[,which(map.types$f)])
    type.map <- c(type.map,rep("f",ncol(f.mat)))
    variable.map <- c(variable.map,colnames(f.mat))
    DATA.out <- cbind(DATA.out,f.mat)
  }
  if(any(map.types$x)){
    x.mat <- as.matrix(DATA[,which(map.types$x)])
    ## for safety...
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
