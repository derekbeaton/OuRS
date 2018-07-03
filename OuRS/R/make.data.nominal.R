## I need a much better way of doing this...
  ## also I need to provide the option of filling NAs in with the barycenters/means.
  ### that should not be the default.
make.data.nominal <- function(datain,impute.NA.to.mean=T){


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
      if(impute.NA.to.mean){
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
  

  return(dataout)
}
