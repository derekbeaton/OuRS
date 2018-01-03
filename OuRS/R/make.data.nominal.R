make.data.nominal <- function(datain){
  
  
  num.nas <- sum(is.na(datain))
  if (num.nas > 0) {
    print(paste("You have ", num.nas, " NAs in your data. makeNominalData automatically imputes NA with the mean of the columns.", sep = ""))
  }
  
  data_dims <- dim(datain)
  var_names <- colnames(datain)
  ind_names <- rownames(datain)
  
  new.col.count <- sum(apply(datain,2,function(x){ uniq <- unique(x); length(uniq) - sum(is.na(uniq))	}))
  dataout <- matrix(0,nrow(datain),new.col.count)
  beginner <- 0
  new_colnames <- matrix(0, 1, 0)
  for (i in 1:data_dims[2]) {
    unique_elements <- unique(datain[, i])
    unique_no_na <- unique_elements[!is.na(unique_elements)]
    mini.mat <- matrix(0, data_dims[1], length(unique_no_na))
    for (j in 1:ncol(mini.mat)) {
      mini.mat[which(datain[, i] == unique_no_na[j]), j] <- 1
      new_colnames <- cbind(new_colnames, paste(var_names[i], ".", unique_no_na[j], sep = ""))
    }
    barycenter <- colSums(mini.mat)/sum(colSums(mini.mat))
    fill_in <- repmat(t(as.matrix(barycenter)), length(which(rowSums(mini.mat) == 0)), 1)
    mini.mat[which(rowSums(mini.mat) == 0), ] <- fill_in
    ender <- beginner + length(unique_no_na)
    dataout[, (beginner + 1):ender] <- mini.mat
    beginner <- ender
  }
  colnames(dataout) <- new_colnames
  rownames(dataout) <- ind_names
  return(as.matrix(dataout))
}