### a data utils file
## put all data transform things here.

## a radically simplified version of "ExPosition::expo.scale" that ensures we have a scale and center returned from scale every time
#' @title Variant of \code{scale}: scaling and centering of matrix-like objects
#' 
#' @description Performs the same functionality as \code{\link{scale}} but always returns the "scale:center" and "scale:scale" attributes. Also allows mean imputation for \code{NA}
#' 
#' @details 
#' If no centering is performed (i.e., set to \code{FALSE}), then \code{attributes()`scale:center` <- rep(0, ncol(DATA))}
#' If no scaling is performed (i.e., set to \code{FALSE}), then \code{attributes()`scale:scale` <- rep(1, ncol(DATA))}
#' If \code{impute_NA_to_mean = TRUE} then all \code{NA} will be replaced with the column-wise mean *after* centering/scaling
#' 
#' @param DATA a numeric matrix
#' @param center logical or numeric (see \code{\link{scale}}). Default is \code{TRUE} which centers the columns (e.g., when \code{TRUE} substract the mean of a column from its respective column)
#' @param scale logical or numeric (see \code{\link{scale}}). Default is \code{TRUE} which scales the columns (e.g., when \code{TRUE} divide a column by its respective standard deviation or scaling factor)
#' @param impute_NA_to_mean a logical (boolean). Default is \code{FALSE}. If \code{TRUE} all \code{NA} will be replaced with the column-wise mean *after* centering and scaling.
#' 
#' @seealso \code{\link{scale}}
#' @author Derek Beaton
#' @export

ours_scale <- function(DATA, center=TRUE, scale=TRUE, impute_NA_to_mean=F){

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

  if(impute_NA_to_mean){
    DATA <- apply(DATA, 2, function(x){ x[is.na(x)] <- mean(x, na.rm = T); x })
  }
  
  DATA

}



#' @title Escofier's approach to "data-doubling" for continuous data
#' 
#' @description Performs data-doubling for continuous data so that the data behave like disjunctive data.
#' This function makes use of \code{\link{ours_scale}} to center and/or scale (presumably) continous data.
#' 
#' @details 
#' Data are returned in a "doubled" way where each column is represented twice as \eqn{[\frac{1-x}{2}, \frac{1+x}{2}]} where \code{x} is a centered and/or scaled vector.
#' If \code{impute_NA_to_mean = TRUE} then all \code{NA} will be replaced with the column-wise mean *after* centering/scaling
#' 
#' @param DATA a numeric matrix
#' @param center logical or numeric (see \code{\link{scale}}). Default is \code{TRUE} which centers the columns (e.g., when \code{TRUE} substract the mean of a column from its respective column)
#' @param scale logical or numeric (see \code{\link{scale}}). Default is \code{TRUE} which scales the columns (e.g., when \code{TRUE} divide a column by its respective standard deviation or scaling factor)
#' @param impute_NA_to_mean a logical (boolean). Default is \code{FALSE}. If \code{TRUE} all \code{NA} will be replaced with the column-wise mean *after* centering and scaling.
#' 
#' @references
#' 
#' @seealso \code{\link{scale}}, \code{\link{disjunctive_coding}}, and \code{\link{thermometer_coding}}
#' @author Derek Beaton
#' @export

escofier_coding <- function(DATA, center=TRUE, scale=TRUE, impute_NA_to_mean=F){

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
  if(impute_NA_to_mean){
    DATA <- apply(DATA, 2, function(x){ x[is.na(x)] <- mean(x, na.rm = T); x })
  }

  DATA

}


#' @title "Thermometer" approach to "data-doubling" for ordinal data
#' 
#' @description Performs data-doubling for ordinal data so that the data behave like disjunctive data.
#' 
#' @details 
#' Data are returned in a "doubled" way where each column is represented twice as \eqn{[\frac{1-x}{2}, \frac{1+x}{2}]} where \code{x} is a centered and/or scaled vector.
#' If \code{impute_NA_to_mean = TRUE} then all \code{NA} will be replaced with the column-wise mean *after* transformation.
#' 
#' @param DATA a numeric matrix
#' @param mins (optional) numeric vector of length \code{ncol(DATA)}. Contains the expected minimum value per column. If none provided then they are computed as \code{apply(DATA, 2, min, na.rm = T)}
#' @param maxs (optional) numeric vector of length \code{ncol(DATA)}. Contains the expected maximum value per column. If none provided then they are computed as \code{apply(DATA, 2, max, na.rm = T)}
#' @param impute_NA_to_mean a logical (boolean). Default is \code{FALSE}. If \code{TRUE} all \code{NA} will be replaced with the column-wise mean *after* transformation.
#' 
#' @references
#' 
#' @seealso \code{\link{escofier_coding}} and \code{\link{disjunctive_coding}}
#' @author Derek Beaton
#' @export

thermometer_coding <- function (DATA, mins, maxs, impute_NA_to_mean=F)
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
  if( any(is.na(maxs)) | any(is.infinite(maxs)) | any(is.nan(maxs)) | any(is.null(maxs)) ){
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
  if( any(is.na(mins)) | any(is.infinite(mins)) | any(is.nan(mins)) | any(is.null(mins)) ){
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
  if(impute_NA_to_mean){
    DATA <- apply(DATA, 2, function(x){ x[is.na(x)] <- mean(x, na.rm = T); x })
  }

  DATA

}





### this one is a bit ugly but it works, so I'm leaving it alone for now.
#' @title Complete disjunctive coding
#' 
#' @description Complete disjunctive coding for categorical (nominal) data
#' 
#' @details If \code{impute_NA_to_mean = TRUE} then all \code{NA} will be replaced with the column-wise mean *after* transformation
#' The presence of a level is indicated with a '1', where the absence of a level is '0'. Imputed values are the mean proportion of 0s & 1s and thus (0,1).
#' 
#' @param DATA a matrix where each column is a variable, and each value within the columns is assumed to be a level within that variable
#' @param impute_NA_to_mean a logical (boolean). Default is \code{FALSE}. If \code{TRUE} all \code{NA} will be replaced with the column-wise mean *after* transformation.
#' 
#' @references
#' 
#' @seealso \code{\link{escofier_coding}} and \code{\link{thermometer_coding}}
#' @author Derek Beaton
#' @export

disjunctive_coding <- function(DATA, impute_NA_to_mean=F){

  if(is.null(colnames(DATA))){
    # warning("'colnames(DATA)' were NULL. Setting to 1:ncol(DATA).")
    colnames(DATA) <- as.character(1:ncol(DATA))
  }

  data_dims <- dim(DATA)
  var_names <- colnames(DATA)
  ind_names <- rownames(DATA)

  new.col.count <- sum(apply(DATA,2,function(x){ uniq <- unique(x); length(uniq) - sum(is.na(uniq))	}))
  dataout <- matrix(0,nrow(DATA),new.col.count)
  beginner <- 0
  variable.map <- new_colnames <- matrix(0, 1, 0)
  for (i in 1:data_dims[2]) {
    unique_elements <- unique(DATA[, i])
    unique_no_na <- unique_elements[!is.na(unique_elements)]
    mini.mat <- matrix(0, data_dims[1], length(unique_no_na))
    for (j in 1:ncol(mini.mat)) {
      mini.mat[which(DATA[, i] == unique_no_na[j]), j] <- 1
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



# categorical = "n"
# continuous, centered only = "c"
# continuous, centered and scaled = "z"
# ordinal = "o"
# do nothing = "x"

### DATA must be a data.frame for things to work generally, else it explodes.
### I generally discourage the use of this function---so much so that I may remove it entirely.
### each column must be correctly formed.


### this one is a bit ugly but it works, so I'm leaving it alone for now.
#' @title Mixed data coding
#' 
#' @description Transforms mixtures of data contained within a data.frame
#' 
#' @details For each type of data, this function cals into the respective transformation function. 
#' \item{"n"} {is categorical data and calls \code{disjunctive_coding}}
#' \item{"c"} {is continuous data that should only be centered and calls \code{escofier_coding} with \code{center = TRUE} and \code{scale = FALSE}}
#' \item{"z"} {is continuous data that should be centered and scaled, and calls \code{escofier_coding} with \code{center = TRUE} and \code{scale = TRUE}}
#' \item{"o"} {is ordinal data and calls \code{thermometer_coding}}
#' \item{"x"} {is 'do nothing' and returns the data exactly as they are}
#' The output (transformed) data are returned in a different order than they are input. The order matches the types as: "n", "c", "z", "o", and then "x"
#' In all cases, values are derived from the data and other parameters of those functions are not available. For examples: you cannot pass \code{min} to \code{thermometer_coding} here. Instead, the observed minimums will be used
#' \code{impute_NA_to_mean} is applied globally and passed into each function exactly as it is used here.
#' If \code{impute_NA_to_mean = TRUE} then all \code{NA} will be replaced with the column-wise mean *after* transformation
#' 
#' @param DATA a matrix where each column is a variable, and each value within the columns is assumed to be a level within that variable
#' @param column.type a character vector of length \code{ncol(DATA)}. Each element indicates the data type of its respective column: "n" is categorical/nominal, "c" is continuous (only center the data), "z" is continuous (center and scale the data) "o" is ordinal, and "x" means 'do nothing to this column'. Each column will be transformed according to \code{disjunctive_coding} (categorical), \code{escofier_coding} (continuous), \code{thermometer_coding} (ordinal), or no transformation.
#' @param impute_NA_to_mean a logical (boolean). Default is \code{FALSE}. If \code{TRUE} all \code{NA} will be replaced with the column-wise mean *after* transformation.
#' 
#' 
#' @seealso \code{\link{disjunctive_coding}}, \code{\link{thermometer_coding}}, and \code{\link{escofier_coding}}
#' @author Derek Beaton
#' @export

mixed_data_coding <- function(DATA, column.type = rep("x",ncol(DATA)), impute_NA_to_mean=F){

  if(!is.data.frame(DATA)){
    stop("mixed_data_coding: 'DATA' must be a data.frame. To recode various types of data see disjunctive_coding(), thermometer_coding(), and escofier_coding().")
  }

  if(is.null(colnames(DATA))){
    colnames(DATA) <- as.character(1:ncol(DATA))
  }

  names(column.type) <- colnames(DATA)

  possible.types <- c("n","c","z","o","x")
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
    n.mat <- disjunctive_coding(n.mat,impute_NA_to_mean=impute_NA_to_mean)
    type.map <- c(type.map,rep("n",ncol(n.mat)))
    variable.map <- c(variable.map,attributes(n.mat)$variable.map)
    DATA.out <- cbind(DATA.out,n.mat)
  }
  if(any(map.types$c)){
    c.mat <- as.matrix(DATA[,which(map.types$z)])
    ## for safety...
    rownames(c.mat) <- rownames(DATA)
    colnames(c.mat) <- colnames(DATA[,which(map.types$c)])
    c.mat <- escofier_coding(c.mat, TRUE, FALSE, impute_NA_to_mean = impute_NA_to_mean)
    type.map <- c(type.map,rep("c",ncol(c.mat)))
    variable.map <- c(variable.map,attributes(c.mat)$variable.map)
    DATA.out <- cbind(DATA.out,c.mat)
  }
  if(any(map.types$z)){
    z.mat <- as.matrix(DATA[,which(map.types$z)])
    ## for safety...
    rownames(z.mat) <- rownames(DATA)
    colnames(z.mat) <- colnames(DATA[,which(map.types$z)])
    z.mat <- escofier_coding(z.mat, TRUE, TRUE, impute_NA_to_mean = impute_NA_to_mean)
    type.map <- c(type.map,rep("z",ncol(z.mat)))
    variable.map <- c(variable.map,attributes(z.mat)$variable.map)
    DATA.out <- cbind(DATA.out,z.mat)
  }
  if(any(map.types$o)){
    o.mat <- as.matrix(DATA[,which(map.types$o)])
    ## for safety...
    rownames(o.mat) <- rownames(DATA)
    colnames(o.mat) <- colnames(DATA[,which(map.types$o)])
    o.mat <- thermometer_coding(o.mat, impute_NA_to_mean = impute_NA_to_mean)
    type.map <- c(type.map,rep("o",ncol(o.mat)))
    variable.map <- c(variable.map,attributes(o.mat)$variable.map)
    DATA.out <- cbind(DATA.out,o.mat)
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


#' @title Correspondence analysis (CA) preprocessing
#' @description Preprocessing/preparing a data matrix for correspondence analysis
#' @details A data matrix processed akin to the way \eqn{\chi^2} is computed (deviations from independence)
#' 
#' @param DATA a numeric matrix
#' 
#' @return a list with 6 elements
#' \item{Ox:} {The observed values}
#' \item{m:} {Row probabilities}
#' \item{w:} {Column probabilities}
#' \item{Ex:} {The expected values}
#' \item{Zx:} {The deviations values}
#' \item{weightedZx:} {The deviations divided by the square root of the row and column probabilities}
#' 
#' @author Derek Beaton
#' @export

ca_preproc <- function(DATA){
  Ox <- DATA/sum(DATA)
  m <- rowSums(Ox)
  w <- colSums(Ox)
  Ex <- m %o% w
  Zx <- Ox - Ex
  weightedZx <- sweep(sweep(Zx,1,sqrt(m),"/"),2,sqrt(w),"/")
  return( list(m=m,w=w,Zx=Zx,Ox=Ox,Ex=Ex,weightedZx=weightedZx) )
}