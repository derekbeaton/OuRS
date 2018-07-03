thermometer.coding <- function(DATA, mins, maxs, norm.to.one = T){
  
  if(missing(mins)){
    mins <- apply(DATA,2,min,na.rm=T)
  }else{
    
    if(length(mins)==ncol(DATA)){
      min.test <- mins > apply(DATA,2,min,na.rm=T)
      if(any(min.test)){
        warning("Some inputted minimums are greater than minimums in the data. We will replace those 'mins' with their respective minimum in the data")
        mins[which(min.test)] <- apply(DATA[,which(min.test)],2,min,na.rm=T)
      }
    }else{
      mins <- apply(DATA,2,min,na.rm=T)
    }
    
  }
  
  
  if(missing(maxs)){
    maxs <- apply(DATA,2,max,na.rm=T)
  }else{
    
    if(length(maxs)==ncol(DATA)){
      max.test <- maxs < apply(DATA,2,max,na.rm=T)
      if(any(max.test)){
        warning("Some inputted maximums are smaller than maximums in the data. We will replace those 'maxs' with their respective maximum in the data")
        maxs[which(max.test)] <- apply(DATA[,which(max.test)],2,max,na.rm=T)
      }
    }else{
      maxs <- apply(DATA,2,max,na.rm=T)
    }
    
  }
  
  dat.col.names <- c(paste0(colnames(DATA),"+"),paste0(colnames(DATA),"-"))
  
  
  from.mins <- sweep(DATA,2,mins,"-")
  
  if(norm.to.one){ ## these should be normed so that the variables = 1.
    DATA <- cbind(sweep( from.mins ,2,maxs,"/"), sweep( sweep(from.mins,2,maxs,"-") * -1,2,maxs,"/"))
  }else{
    DATA <- cbind( from.mins ,  sweep(from.mins,2,maxs,"-") * -1)
  }
  colnames(DATA) <- dat.col.names
  
  
  DATA <- as.matrix(DATA)
  attributes(DATA)$variable.map <- gsub("\\-","",gsub("\\+","",dat.col.names))
  
  return(DATA)
}
