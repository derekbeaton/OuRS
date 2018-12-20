

cont.multi.boxplot <- function(
  data,
  center=T, scale=F,
  IQR = 1.5
){

  # scaling the data and creating rownames, if necessary
  mydat <- scale(data,center=center,scale=scale)
  if(is.null(rownames(mydat))){rownames(mydat) <- 1:nrow(mydat)}

  # identifying points beyond the whiskers (= outliers)
  mydat.find.outliers <- apply(mydat,2,function(i){boxplot.stats(i,coef=IQR)$out})

  if(length(mydat.find.outliers) == 0){ # no outliers for any variable

    outliers.formatted <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("SUBJECT","VARIABLE","OUT.VALUE"))
    outliers.by.subject <- data.frame(SUBJECT=rownames(mydat),COUNT=0)

  }else{ # at least one variable

    # formatting output to get clean data.frame (from list - and remove variables without outliers)

    if(any(sapply(mydat.find.outliers,function(i){length(i)==0}))){
      outliers.formatted <- lapply(mydat.find.outliers[-which(sapply(mydat.find.outliers,function(i){length(i)==0}))],
                                   function(i){data.frame(SUBJECT=names(i),OUT.VALUE=i)})

    }else{ # no variables without outliers
      outliers.formatted <- lapply(mydat.find.outliers,function(i){data.frame(SUBJECT=names(i),OUT.VALUE=i)})
    }
    outliers.formatted <- mapply(`[<-`, outliers.formatted, 'VARIABLE', value = names(outliers.formatted), SIMPLIFY = FALSE)



    # transforming into long format
    outliers.formatted.df <- data.frame()
    for(i in outliers.formatted){ # any chance it's not a list?
      outliers.formatted.df <- rbind(outliers.formatted.df,i)
    }
    outliers.formatted <- outliers.formatted.df


    # number of outliers per obs
    outliers.by.subject <- as.data.frame(table(outliers.formatted$SUBJECT))   # should add back all the 0s
    colnames(outliers.by.subject) <- c("SUBJECT","COUNT")
    outliers.by.subject <- merge(outliers.by.subject,data.frame(SUBJECT=rownames(mydat)),all=T)
    outliers.by.subject$COUNT <-ifelse(is.na(outliers.by.subject$COUNT),0,outliers.by.subject$COUNT)
    outliers.by.subject <- outliers.by.subject[order(outliers.by.subject$COUNT,decreasing = T),]


  } # end of at least one variable

  rownames(outliers.by.subject) <- outliers.by.subject$SUBJECT
  outliers.by.subject <- outliers.by.subject[,-1,drop=F]

  # transforming mydat into long format - to merge with outliers
  mydat.long <- data.frame()
  for(i in colnames(mydat)){
    mydat.long <- rbind(mydat.long,data.frame(SUBJECT=rownames(mydat),VARIABLE=i,VALUE=mydat[,i]))
  }
  rownames(mydat.long) <- NULL

  # merging long mydat with outliers and creating label variable
  mydat.merged <- merge(mydat.long,outliers.formatted,all=T)

  # reformatting outliers.formatted for output
  outliers.formatted <- outliers.formatted[order(as.character(outliers.formatted$SUBJECT),-outliers.formatted$OUT.VALUE),c("SUBJECT","VARIABLE","OUT.VALUE")]

  outliers.out <- list(input.parameters = list(center=center, scale=scale, IQR = IQR),
                       set.outliers = outliers.formatted,
                       outliers.by.subject = outliers.by.subject,
                       data.with.outliers = mydat.merged)

  return(outliers.out)
  }
