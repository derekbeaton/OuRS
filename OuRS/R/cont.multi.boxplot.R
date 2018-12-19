

########################################################################################################################################
########################################################################################################################################
# Univariate Assessment for OuRS package:
# A series of boxplots at some IQR; outliers are those beyond the whiskers; list number of outliers for each participant
########################################################################################################################################
########################################################################################################################################



### make lower margin adjust with length of variable names!?



cont.multi.boxplot <- function(
  data,
  center=T, scale=F,
  IQR = 1.5,
  plot.box=T
){
  
  # scaling the data, if necessary
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
    
    
    
    # need to assume some will be empty... but for now
    outliers.formatted.df <- data.frame()
    for(i in outliers.formatted){ # any chance it's not a list?
      outliers.formatted.df <- rbind(outliers.formatted.df,i)
    }
    outliers.formatted <- outliers.formatted.df

    # number of outliers per participant
    
    outliers.by.subject <- as.data.frame(table(outliers.formatted$SUBJECT))   # should add back all the 0s
    colnames(outliers.by.subject) <- c("SUBJECT","COUNT")
    outliers.by.subject <- merge(outliers.by.subject,data.frame(SUBJECT=rownames(mydat)),all=T)
    outliers.by.subject$COUNT <-ifelse(is.na(outliers.by.subject$COUNT),0,outliers.by.subject$COUNT)
    outliers.by.subject <- outliers.by.subject[order(outliers.by.subject$COUNT,decreasing = T),]
    
    
  } # end of at least one variable
  
  
  #####

  all.outliers <- outliers.formatted[order(as.character(outliers.formatted$SUBJECT),-outliers.formatted$OUT.VALUE),c("VARIABLE","OUT.VALUE")]
  
  #####
  outliers.by.count <- as.data.frame(table(outliers.by.subject$COUNT))
  colnames(outliers.by.count) <- c("OUTLIERS.PER.OBS","COUNT")

  rownames(outliers.by.subject) <- outliers.by.subject$SUBJECT
  outliers.by.subject <- outliers.by.subject[,-1,drop=F]
  
  
  outliers.out <- list(set.outliers = all.outliers,
                       outliers.by.subject = outliers.by.subject,
                       outliers.by.count = outliers.by.count)
  
  
  # I think we are officially preparing for plotting now...
if(plot.box){  

  # transposing mydat into long form
  mydat.long <- data.frame()
  for(i in colnames(mydat)){
    mydat.long <- rbind(mydat.long,data.frame(SUBJECT=rownames(mydat),VARIABLE=i,VALUE=mydat[,i]))
  }
  rownames(mydat.long) <- NULL
  
  # merging long mydat with outliers and creating label variable
  mydat.merged <- merge(mydat.long,outliers.formatted,all=T)
  mydat.merged$OUT.LABEL <- ifelse(is.na(mydat.merged$OUT.VALUE),"",as.character(mydat.merged$SUBJECT))
  
  
  # just going to print all variables.  Users will plot = F if too many variables!?

  par(mar=c(11.1,4.1,4.1,2.1))
  
  boxplot(VALUE ~ VARIABLE,data=mydat.merged,range=IQR,
          col="grey",pch=20, whisklty="solid", staplewex=0, # end of whisker is the "staple", las=2,
          main=sprintf("Univariate Boxplots with whiskers = %.2f x IQR",IQR),xaxt="n",
          ylab=ifelse(center,ifelse(scale,"Centered & Standardized Value", "Centered Value"),ifelse(scale,"Standardized Value","Value")))
  axis(side=1,at=1:ncol(mydat),labels=colnames(mydat),las=2,tck= -0.01,cex.axis=0.75)
  
  
  par(mar=c(5.1,4.1,4.1,2.1))

    # 
} # end of if(plot.box)
  return(outliers.out)
  }








