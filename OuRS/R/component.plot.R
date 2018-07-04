component.plot <- function(scores, axes=c(1,2), pch=20, col="mediumorchid4", line.col="grey80", lty=2, lwd=2,
                              main="Component scores",
                              xlab=paste0("Component ",axes[1]),
                              ylab=paste0("Component ",axes[2]),
                              xlim=c(-max(abs(scores[,axes])),max(abs(scores[,axes])))*1.2,
                              ylim=c(-max(abs(scores[,axes])),max(abs(scores[,axes])))*1.2,
                              asp=1, pos=3, display_names=T,cex=1,text.cex=1,
                              ...){

  plot(0, type="n", xlim=xlim, ylim=ylim, main=main, xlab=xlab, ylab=ylab, axes=F, asp=asp)
  abline(h=0,lty=2,lwd=2, col="grey60")
  abline(v=0,lty=2,lwd=2, col="grey60")
  points(scores[,axes], col=col, pch=pch, cex=cex, ...)
    ## will try to employ a "repel" later.
  
  ### actually, display names should be a vector
  # if(length(display_names)==1){
  #  display_names <- rep(display_names,nrow(scores))
  # }else if(length(display_names)!=nrow(scores)){
  #  ## a default.
  #  warning("length(display_names) is not 1 or nrow(scores). Setting to default.")
  #  display_names <- rep(T,nrow(scores))
  # }
    
  # if( !all(unlist(lapply(display_names,is.logical))) | any(is.na(display_names)) ){
  #  warning("Not all items in 'display_names' were TRUE or FALSE")
  #  display_names <- rep(T,nrow(scores))
  # }
  # test if display_names is all logical else set to length of vector as false
  
    ## this might get weird if there is only 1.
  #if (sum(display_names) > 0) {
  if (display_names) {
    text(scores[, axes], labels = rownames(scores), 
         pos = pos, col = col, cex = text.cex)
  }

}
