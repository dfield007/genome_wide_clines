plotHZpedigree <- function(HZdata,xMin,xMax,yMin,yMax,gridAdd,adjustMin,TitleStart) {
  # TitleStart <- "Planoles"
  # HZdata <- Antspec2016
  # xMin <- 418500.0; xMax <- 427000.0; yMin <- 4685500.0; yMax <- 4687500.0
  # xMin <- 410335.1; xMax <- 434653; yMin <- 4683600; yMax <- 4691189
  # gridAdd <- TRUE
  HZdata[,"Easting"] <- as.vector(HZdata[,"Easting"])
  HZdata[,"Northing"] <- as.vector(HZdata[,"Northing"])
  HZdata <- HZdata[!is.na(HZdata[,"Easting"]),]
  HZdata <- HZdata[HZdata[,"Easting"]!="NA",]
  HZdata[,"Easting"] <- as.numeric(as.vector(HZdata[,"Easting"]))
  HZdata[,"Northing"] <- as.numeric(as.vector(HZdata[,"Northing"]))
  #HZdata <- HZdata[order(HZdata[,"Easting"]),]
  HZtoPlot_X <- as.vector(HZdata[,"Easting"])
  HZtoPlot_Y <- as.vector(HZdata[,"Northing"])
  HZtoPlot_XY <- cbind(HZtoPlot_X,HZtoPlot_Y)
  
  # adjustMin <- FALSE
  if (adjustMin==TRUE) {
    for (col in 1:ncol(HZtoPlot_XY)) {
      HZtoPlot_XY[,col] <- as.numeric(as.vector(HZtoPlot_XY[,col]))
    }
    HZtoPlot_XY <- data.table(HZtoPlot_XY)
    minAdjust1 <- min(as.numeric(HZtoPlot_XY[[1]]))
    minAdjust2 <- min(as.numeric(HZtoPlot_XY[[2]]))
    HZtoPlot_XY[[1]] <- as.numeric(HZtoPlot_XY[[1]])-minAdjust1
    HZtoPlot_XY[[2]] <- as.numeric(HZtoPlot_XY[[2]])-minAdjust2
    HZtoPlot_XY <- cbind(HZtoPlot_XY[[1]],HZtoPlot_XY[[2]])
  }
  HZdata[,"Easting"] <- HZtoPlot_XY[,1]
  HZdata[,"Northing"] <- HZtoPlot_XY[,2]
  #HZdata <- HZdata[!is.na(HZdata[,"phenoCat_final"]=="Y" | HZdata[,"phenoCat_final"]=="W" | HZdata[,"phenoCat_final"]=="FO" | HZdata[,"phenoCat_final"]=="WO" | HZdata[,"phenoCat_final"]=="FR" | HZdata[,"phenoCat_final"]=="WR"),]
  #   if (DataSet == "Clines_Planoles") {
  #     MainTitle <- paste(TitleStart,paste("n = ",SampleSizePheno,sep=""),sep=", ")
  #   }
  #   if (DataSet == "Clines_Cadi") {
  #     MainTitle <- paste("Cadi Hybrid zone",paste("n = ",SampleSizePheno,sep=""),sep=", ")
  #   }
  #   if (DataSet == "Planoles") {
  #     MainTitle <- paste("Planoles Hybrid zone",paste("n = ",SampleSizePheno,sep=""),sep=", ")
  #   }
  #   if (DataSet == "Cadi") {
  #     MainTitle <- paste("Cadi Hybrid zone",paste("n = ",SampleSizePheno,sep=""),sep=", ")
  #   }
  HZdataCut <- HZdata[HZdata[,"Easting"]>=xMin & HZdata[,"Easting"]<=xMax & HZdata[,"Northing"]>=yMin & HZdata[,"Northing"]<=yMax,]
  SampleSizePheno <- nrow(HZdataCut)
  # TitleStart <- "test"
  MainTitle <- paste(TitleStart,paste("n = ",SampleSizePheno,sep=""),sep=", ")
  plot(cbind(HZdataCut[,"Easting"],HZdataCut[,"Northing"]), 
       main = MainTitle, xlim = c(xMin,xMax), ylim = c(yMin,yMax), 
       xlab = "Easting (m)",ylab = "Northing (m)",col = "black", bg = "black",type="n",xaxs = "i", yaxs = "i")
  # gridAdd <- T
  if (gridAdd == TRUE) {
    XvaluesLines1 <- sort(rep(seq(xMin,xMax,200),2))
    YvalueLines1 <- rep(c(yMin,yMax),(length(XvaluesLines1)/2))
    Xcords1 <- matrix(XvaluesLines1,(length(XvaluesLines1)/2),2,byrow=TRUE)
    Ycords1 <- matrix(YvalueLines1,(length(YvalueLines1)/2),2,byrow=TRUE)
    
    YvaluesLines2 <- sort(rep(seq(yMin,yMax,200),2))
    XvalueLines2 <- rep(c(xMin,xMax),(length(YvaluesLines2)/2))
    Xcords2 <- matrix(XvalueLines2,(length(XvalueLines2)/2),2,byrow=TRUE)
    Ycords2 <- matrix(YvaluesLines2,(length(YvaluesLines2)/2),2,byrow=TRUE)
    for (thisRow in 1:nrow(Xcords1)) {
      # thisRow <- 1
      lines(Xcords1[thisRow,], Ycords1[thisRow,], type="l", pch=22, col="darkgrey",lty=1,lwd=0.4)
    }
    for (thisRow in 1:nrow(Xcords2)) {
      lines(Xcords2[thisRow,], Ycords2[thisRow,], type="l", pch=22, col="darkgrey",lty=1,lwd=0.4)
    }
  }
  
    # Note: if transparent points are preffered, change the 4th number in rgb tp <100
    points(cbind(as.vector(HZdataCut[HZdataCut[,"phenoCat_final"]==-9 | HZdataCut[,"phenoCat_final"]==0,"Easting"]),
                 as.vector(HZdataCut[HZdataCut[,"phenoCat_final"]==-9 | HZdataCut[,"phenoCat_final"]==0,"Northing"])), 
           pch = 21, cex=0.7,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(255,255,255,255,maxColorValue=255))
    points(cbind(as.vector(HZdataCut[HZdataCut[,"phenoCat_final"]=="Y","Easting"]),
                 as.vector(HZdataCut[HZdataCut[,"phenoCat_final"]=="Y","Northing"])), 
           pch = 21, cex=0.7,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(255,255,0,255,maxColorValue=255))
    points(cbind(as.vector(HZdataCut[HZdataCut[,"phenoCat_final"]=="FR","Easting"]),
                 as.vector(HZdataCut[HZdataCut[,"phenoCat_final"]=="FR","Northing"])), 
           pch = 21, cex=0.7,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(255,0,0,255,maxColorValue=255))
    points(cbind(as.vector(HZdataCut[HZdataCut[,"phenoCat_final"]=="FO","Easting"]),
                 as.vector(HZdataCut[HZdataCut[,"phenoCat_final"]=="FO","Northing"])), 
           pch = 21, cex=0.7,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(255,165,0,255,maxColorValue=255))
    points(cbind(as.vector(HZdataCut[HZdataCut[,"phenoCat_final"]=="WR","Easting"]),
                 as.vector(HZdataCut[HZdataCut[,"phenoCat_final"]=="WR","Northing"])), 
           pch = 21, cex=0.7,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(255,195,203,255,maxColorValue=255))
    points(cbind(as.vector(HZdataCut[HZdataCut[,"phenoCat_final"]=="W","Easting"]),
                 as.vector(HZdataCut[HZdataCut[,"phenoCat_final"]=="W","Northing"])), 
           pch = 21, cex=0.7,col="black",bg=rgb(255,255,255,255,maxColorValue=255))
    points(cbind(as.vector(HZdataCut[HZdataCut[,"phenoCat_final"]=="WO","Easting"]),
                 as.vector(HZdataCut[HZdataCut[,"phenoCat_final"]=="WO","Northing"])), 
           pch = 21, cex=0.7,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(238,216,174,255,maxColorValue=255))

}
