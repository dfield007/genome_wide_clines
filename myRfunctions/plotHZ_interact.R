plotHZ_designateRegions <- function(HZdata,xMin,xMax,yMin,yMax,gridAdd,adjustMin,TitleStart) {
  # TitleStart <- "Planoles"
  # HZdata <- Antspec2016
  # xMin <- 418500.0; xMax <- 427000.0; yMin <- 4685500.0; yMax <- 4687500.0
  # xMin <- 410335.1; xMax <- 434653; yMin <- 4683600; yMax <- 4691189
  # gridAdd <- TRUE
  HZdata[,"Easting"] <- as.numeric(HZdata[,"Easting"])
  HZdata[,"Northing"] <- as.numeric(HZdata[,"Northing"])
  #HZdata <- HZdata[order(HZdata[,"Easting"]),]
  HZdata <- HZdata[!is.na(HZdata[,"Easting"]),]
  HZtoPlot_X <- as.vector(HZdata[,"Easting"])
  HZtoPlot_X <- HZtoPlot_X[!is.na(HZtoPlot_X)]
  HZtoPlot_X <- as.numeric(HZtoPlot_X)
  HZtoPlot_Y <- as.vector(HZdata[,"Northing"])
  HZtoPlot_Y <- HZtoPlot_Y[!is.na(HZtoPlot_Y)]
  HZtoPlot_Y <- as.numeric(HZtoPlot_Y)
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
  #HZdata <- HZdata[!is.na(HZdata[,"PhenoCat"]=="Y" | HZdata[,"PhenoCat"]=="W" | HZdata[,"PhenoCat"]=="FO" | HZdata[,"PhenoCat"]=="WO" | HZdata[,"PhenoCat"]=="FR" | HZdata[,"PhenoCat"]=="WR"),]
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
  #HZdataCut <- HZdata[HZdata[,"Easting"]>=xMin & HZdata[,"Easting"]<=xMax & HZdata[,"Northing"]>=yMin & HZdata[,"Northing"]<=yMax,]
  SampleSizePheno <- nrow(HZdataCut)
  HZdata_xy <- cbind(HZdata["Easting"], HZdata["Northing"])
  # TitleStart <- "test"
  MainTitle <- paste(TitleStart,paste("n = ",SampleSizePheno,sep=""),sep=", ")
  X11()
  plot(cbind(HZdata_xy[,"Easting"],HZdata_xy[,"Northing"]), 
       main = MainTitle, xlim = c(xMin,xMax), ylim = c(yMin,yMax), 
       xlab = "Easting (m)",ylab = "Northing (m)",pch = 21, cex=0.7,col="grey",
       bg="white",xaxs = "i", yaxs = "i")
  # 
  # # gridAdd <- T
  # if (gridAdd == TRUE) {
  #   XvaluesLines1 <- sort(rep(seq(xMin,xMax,200),2))
  #   YvalueLines1 <- rep(c(yMin,yMax),(length(XvaluesLines1)/2))
  #   Xcords1 <- matrix(XvaluesLines1,(length(XvaluesLines1)/2),2,byrow=TRUE)
  #   Ycords1 <- matrix(YvalueLines1,(length(YvalueLines1)/2),2,byrow=TRUE)
  #   
  #   YvaluesLines2 <- sort(rep(seq(yMin,yMax,200),2))
  #   XvalueLines2 <- rep(c(xMin,xMax),(length(YvaluesLines2)/2))
  #   Xcords2 <- matrix(XvalueLines2,(length(XvalueLines2)/2),2,byrow=TRUE)
  #   Ycords2 <- matrix(YvaluesLines2,(length(YvaluesLines2)/2),2,byrow=TRUE)
  #   for (thisRow in 1:nrow(Xcords1)) {
  #     # thisRow <- 1
  #     lines(Xcords1[thisRow,], Ycords1[thisRow,], type="l", pch=22, col="darkgrey",lty=1,lwd=0.4)
  #   }
  #   for (thisRow in 1:nrow(Xcords2)) {
  #     lines(Xcords2[thisRow,], Ycords2[thisRow,], type="l", pch=22, col="darkgrey",lty=1,lwd=0.4)
  #   }
  # }
  # # Note: if transparent points are preffered, change the 4th number in rgb tp <100
  # points(cbind(as.vector(HZdataCut[HZdataCut[,"phenoCat"]==-9 | HZdataCut[,"phenoCat"]==0,"Easting"]),
  #              as.vector(HZdataCut[HZdataCut[,"phenoCat"]==-9 | HZdataCut[,"phenoCat"]==0,"Northing"])), 
  #        pch = 21, cex=0.7,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(255,255,255,255,maxColorValue=255))
  # points(cbind(as.vector(HZdataCut[HZdataCut[,"phenoCat"]=="Y","Easting"]),
  #              as.vector(HZdataCut[HZdataCut[,"phenoCat"]=="Y","Northing"])), 
  #        pch = 21, cex=0.7,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(255,255,0,255,maxColorValue=255))
  # points(cbind(as.vector(HZdataCut[HZdataCut[,"phenoCat"]=="FR","Easting"]),
  #              as.vector(HZdataCut[HZdataCut[,"phenoCat"]=="FR","Northing"])), 
  #        pch = 21, cex=0.7,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(255,0,0,255,maxColorValue=255))
  # points(cbind(as.vector(HZdataCut[HZdataCut[,"phenoCat"]=="FO","Easting"]),
  #              as.vector(HZdataCut[HZdataCut[,"phenoCat"]=="FO","Northing"])), 
  #        pch = 21, cex=0.7,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(255,165,0,255,maxColorValue=255))
  # points(cbind(as.vector(HZdataCut[HZdataCut[,"phenoCat"]=="WR","Easting"]),
  #              as.vector(HZdataCut[HZdataCut[,"phenoCat"]=="WR","Northing"])), 
  #        pch = 21, cex=0.7,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(255,195,203,255,maxColorValue=255))
  # points(cbind(as.vector(HZdataCut[HZdataCut[,"phenoCat"]=="W","Easting"]),
  #              as.vector(HZdataCut[HZdataCut[,"phenoCat"]=="W","Northing"])), 
  #        pch = 21, cex=0.7,col="black",bg=rgb(255,255,255,255,maxColorValue=255))
  # points(cbind(as.vector(HZdataCut[HZdataCut[,"phenoCat"]=="WO","Easting"]),
  #              as.vector(HZdataCut[HZdataCut[,"phenoCat"]=="WO","Northing"])), 
  #        pch = 21, cex=0.7,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(238,216,174,255,maxColorValue=255))
  # core LR
  cat("\nCore Lower Road, draw area...")
  selectedPoints_LR <- fhs(HZdataCut_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_LR <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_LR)),]
  points(selectedPoints_LR[,"Easting"],selectedPoints_LR[,"Northing"],col="orange")
  # core UpperRoad
  cat("\nCore Upper Road, draw area...")
  selectedPoints_UR <- fhs(HZdataCut_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_UR <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_UR)),]
  points(selectedPoints_UR[,"Easting"],selectedPoints_UR[,"Northing"],col="orange")
  
  # yellow flank 1 LowerRoad
  cat("\nYellow flank 1 Lower Road, draw area...")
  selectedPoints_YF1_LR <- fhs(HZdataCut_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_YF1_LR <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_YF1_LR)),]
  points(selectedPoints_YF1_LR[,"Easting"],selectedPoints_YF1_LR[,"Northing"],col="yellow")
  # yellow flank 2 
  cat("\nYellow flank 2 Lower Road, draw area...")
  selectedPoints_YF2_LR <- fhs(HZdataCut_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_YF2_LR <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_YF2_LR)),]
  points(selectedPoints_YF2_LR[,"Easting"],selectedPoints_YF2_LR[,"Northing"],col="wheat")
  
  # yellow flank 1 UpperRoad
  cat("\nYellow flank 1 Upper Road, draw area...")
  selectedPoints_YF1_UR <- fhs(HZdataCut_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_YF1_UR <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_YF1_UR)),]
  points(selectedPoints_YF1_UR[,"Easting"],selectedPoints_YF1_UR[,"Northing"],col="yellow")
  # yellow flank 2 UpperRoad
  cat("\nYellow flank 2 Upper Road, draw area...")
  selectedPoints_YF2_UR <- fhs(HZdataCut_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_YF2_UR <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_YF2_UR)),]
  points(selectedPoints_YF2_UR[,"Easting"],selectedPoints_YF2_UR[,"Northing"],col="wheat")
  
  # yellow diagonal road
  cat("\nYellow Diagonal Road, draw area...")
  selectedPoints_YF_Diag <- fhs(HZdataCut_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_YF_Diag <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_YF_Diag)),]
  points(selectedPoints_YF_Diag[,"Easting"],selectedPoints_YF_Diag[,"Northing"],col="yellow")
  
  # magenta flank 1 LowerRoad
  cat("\nMagenta flank 1 Lower Road, draw area...")
  selectedPoints_MF1_LR <- fhs(HZdataCut_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_MF1_LR <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_MF1_LR)),]
  points(selectedPoints_MF1_LR[,"Easting"],selectedPoints_MF1_LR[,"Northing"],col="red")
  # magenta flank 2 LowerRoad
  cat("\nMagenta flank 2 Lower Road, draw area...")
  selectedPoints_MF2_LR <- fhs(HZdataCut_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_MF2_LR <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_MF2_LR)),]
  points(selectedPoints_MF2_LR[,"Easting"],selectedPoints_MF2_LR[,"Northing"],col="magenta")
  # magenta flank 3 LowerRoad
  cat("\nMagenta flank 3 Lower Road, draw area...")
  selectedPoints_MF3_LR <- fhs(HZdataCut_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_MF3_LR <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_MF3_LR)),]
  points(selectedPoints_MF3_LR[,"Easting"],selectedPoints_MF3_LR[,"Northing"],col="pink")
  
  # magenta flank 1 upperRoad
  cat("\nMagenta flank 1 Upper Road, draw area...")
  selectedPoints_MF1_UR <- fhs(HZdataCut_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_MF1_UR <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_MF1_UR)),]
  points(selectedPoints_MF1_UR[,"Easting"],selectedPoints_MF1_UR[,"Northing"],col="red")
  # magenta flank 2 LowerRoad
  cat("\nMagenta flank 2 Upper Road, draw area...")
  selectedPoints_MF2_UR <- fhs(HZdataCut_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_MF2_UR <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_MF2_UR)),]
  points(selectedPoints_MF2_LR[,"Easting"],selectedPoints_MF2_LR[,"Northing"],col="magenta")
  # magenta flank 3 UpperRoad
  cat("\nMagenta flank 3 Upper Road, draw area...")
  selectedPoints_MF3_UR <- fhs(HZdataCut_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_MF3_UR <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_MF3_UR)),]
  points(selectedPoints_MF3_UR[,"Easting"],selectedPoints_MF3_UR[,"Northing"],col="pink")
  return(list(selectedPoints_LR,selectedPoints_UR,selectedPoints_YF1_LR,selectedPoints_YF2_LR,
              selectedPoints_YF1_UR,selectedPoints_YF2_UR,selectedPoints_YF_Diag,
              selectedPoints_MF1_LR,selectedPoints_MF2_LR,selectedPoints_MF3_LR,
              selectedPoints_MF1_UR,selectedPoints_MF2_UR,selectedPoints_MF3_UR))
}
