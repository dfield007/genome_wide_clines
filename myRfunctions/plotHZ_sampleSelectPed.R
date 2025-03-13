plotHZ_SampleSelectPed <- function(HZdata,xMin,xMax,yMin,yMax,gridAdd,adjustMin,TitleStart) {
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
  
  SampleSizePheno <- nrow(HZdata)
  HZdata_xy <- cbind(HZdata["Easting"], HZdata["Northing"])
  # TitleStart <- "test"
  MainTitle <- paste(TitleStart,paste("n = ",SampleSizePheno,sep=""),sep=", ")
  X11()
  plot(cbind(HZdata_xy[,"Easting"],HZdata_xy[,"Northing"]), 
       main = MainTitle, xlim = c(xMin,xMax), ylim = c(yMin,yMax), 
       xlab = "Easting (m)",ylab = "Northing (m)",pch = 21, cex=0.7,col="grey",
       bg="white",xaxs = "i", yaxs = "i")
  
  # core LR
  cat("\nCore, draw area...")
  selectedPoints_Core <- fhs(HZdata_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_Core <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_Core)),]
  points(selectedPoints_Core[,"Easting"],selectedPoints_Core[,"Northing"],col="orange")
  
  # yellow flanks 
  cat("\nYellow flanks , draw area...")
  selectedPoints_Yellow <- fhs(HZdata_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_Yellow <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_Yellow)),]
  points(selectedPoints_Yellow[,"Easting"],selectedPoints_Yellow[,"Northing"],col="yellow")
  # magenta flanks
  cat("\nMagenta flanks, draw area...")
  selectedPoints_Magenta <- fhs(HZdata_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_Magenta <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_Magenta)),]
  points(selectedPoints_Magenta[,"Easting"],selectedPoints_Magenta[,"Northing"],col="red")
  selectedPoints <- list()
  selectedPoints[["Core"]] <- selectedPoints_Core
  selectedPoints[["Yellow"]] <- selectedPoints_Yellow
  selectedPoints[["Magenta"]] <- selectedPoints_Magenta
  return(selectedPoints)
}
