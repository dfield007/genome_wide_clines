# Note 21/9/2019
# updated code to also return the x y coordinates of each polygon. 
# This could be reused in another function to set the sections of the hybrid zone for consistency (i.e. pedigree and tag assay)
plotHZ_SampleSelect <- function(HZdata,xMin,xMax,yMin,yMax,gridAdd,adjustMin,TitleStart) {
  # TitleStart <- "Planoles"
  # HZdata <- Antspec2017
  # xMin <- 421000; xMax <- 426000; yMin <- 4685000; yMax <- 4688000
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
  SampleSizePheno <- nrow(HZdata)
  HZdata_xy <- cbind(HZdata["Easting"], HZdata["Northing"])
  # TitleStart <- "test"
  MainTitle <- paste(TitleStart,paste("n = ",SampleSizePheno,sep=""),sep=", ")
  X11()
  plot(cbind(HZdata_xy[,"Easting"],HZdata_xy[,"Northing"]), 
       main = MainTitle, xlim = c(xMin,xMax), ylim = c(yMin,yMax), 
       xlab = "Easting (m)",ylab = "Northing (m)",pch = 21, cex=0.7,col="grey",
       bg="white",xaxs = "i", yaxs = "i")
 
  # setup list to receive data
  outputList <- list()
  outputList[["selectedPoints_Core"]] <- list()
  # core LR
  cat("\nCore, draw area...")
  selectedPoints_Core_raw <- fhs(HZdata_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_Core <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_Core_raw)),]
  outputList[["selectedPoints_Core"]] <- selectedPoints_Core
  points(selectedPoints_Core[,"Easting"],selectedPoints_Core[,"Northing"],col="orange")
  xyCoordinates_Core <- attributes(selectedPoints_Core_raw)
  xyCoordinates_Core <- xyCoordinates_Core$gate
  # yellow flanks 
  cat("\nYellow flanks , draw area...")
  outputList[["selectedPoints_Yellow"]] <- list()
  selectedPoints_Yellow_raw <- fhs(HZdata_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_Yellow <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_Yellow_raw)),]
  outputList[["selectedPoints_Yellow"]] <- selectedPoints_Yellow
  points(selectedPoints_Yellow[,"Easting"],selectedPoints_Yellow[,"Northing"],col="yellow")
  xyCoordinates_Yellow <- attributes(selectedPoints_Yellow_raw)
  xyCoordinates_Yellow <- xyCoordinates_Yellow$gate
  # magenta flanks
  cat("\nMagenta flanks, draw area...")
  outputList[["selectedPoints_Magenta"]] <- list()
  selectedPoints_Magenta_raw <- fhs(HZdata_xy)
  #selectedPoints <- sapply(list(HZdataCut[,"Easting"],HZdataCut[,"Northing"]),"[",identify(HZdataCut[,"Easting"],HZdataCut[,"Northing"]))
  selectedPoints_Magenta <- HZdata[rownames(HZdata)%in%as.numeric(as.vector(selectedPoints_Magenta_raw)),]
  outputList[["selectedPoints_Magenta"]] <- selectedPoints_Magenta
  points(selectedPoints_Magenta[,"Easting"],selectedPoints_Magenta[,"Northing"],col="red")
  xyCoordinates_Magenta <- attributes(selectedPoints_Magenta_raw)
  xyCoordinates_Magenta <- xyCoordinates_Magenta$gate
  # Sending polygon data to list
  outputList[["xyCoordinatesPolygon_Core"]] <- list()
  outputList[["xyCoordinatesPolygon_Yellow"]] <- list()
  outputList[["xyCoordinatesPolygon_Magenta"]] <- list()
  outputList[["xyCoordinatesPolygon_Core"]] <- xyCoordinates_Core
  outputList[["xyCoordinatesPolygon_Yellow"]] <- xyCoordinates_Yellow
  outputList[["xyCoordinatesPolygon_Magenta"]] <- xyCoordinates_Magenta
  return(outputList)
}
