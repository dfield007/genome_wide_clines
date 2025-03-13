plotTerrainPiePlot <- function(DataSetName,DEM_merged,SNPsFullData,DemeData,locusToPlot,
                               SNPlociListUpdated,minNum,xRange,yRange,
                               alphaTerm,colPallete,colPalleteNum,poolsAdd,
                               demeSize,plotPies,pieRadius,pointsInd,ptSize,gridAdd,landMarks,fileFormat,
                               lowerRoadData,upperRoadData,markRoads=T) {
  # SNPlociListUpdated <- summaryTableHZ
  # xRange=c(410000,435572.1); yRange=c(4684000,4691800)
  # xRange=c(418000,429572.1); yRange=c(4684000,4689800)
  # range for Cadi
  # xRange=c(404550,408800); yRange=c(4677200,4684600)
  # DemeData <- DemeData_100m_Cadi
  # DemeData <- HZ_Planoles_200mDemes
  
  # names(DataSet)
  # locusToPlot <- "ros_assembly_543443"
  # alphaTerm <- 0.8
  # colPallete <- "land.pal"
  # colPalleteNum <- 50
  # SNPsFullData <- HZ_Planoles_Clines
  # pieRadius <- 120
  # demeSize <- 200
  # pointsInd <- "none"
  # gridAdd <- FALSE
  # landMarks <- FALSE
  # fileFormat <- "pdf"
  # thislocusName <- "ros_assembly_543443"
  # DataSet <- "Clines_Planoles"
  # minNum <- 10
  # markRoads <- T
  land.pal <- colorRampPalette(
    c("#336600", "#F3CA89", "#D9A627", 
      "#A49019", "#9F7B0D", "#996600", "#B27676", "#C2B0B0", "#E5E5E5", 
      "#FFFFFF"))
  thisLocusName <- locusToPlot
  if (fileFormat=="eps") {
    filenamePiePlot <- paste(thisLocusName,"_PiePlot.eps")
  }
  if (fileFormat=="pdf") {
    filenamePiePlot <- paste(thisLocusName,"_PiePlot.pdf")
  }
  if (DataSetName=="Clines_Planoles") {
    if (fileFormat=="eps") {
      cairo_ps(filenamePiePlot,width = 24, height = 9.8, pointsize = 14,antialias="default")
      par(mar=c(5.5,4.5,4,2) + 0.1) 
    }
    if (fileFormat=="pdf") {
      pdf(filenamePiePlot,width = 24, height = 9.8, pointsize = 12)
      par(mfrow=c(1,1))
      par(mar=c(5.5,4.5,4,2) + 0.1) 
      #par(mar=c(2,2,2,2) + 0.1) 
      
    }
  }
  if (DataSetName=="Clines_Cadi") {
    if (fileFormat=="eps") {
      cairo_ps(filenamePiePlot,width = 12, height = 16, pointsize = 14,antialias="default")
      par(mar=c(5.5,4.5,4,2) + 0.1)
    }
    if (fileFormat=="pdf") {
      pdf(filenamePiePlot,width = 12, height = 16, pointsize = 14)
      par(mar=c(5.5,4.5,4,2) + 0.1)
      #par(mar=c(1,1,1,1) + 0.1)
      
    }
  }
  markerLabelPieChart <- function(thisLocusName,SNPlociListUpdated) {
    # SampleSize <- 1514
    print("debug1")
    ThisLGDataSetHZ <- summaryTableHZ
    markerLabelMake <- NULL
    SampleLabel <- NULL
    colnames(SNPlociListUpdated)
    linkageGroup <- SNPlociListUpdated[SNPlociListUpdated[,"LocusName"]==thisLocusName,"LG"]
    cM <- SNPlociListUpdated[SNPlociListUpdated[,"LocusName"]==thisLocusName,"cM"]
    Note <- as.vector(SNPlociListUpdated[SNPlociListUpdated[,"LocusName"]==thisLocusName,"Type"])
    fullLabel <- paste(thisLocusName,
                       paste(paste("LG",linkageGroup,sep=" "),
                             paste(paste("cM",cM,sep=" "),Note,sep=", "),sep=", "),sep=", ")
    return(fullLabel)
  }
  chartLabel <- locusToPlot
  par(mfrow=c(1,1))
  par(mar=c(2,1,2,1))
  par(oma=c(2,1,2,1))
  plot(DEM_merged,main=chartLabel,ylim=yRange,xlim=xRange,alpha=alphaTerm,
       axes=FALSE,yaxt="n",xaxt="n",ann=F,cex=2,
       col=land.pal(colPalleteNum))
  thisLocusData <- DemeData[[thisLocusName]]
  # Get allele frequency data
  thisDemeMatrix2D_pCount <- thisLocusData[["DemeMatrix2D_pCount"]]
  thisDemeMatrix2D_pqCountNoZero <- thisLocusData[["DemeMatrix2D_pqCountNoZero"]]
  thisDemeMatrix2D_pqCountNoZero <- thisDemeMatrix2D_pqCountNoZero[thisDemeMatrix2D_pqCountNoZero[,"p"]!="NaN",]
  sum(thisDemeMatrix2D_pqCountNoZero[,"pCount"])+ sum(thisDemeMatrix2D_pqCountNoZero[,"qCount"])
  # mark roads?
  if (markRoads==T) {
    plot(lowerRoad_UTM,add=T,col='black',lwd=2)
    plot(upperRoad_UTM,add=T,col='black',lwd=2)
    #plot(diagonal_UTM,add=T,col='black',lwd=2)
  }
  if (pointsInd=="gray")  {
    SNPsFullDataPheno <- SNPsFullData[SNPsFullData[,"phenoCat"]!=-9,]
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="Y","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="Y","Northing"])), 
           pch = 21, cex=ptSize,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(255,255,0,255,maxColorValue=255))
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="FR","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="FR","Northing"])), 
           pch = 21, cex=ptSize,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(255,0,0,255,maxColorValue=255))
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="FO","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="FO","Northing"])), 
           pch = 21, cex=ptSize,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(255,165,0,255,maxColorValue=255))
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="WR","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="WR","Northing"])), 
           pch = 21, cex=ptSize,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(255,195,203,255,maxColorValue=255))
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="W","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="W","DistNorthOfCentre"])), 
           pch = 21, cex=ptSize,col="black",bg=rgb(255,255,255,255,maxColorValue=255))
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="WO","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="WO","Northing"])), 
           pch = 21, cex=ptSize,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(238,216,174,255,maxColorValue=255))
  }
  if (pointsInd=="colour") {
    #write.csv(SNPsFullData,"SNPsFullData.csv")
    #SNPsFullData
    SNPsFullDataPheno <- SNPsFullData[SNPsFullData[,"phenoCat"]!=-9,]
    #SNPsFullDataPheno[SNPsFullDataPheno[,"ID"]=="M3001","phenoCat"]
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="Y","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="Y","Northing"])), 
           pch = 21, cex=ptSize,col=colors()[312],bg=colors()[312])
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="FR","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="FR","Northing"])), 
           pch = 21, cex=ptSize,col=colors()[312],bg=colors()[312])
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="FO","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="FO","Northing"])), 
           pch = 21, cex=ptSize,col=colors()[312],bg=colors()[312])
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="WR","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="WR","Northing"])), 
           pch = 21, cex=ptSize,col=colors()[312],bg=colors()[312])
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="W","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="W","Northing"])), 
           pch = 21, cex=ptSize,col=colors()[312],bg=colors()[312])
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="WO","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat"]=="WO","Northing"])), 
           pch = 21, cex=ptSize,col=colors()[312],bg=colors()[312])
  }  
  
  # Allele freqs - Change polarity where necessary
  # last check for NaN
  thisDemeMatrix2D_pqCountNoZero <- thisDemeMatrix2D_pqCountNoZero[thisDemeMatrix2D_pqCountNoZero[,"p"]!="NaN",]
  if (mean(thisDemeMatrix2D_pqCountNoZero[1:4,"p"])>0.60 & mean(thisDemeMatrix2D_pqCountNoZero[(nrow(thisDemeMatrix2D_pqCountNoZero)-5):(nrow(thisDemeMatrix2D_pqCountNoZero)-1),"p"]) < 0.40) {
    #cat(c("..switching polarity"))
    thisDemeMatrix2D_pqCountNoZero_temp <- thisDemeMatrix2D_pqCountNoZero
    thisDemeMatrix2D_pqCountNoZero_temp[,"q"] <- thisDemeMatrix2D_pqCountNoZero[,"p"]
    thisDemeMatrix2D_pqCountNoZero_temp[,"p"] <- thisDemeMatrix2D_pqCountNoZero[,"q"]
    thisDemeMatrix2D_pqCountNoZero_temp[,"qCount"] <- thisDemeMatrix2D_pqCountNoZero[,"pCount"]
    thisDemeMatrix2D_pqCountNoZero_temp[,"pCount"] <- thisDemeMatrix2D_pqCountNoZero[,"qCount"]
    #thisDataLoc_good <- thisDataLoc
    thisDemeMatrix2D_pqCountNoZero <- thisDemeMatrix2D_pqCountNoZero_temp
  }
  # a few stragglers not picked up
  if (thisLocusName=="s678C8923637_374761" | thisLocusName=="s678C8923637_490488" | thisLocusName=="s117reverse_18263" | thisLocusName=="s117reverse_131947") {
    #cat(c("..switching polarity"))
    thisDemeMatrix2D_pqCountNoZero_temp <- thisDemeMatrix2D_pqCountNoZero
    thisDemeMatrix2D_pqCountNoZero_temp[,"q"] <- thisDemeMatrix2D_pqCountNoZero[,"p"]
    thisDemeMatrix2D_pqCountNoZero_temp[,"p"] <- thisDemeMatrix2D_pqCountNoZero[,"q"]
    thisDemeMatrix2D_pqCountNoZero_temp[,"qCount"] <- thisDemeMatrix2D_pqCountNoZero[,"pCount"]
    thisDemeMatrix2D_pqCountNoZero_temp[,"pCount"] <- thisDemeMatrix2D_pqCountNoZero[,"qCount"]
    #thisDataLoc_good <- thisDataLoc
    thisDemeMatrix2D_pqCountNoZero <- thisDemeMatrix2D_pqCountNoZero_temp
  }
  rownames(thisDemeMatrix2D_pqCountNoZero) <- 1:nrow(thisDemeMatrix2D_pqCountNoZero)
  
  # find axes in metres
  minX <- min(thisDemeMatrix2D_pqCountNoZero[,"X"])
  XinMetres <- thisDemeMatrix2D_pqCountNoZero[,"X"]-minX
  minY <- min(thisDemeMatrix2D_pqCountNoZero[,"Y"])
  YinMetres <- sort(thisDemeMatrix2D_pqCountNoZero[,"Y"]-minY)
  max(thisDemeMatrix2D_pqCountNoZero[,"Y"])
  # metre range within plotting range
  Xkeep <- thisDemeMatrix2D_pqCountNoZero[,"X"]>xRange[1] & thisDemeMatrix2D_pqCountNoZero[,"X"]<xRange[2]
  Ykeep <- thisDemeMatrix2D_pqCountNoZero[,"Y"]>yRange[1] & thisDemeMatrix2D_pqCountNoZero[,"Y"]<yRange[2]
  Easting <- thisDemeMatrix2D_pqCountNoZero[,"X"]
  Northing <- thisDemeMatrix2D_pqCountNoZero[,"Y"]
  Easting_keep <- sort(Easting[Xkeep & Ykeep])
  Northing_keep <- sort(Northing[Ykeep])
  
  # XinMetres_keep <- sort(XinMetres[Xkeep & Ykeep])
  # YinMetres_keep <- sort(YinMetres[Ykeep & Ykeep])
  # labelPos_ticks <- seq(minY,max(Easting_keep),1000)
  # 
  # labelPos_Xeasting <- seq(min(Easting_keep),max(Easting_keep),1000)
  # labelPos_Xmetres <- seq(min(XinMetres_keep),max(XinMetres_keep),1000)
  # labelPos_Xeasting <- seq(min(Easting_keep),max(Easting_keep),1000)
  # labelPos_Xmetres <- seq(min(XinMetres_keep),max(XinMetres_keep),1000)
  # labelPos_Xeasting <- seq(min(Easting_keep),max(Easting_keep),1000)
  # labelPos_Xmetres <- seq(min(XinMetres_keep),max(XinMetres_keep),1000)
  
  labelPos_Xeasting <- seq(min(Easting),max(Easting),1000)
  labelPos_Xmetres <- seq(min(XinMetres),max(XinMetres),1000)
  labelPos_Ynorthing <- seq(min(Northing),max(Northing)+600,1000)
  labelPos_Ymetres <- seq(min(YinMetres),max(YinMetres)+600,1000)
  #rect(xleft=-0.2, ybottom=39, xright=4.2, ytop=41, col='white',density = NULL, border = NA)
  #rect(xleft=-0.2, ybottom=44, xright=4.2, ytop=45, col='white',density = NULL, border = NA)
  #axis(2, at=c(41,41.5,42,42.5,43,43.5,44), labels=NA,col.axis="black", las=1,cex.axis=1.5,lwd=2.7,pos=0, tck = -0.01)
  #axis(2, at=c(41.5,42.5,43.5), labels=NA,col.axis="black", las=1,cex.axis=1.5,lwd=2.7,pos=0, tck = -0.005)
  #yLabels <-c(expression(paste(41^o,"N",sep=" ")),expression(paste(42^o,"N",sep=" ")),expression(paste(43^o,"N",sep=" ")),
  #            expression(paste(44^o,"N",sep=" ")))
  #xLabels <-c(expression(paste(0^o,"E",sep=" ")),expression(paste(1^o,"E",sep=" ")),expression(paste(2^o,"E",sep=" ")),
  #            expression(paste(3^o,"E",sep=" ")),expression(paste(4^o,"E",sep=" ")))
  #axis(2, at=c(41,42,43,44), labels=yLabels,col.axis="black", las=1,cex.axis=1.7,lwd=0,pos=0.09, tck = 0)
  axis(1, at=labelPos_Xeasting, labels=NA,col.axis="black", las=1,cex.axis=1.5,lwd=2.7,pos=min(Northing_keep), tck = -0.01)
  axis(1, at=labelPos_Xeasting, labels=labelPos_Xmetres/1000,col.axis="black", las=1,cex.axis=1.8,lwd=0,pos=min(Northing_keep), tck = -0.01)
  
  if (DataSetName=="Clines_Cadi") {
    axis(2, at=labelPos_Ynorthing, labels=NA,col.axis="black", las=1,cex.axis=1.5,lwd=2.7,pos=min(Easting_keep), tck = -0.01)
    axis(2, at=labelPos_Ynorthing, labels=labelPos_Ymetres/1000,col.axis="black", las=1,cex.axis=1.8,lwd=0,pos=min(Easting_keep), tck = -0.01)
  }
  # plotting pie charts
  totalN <- (thisDemeMatrix2D_pqCountNoZero[,"pCount"]+thisDemeMatrix2D_pqCountNoZero[,"qCount"])/2 
  keep <- which(totalN>=minNum)
  thisDemeMatrix2D_pqCountNoZero <- thisDemeMatrix2D_pqCountNoZero[keep,]
  
  
  for (thisDeme2D in 1:nrow(thisDemeMatrix2D_pqCountNoZero)) {
    # thisDeme2D <- 1
    thisXraw <- as.numeric(as.vector(thisDemeMatrix2D_pqCountNoZero[thisDeme2D,"X"]))
    thisYraw <- as.numeric(as.vector(thisDemeMatrix2D_pqCountNoZero[thisDeme2D,"Y"]))
    #     # split X and Y
    #     thisXrange <- as.numeric(unlist(strsplit(thisXraw,"_")))
    #     thisXmid <- min(thisXrange)+((max(thisXrange)-min(thisXrange))/2)
    #     thisYrange <- as.numeric(unlist(strsplit(thisYraw,"_")))
    #     thisYmid <- min(thisYrange)+((max(thisYrange)-min(thisYrange))/2)
    thisP <- as.numeric(as.vector(thisDemeMatrix2D_pqCountNoZero[thisDeme2D,"pCount"]))
    thisQ <- as.numeric(as.vector(thisDemeMatrix2D_pqCountNoZero[thisDeme2D,"qCount"]))
    cat("\n .",thisDeme2D)
    cat("..P/Q:",c(thisP,thisQ))
    #add.pie(z=c(thisP,thisQ), x=thisX, y=thisY, radius=40,col=c("red", "yellow"),labels="")
    # pie plot size has to be tweaked depending on deme size
    # small pie plots, original was radius 69 for deme size 200
    #     if (sum(thisP,thisQ)>minNum) {
    #       add.pie(z=c(thisP,thisQ), x=thisXmid, y=thisYmid, radius=69,col=c((rgb(255,0,0,255,maxColorValue=255)),(rgb(255,255,0,255,maxColorValue=255))),labels="")
    #     }
    # larger pie plots - 120 works okay for deme size 400
    totalN <- thisP + thisQ
    #if ((totalN < minNum)==F) {
    add.pie(z=c(thisP,thisQ), x=thisXraw, y=thisYraw, radius=pieRadius,col=c((rgb(255,0,0,255,maxColorValue=255)),(rgb(255,255,0,255,maxColorValue=255))),labels="")
    #}
  }
  
  # mark pools?
  if (poolsAdd==TRUE) {
    points(x=411604.6046997,y=4690386.0014007,pch=21,bg='blue',col='black') # YP4 - pool 1
    points(x=421903.0320033,y=4686484.3909220,pch=21,bg='blue',col='black') # YP - pool 3
    points(x=422119.7797325,y=4686390.2520109,pch=21,bg='blue',col='black') # YP - pool 4
    points(x=424462.3465657,y=4686112.3149515,pch=21,bg='blue',col='black') # MP - pool 5
    points(x=425146.6755758,y=4686021.5773514,pch=21,bg='blue',col='black') # MP - pool 6
    points(x=431621.6848923,y=4686847.1448457,pch=21,bg='blue',col='black') # MP - pool 6
  }
  dev.off()  
}
