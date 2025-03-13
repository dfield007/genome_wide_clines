# DemeSetup2D_LD
# Author: David. L. Field
# 14/12/2019

# This is a function makes a 2D terrain (altitude) plot with pie charts for allele frequencies on top.

plotTerrainPiePlot <- function(DataSetName,DEM_merged,SNPsFullData,DemeData,locusToPlot,SNPlociListUpdated,minNum,xRange,
                               yRange,alphaTerm,colPallete,colPalleteNum,poolsAdd,demeSize,plotPies,pieRadius,
                               pointsInd,ptSize,gridAdd,landMarks,fileFormat,markRoads=T) {
  # DataSetName <- "Planoles"
  # SNPlociListUpdated <- SNPlociListUpdated
  # xRange=c(410000,435572.1); yRange=c(4684000,4691800)
  # xRange=c(418000,429572.1); yRange=c(4684000,4689800)
  # range for Cadi
  # xRange=c(404550,408800); yRange=c(4677200,4684600)
  # DemeData <- HZ_Planoles_100mDemes
  # SNPsFullData <- GenoEcol_cut
  # names(DataSet)
  # locusToPlot <- "ros_assembly_543443"
  # alphaTerm <- 0.8; ptSize <- 0.7; colPallete <- land.pal; colPalleteNum <- 50; pieRadius <- 80; demeSize <- 100
  # pointsInd <- "colour"; gridAdd <- FALSE; landMarks <- FALSE; fileFormat <- "pdf"; thislocusName <- "ros_assembly_543443"
  # DataSet <- "Clines_Cadi"
  # minNum <- 10; markRoads <- T; plotPies <- T; poolsAdd <- T
  thisLocusName <- locusToPlot
  if (fileFormat=="eps") {
    filenamePiePlot <- paste(thisLocusName,"_TopoPiePlot.eps")
  }
  if (fileFormat=="pdf") {
    filenamePiePlot <- paste(thisLocusName,"_TopoPiePlot.pdf")
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
  plot(DEM_merged,main=chartLabel,ylim=yRange,xlim=xRange,alpha=alphaTerm,
       axes=FALSE,yaxt="n",xaxt="n",ann=F,cex=2,
       col=colPallete(colPalleteNum))
  thisLocusData <- DemeData[[thisLocusName]]
  # Get allele frequency data
  thisDemeMatrix2D_pCount <- thisLocusData[["DemeMatrix2D_pCount"]]
  #names(thisLocusData)
  thisDemeMatrix2D_pqCountNoZero <- thisLocusData[["DemeMatrix2D_combine_noZero"]]
  thisDemeMatrix2D_pqCountNoZero <- thisDemeMatrix2D_pqCountNoZero[thisDemeMatrix2D_pqCountNoZero[,"p"]!="NaN",]
  # sum(thisDemeMatrix2D_pqCountNoZero[,"pCount"])+ sum(thisDemeMatrix2D_pqCountNoZero[,"qCount"])
  # mark roads?
  if (markRoads==T) {
    plot(lowerRoad_UTM,add=T,col='black',lwd=2)
    plot(lowerRoad_new_UTM,add=T,col='black',lwd=2)
    plot(lowerRoad_1_UTM,add=T,col='black',lwd=2)
    plot(lowerRoad_2_UTM,add=T,col='black',lwd=2)
    plot(lowerRoad_3_UTM,add=T,col='black',lwd=2)
    plot(upperRoad_new_UTM,add=T,col='black',lwd=2)
    plot(upperRoad_UTM,add=T,col='black',lwd=2)
    # plot(upperRoad_1_UTM,add=T,col='black',lwd=2)
    # plot(upperRoad_UTM,add=T,col='black',lwd=2)
    
    plot(upperRoad_2_UTM,add=T,col='black',lwd=2)
    plot(upperRoad_3_UTM,add=T,col='black',lwd=2)
    plot(upperRoad_4_UTM,add=T,col='black',lwd=2)
    plot(upperRoad_5_UTM,add=T,col='black',lwd=2)
    plot(upperRoad_6_UTM,add=T,col='black',lwd=2)
    plot(upperRoad_7_UTM,add=T,col='black',lwd=2)
    plot(upperRoad_8_UTM,add=T,col='black',lwd=2)
    plot(upperRoad_9_UTM,add=T,col='black',lwd=2)
    plot(upperRoad_10_UTM,add=T,col='black',lwd=2)
    
    #plot(diagonal_UTM,add=T,col='black',lwd=2)
  }
  if (pointsInd=="gray")  {
    SNPsFullDataPheno <- SNPsFullData[SNPsFullData[,"phenoCat_final"]!=-9,]
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="Y","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="Y","Northing"])), 
           pch = 21, cex=ptSize,col=rgb(0,0,0,100,maxColorValue=255),bg="gray")
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="FR","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="FR","Northing"])), 
           pch = 21, cex=ptSize,col=rgb(0,0,0,100,maxColorValue=255),bg="gray")
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="FO","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="FO","Northing"])), 
           pch = 21, cex=ptSize,col=rgb(0,0,0,100,maxColorValue=255),bg="gray")
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="WR","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="WR","Northing"])), 
           pch = 21, cex=ptSize,col=rgb(0,0,0,100,maxColorValue=255),bg="gray")
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="W","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="W","DistNorthOfCentre"])), 
           pch = 21, cex=ptSize,col="black",bg="gray")
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="WO","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="WO","Northing"])), 
           pch = 21, cex=ptSize,col=rgb(0,0,0,100,maxColorValue=255),bg="gray")
  }
  if (pointsInd=="colour") {
    #write.csv(SNPsFullData,"SNPsFullData.csv")
    #SNPsFullData
    SNPsFullDataPheno <- SNPsFullData[SNPsFullData[,"phenoCat_final"]!=-9,]
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="Y","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="Y","Northing"])),
           pch = 21, cex=ptSize,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(255,255,0,255,maxColorValue=255))
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="FR","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="FR","Northing"])),
           pch = 21, cex=ptSize,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(255,0,0,255,maxColorValue=255))
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="FO","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="FO","Northing"])),
           pch = 21, cex=ptSize,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(255,165,0,255,maxColorValue=255))
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="WR","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="WR","Northing"])),
           pch = 21, cex=ptSize,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(255,195,203,255,maxColorValue=255))
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="W","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="W","Northing"])),
           pch = 21, cex=ptSize,col="black",bg=rgb(255,255,255,255,maxColorValue=255))
    points(cbind(as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="WO","Easting"]),
                 as.vector(SNPsFullDataPheno[SNPsFullDataPheno[,"phenoCat_final"]=="WO","Northing"])),
           pch = 21, cex=ptSize,col=rgb(0,0,0,100,maxColorValue=255),bg=rgb(238,216,174,255,maxColorValue=255))
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
  
  rownames(thisDemeMatrix2D_pqCountNoZero) <- 1:nrow(thisDemeMatrix2D_pqCountNoZero)
  
  # find axes in metres
  minX <- min(thisDemeMatrix2D_pqCountNoZero[,"X"])
  XinMetres <- thisDemeMatrix2D_pqCountNoZero[,"X"]-minX
  minY <- min(thisDemeMatrix2D_pqCountNoZero[,"Y"])
  YinMetres <- sort(thisDemeMatrix2D_pqCountNoZero[,"Y"]-minY)
  # max(thisDemeMatrix2D_pqCountNoZero[,"Y"])
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
  # easting axis
  axis(1,at=seq(410000,436000,2000),labels=NA,outer=F,col.axis="black",las=1,cex.axis=1.5,lwd=2.7,pos=max(Northing)+700, tck = 0.01)
  axis(1,at=seq(410000,436000,2000),labels=seq(410000,436000,2000),col.axis="black", las=1,cex.axis=1.8,lwd=0,pos=max(Northing)+1500)
  # northing axis
  axis(2,at=seq(4682000,4694000,2000),labels=NA,col.axis="black",las=1,cex.axis=1.5,lwd=2.7,pos=min(Easting)-200,tck = 0.01)
  axis(2,at=seq(4682000,4694000,2000),labels=NA,col.axis="black",las=1,cex.axis=1.5,lwd=2.7,pos=max(Easting)+200,tck = 0.01)
  axis(2,at=seq(4682000,4694000,2000),labels=seq(4682000,4694000,2000),col.axis="black", las=1,cex.axis=1.8,lwd=0,pos=max(Easting)+2200)
  #rect(xleft=-0.2, ybottom=39, xright=4.2, ytop=41, col='white',density = NULL, border = NA)
  #rect(xleft=-0.2, ybottom=44, xright=4.2, ytop=45, col='white',density = NULL, border = NA)
  #axis(2, at=c(41,41.5,42,42.5,43,43.5,44), labels=NA,col.axis="black", las=1,cex.axis=1.5,lwd=2.7,pos=0, tck = -0.01)
  #axis(2, at=c(41.5,42.5,43.5), labels=NA,col.axis="black", las=1,cex.axis=1.5,lwd=2.7,pos=0, tck = -0.005)
  #yLabels <-c(expression(paste(41^o,"N",sep=" ")),expression(paste(42^o,"N",sep=" ")),expression(paste(43^o,"N",sep=" ")),
  #            expression(paste(44^o,"N",sep=" ")))
  #xLabels <-c(expression(paste(0^o,"E",sep=" ")),expression(paste(1^o,"E",sep=" ")),expression(paste(2^o,"E",sep=" ")),
  #            expression(paste(3^o,"E",sep=" ")),expression(paste(4^o,"E",sep=" ")))
  #axis(2, at=c(41,42,43,44), labels=yLabels,col.axis="black", las=1,cex.axis=1.7,lwd=0,pos=0.09, tck = 0)
  axis(1, at=labelPos_Xeasting, labels=NA,col.axis="black", las=1,cex.axis=1.5,lwd=2.7,pos=min(Northing_keep)-210, tck = -0.01)
  axis(1, at=labelPos_Xeasting, labels=labelPos_Xmetres/1000,col.axis="black", las=1,cex.axis=1.8,lwd=0,pos=min(Northing_keep)-410, tck = -0.01)
  
  if (DataSetName=="Clines_Cadi") {
    axis(2, at=labelPos_Ynorthing, labels=NA,col.axis="black", las=1,cex.axis=1.5,lwd=2.7,pos=min(Easting_keep), tck = -0.01)
    axis(2, at=labelPos_Ynorthing, labels=labelPos_Ymetres/1000,col.axis="black", las=1,cex.axis=1.8,lwd=0,pos=min(Easting_keep), tck = -0.01)
  }
  # plotting pie charts
  totalN <- (thisDemeMatrix2D_pqCountNoZero[,"pCount"]+thisDemeMatrix2D_pqCountNoZero[,"qCount"])/2 
  keep <- which(totalN>=minNum)
  thisDemeMatrix2D_pqCountNoZero <- thisDemeMatrix2D_pqCountNoZero[keep,]
  
  if (plotPies==T) {
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
      #cat("\n .",thisDeme2D)
      #cat("..P/Q:",c(thisP,thisQ))
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
  }
  # mark pools?
  if (poolsAdd==TRUE) {
    points(x=411635.27,y=4690296.66,pch=23,col='black',bg=c(rgb(255,255,255,230,maxColorValue=255)),lwd=1.5,cex=2) # YP4 - pool 1
    points(x=422120.82,y=4686387.46,pch=23,col='black',bg=c(rgb(255,255,255,230,maxColorValue=255)),lwd=1.5,cex=2) # YP - pool 2
    points(x=421968.27,y=4686511.65,pch=23,col='black',bg=c(rgb(255,255,255,230,maxColorValue=255)),lwd=1.5,cex=2) # YP - pool 3
    points(x=424440.68,y=4686082.79,pch=23,col='black',bg=c(rgb(255,255,255,230,maxColorValue=255)),lwd=1.5,cex=2) # MP - pool 4
    points(x=425130.38,y=4685954.23,pch=23,col='black',bg=c(rgb(255,255,255,230,maxColorValue=255)),lwd=1.5,cex=2) # MP - pool 5
    points(x=431641.97,y=4686865.38,pch=23,col='black',bg=c(rgb(255,255,255,230,maxColorValue=255)),lwd=1.5,cex=2) # MP - pool 6
  }
  dev.off()  
}