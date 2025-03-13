# plotting different transects for cline fitting

plotHZ_clineTransect <- function(demeData,demeSize,lnLLsummary_res,locusPlot,outputFormat,plotX,plotY,plotAllTransects) {
  # demeData <- HZ_Planoles_200mDemes
  # locusPlot <- "ros_assembly_543443"
  # outputFormat <- "pdf"
  # plotX <- c(-3000,25000)
  # plotY <- c(-1000,7400)
  # plotAllTransects <- T
  
  par(mfrow=c(1,1))
  par(mar=c(3,3,3,3))
  par(oma=c(3,3,3,3))
  tempData <- demeData
  plot(tempData[,"X_adj"],tempData[,"Y_adj"],type='n',
       main="",cex.main=1,bty="n",yaxt="n",xaxt="n",
       xlab = "",ylab = "",col = "black", bg = "black",xaxs = "i", yaxs = "i",
       xlim=c(xMin-2000,24000),
       ylim=c(yMin-4000,10000))
  axis(1,at=seq(xMin-2000,xMax,1000),labels=NA,col.axis="black",las=2,cex.axis=1.1,lwd=1.6,pos=yMin-4000, tck = -.02)
  axis(1,at=seq(xMin-2000,xMax,2000),labels=seq(xMin-2000,xMax,2000)/1000,col.axis="black",tck=-0.3,las=1,cex.axis=1.1,lwd=0,pos=yMin-4000)
  axis(2,at=seq(yMin-4000,yMax+3800,1000),labels=NA,col.axis="black",las=1,cex.axis=1.1,lwd=1.6,pos=xMin-2000,tck = -.02)
  axis(2,at=seq(yMin-4000,yMax+3800,2000),labels=seq(yMin-4000,yMax+3800,2000)/1000,col.axis="black",tck=-.03,las=1,cex.axis=1.1,lwd=0,pos=xMin-2000.5)
  xMinLine <- min(seq(xMin-2000,xMax,2000))
  xMaxLine <- max(seq(xMin-2000,24000,2000))
  yMinLine <- min(seq(yMin-4000,yMax+3800,1000))
  yMaxLine <- max(seq(yMin-4000,yMax+3800,1000))
  lines(x=c(xMaxLine,xMaxLine),y=c(yMinLine,yMaxLine),lwd=1.6)
  lines(x=c(xMinLine-2000,xMaxLine),y=c(yMaxLine,yMaxLine),lwd=1.6)
  lines(x=c(xMaxLine-2000,xMaxLine),y=c(yMinLine,yMinLine),lwd=1.8)
  
  title(paste(locusPlot,paste(demeSize,"m demes",sep=""),sep=", "),line=2.3,cex.main=0.9)
  mtext("East (km) / Easting",side=1,cex=1,line=2.5)
  mtext("North (km) / Northing",side=2,cex=1,line=2.5)
  axis(1,at=seq(-435,24000,2000),labels=NA,outer=F,col.axis="black",las=2,cex.axis=1.1,lwd=1.6,pos=yMaxLine, tck = 0.02)
  axis(1,at=seq(-435,24000,2000),labels=seq(410000,435435,2000),col.axis="black",tck=-0.3,las=1,cex.axis=1.1,lwd=0,pos=yMaxLine+1050)
  axis(2,at=seq(-2246,10000,2000),labels=NA,col.axis="black",las=1,cex.axis=1.1,lwd=1.6,pos=xMaxLine,tck = 0.02)
  axis(2,at=seq(-2246,10000,2000),labels=seq(4682000,4694000,2000),col.axis="black",tck=-.03,las=1,cex.axis=1.1,lwd=0,pos=xMaxLine+2500)
  
  if (plotAllTransects==T) {
    for (thisCurve in 1:nrow(lnLLsummary_res)) {
      # thisCurve <- 1
      transectGradient <- lnLLsummary_res[thisCurve,7]
      transectIntercept <- lnLLsummary_res[thisCurve,8]
      curve(((transectGradient*x)+(transectIntercept)), from=-3000,
            to=25000, col="lightgray",add = TRUE,type = "l",lwd=2,lty=1)
    }
  }
  
  # Identify best fitting and highlight
  rowMaxLL <- which(lnLLsummary_res[,"maxLL"]==max(lnLLsummary_res[,"maxLL"]))
  
  transectGradient <- lnLLsummary_res[rowMaxLL,7]
  transectIntercept <- lnLLsummary_res[rowMaxLL,8]
  curve(((transectGradient*x)+(transectIntercept)), from=-3000,
        to=25000, col="black",add = TRUE,type = "l",lwd=2.5,lty=1)
  
  # plotting pie charts
  totalN <- (tempData[,"pCount"]+tempData[,"qCount"])/2
  keep <- which(totalN>=10)
  tempData <- tempData[keep,]
  allthisP <- as.numeric(as.vector(tempData[,"pCount"]))
  allthisQ <- as.numeric(as.vector(tempData[,"qCount"]))
  allSum <- allthisP+allthisQ
  allFreq <- allthisP/allSum
  
  if (mean(allFreq[1:4])>0.6) {
    allthisP_temp <- allthisQ
    allthisQ_temp <- allthisP
    allthisQ <- allthisP_temp
    allthisP <- allthisQ_temp
  }
  for (thisDeme2D in 1:nrow(tempData)) {
    # thisDeme2D <- 1
    thisXraw <- as.numeric(as.vector(tempData[thisDeme2D,"X_adj"]))
    thisYraw <- as.numeric(as.vector(tempData[thisDeme2D,"Y_adj"]))
    thisP <- as.numeric(allthisP[thisDeme2D])
    thisQ <- as.numeric(allthisQ[thisDeme2D])
    # pie plot size has to be tweaked depending on deme size
    # small pie plots, original was radius 69 for deme size 200
    #     if (sum(thisP,thisQ)>minNum) {
    #       add.pie(z=c(thisP,thisQ), x=thisXmid, y=thisYmid, radius=69,col=c((rgb(255,0,0,255,maxColorValue=255)),(rgb(255,255,0,255,maxColorValue=255))),labels="")
    #     }
    # larger pie plots - 120 works okay for deme size 400
    totalN <- thisP + thisQ
    #if ((totalN < minNum)==F) {
    if (demeSize==200) {
      add.pie(z=c(thisP,thisQ), x=thisXraw, y=thisYraw, radius=80, col=c((rgb(255,0,0,255,maxColorValue=255)), (rgb(255,225,0,255,maxColorValue=255))),labels="")
    }
    if (demeSize==100 | demeSize==50) {
      add.pie(z=c(thisP,thisQ), x=thisXraw, y=thisYraw, radius=60, col=c((rgb(255,0,0,255,maxColorValue=255)), (rgb(255,225,0,255,maxColorValue=255))),labels="")
    }
  }
  
    }
