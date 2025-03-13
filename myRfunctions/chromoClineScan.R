# Chromosome Scan figures for KASP ----
genomeScanPlotWidth_KASP <- function(thisData,deltaP,plotPos1,plotPos2,startPlot,minPlot,maxPlot,widthMax,
                                      stepPlot,posData,LG,textPt,thisData_ref,plotColGenes) {
  # thisData <- clineFits_LG6_ros_assem; thisData_ref <- clineFits_ref
  # plotPos1=c(0.03,0.9,0.18,0.35);plotPos2=c(0.86,0.98,0.17,0.35); startPlot=T; LG=6
  # minPlot=0;maxPlot=0.9;widthMax=14;stepPlot=0.1; posData="position_thisScaff"
  # deltaP <- 0.6
  thisData <- thisData[thisData[,"LG"]==LG,]
  cat("\n---------------------------")
  
  # Sulf 34987000-34802000
  if (LG==4) {
    thisData_sulf <- thisData[thisData[,posData]>34802000 & thisData[,posData]<34987000,]
    #cat("\n sulf range positions: ",c(min(thisData_sulf[,"position"]),max(thisData_sulf[,"position"])))
    cat("\nClines sulf region total, n =",nrow(thisData_sulf))
    cat("\nClines thisData total, n =",nrow(thisData))
  }
  if (LG==6) {
    thisData_RosEl <- thisData[thisData[,"scaffold"]=="ros_assembly",]
    #cat("\n sulf range positions: ",c(min(thisData_sulf[,"position"]),max(thisData_sulf[,"position"])))
    cat("\nClines RosEl region total, n =",nrow(thisData_RosEl))
    cat("\nClines thisData total, n =",nrow(thisData))
  }
  Xvals <- seq(minPlot,maxPlot,stepPlot)
  Xvals_smallTick <- seq(minPlot,maxPlot,stepPlot/10)
  #yOrigin <- floor((max(thisData[,"width"])+3000)/1000)
  # plot
  par(fig = plotPos1, mar=c(1,1,1,1),new=((startPlot==TRUE)==FALSE))
  plot(NA, xlab="", xaxt="n",
       ylab="", main="", las=2, bty="n",cex=0.5,yaxt="n",xaxt="n",type='n',
       xlim=c(minPlot,maxPlot), ylim=rev(c(-1.5,widthMax+3)),
       cex.axis=1.4,cex.lab=1.4,col=rgb(100,100,100,120,maxColorValue=255), pch=16)
  # axis(2, at=rev(c(0,2,4,6,8,10,12,14,16,18,20)), labels=NA,col.axis="black", las=1,cex.axis=1,lwd=2,pos=min(Xvals), tck = -.05)
  # axis(2, at=rev(c(0,4,8,12,16)), labels=rev(c(0,4,8,12,16)),col.axis="black", las=1,cex.axis=1,lwd=0,pos=min(Xvals))
  # axis(2, at=18, labels="n.f.",col.axis="black", las=1,cex.axis=1,lwd=0,pos=min(Xvals))
  
  axis(2, at=rev(seq(0,widthMax+2,2)), labels=NA,col.axis="black", las=1,cex.axis=1,lwd=2,pos=min(Xvals), tck = -.05)
  axis(2, at=rev(seq(0,widthMax-2,2)), labels=rev(seq(0,widthMax-2,2)),col.axis="black", las=1,cex.axis=1,lwd=0,pos=min(Xvals))
  axis(2, at=widthMax, labels="n.f.",col.axis="black", las=1,cex.axis=1,lwd=0,pos=min(Xvals))
  axis(1, at=Xvals,labels=NA,col.axis="black", las=1,cex.axis=1.2,lwd=2,pos=widthMax+2, tck = -.05)
  axis(1, at=Xvals_smallTick,labels=NA,col.axis="black", las=1,cex.axis=1.2,lwd=1.4,pos=widthMax+2, tck = -.02)
  axis(1, at=Xvals,labels=Xvals,col.axis="black", las=1,cex.axis=1,lwd=0,pos=widthMax+2+0.05)
  mtext("width",side=2,cex=1.4,line=0.8)
  
  # WidthQuantiles <- quantile(thisData[,"width_est2"], na.rm=T, c(0.5, 0.05))
  # lines(x=c(0,max(range(thisData[,"position"]))),y=c(as.numeric(WidthQuantiles[1]),as.numeric(WidthQuantiles[1]))/1000,col = "black", lwd=1.5,lty = 3)
  # lines(x=c(0,max(range(thisData[,"position"]))),y=c(as.numeric(WidthQuantiles[2]),as.numeric(WidthQuantiles[2]))/1000,col = "black", lwd=1.5,lty = 3)
  # 
  if (plotColGenes==T) {
    # ROS1
    rect(xleft=(541738/1000000),xright=(544102/1000000),ybottom=0,ytop=widthMax+1.5, col='pink',lty = "solid", lwd=2, border='pink')
    # ROS2
    rect(xleft=(549952/1000000),xright=(567141/1000000),ybottom=0,ytop=widthMax+1.5, col='pink',lty = "solid", lwd=2, border='pink')
    # ROS3
    rect(xleft=(573791/1000000),xright=(576279/1000000),ybottom=0,ytop=widthMax+1.5, col='pink',lty = "solid", lwd=2, border='pink')
    # Eluta
    rect(xleft=(705682/1000000),xright=(716427/1000000),ybottom=0,ytop=widthMax+1.5, col='pink',lty = "solid", lwd=2, border='pink')
  }
    
  if (nrow(thisData_ref)>0) {
    lines(x=c(min(Xvals),max(Xvals)),y=c(as.numeric(thisData_ref[1,"width"])/1000,as.numeric(thisData_ref[1,"width"])/1000)
          ,col = "darkgray", lwd=3,lty = 2)
  }
  # add points

  for (thisLoc in 1:nrow(thisData)) {
    # thisLoc <- 7
    if (thisData[thisLoc,"width"]!=0) {
      lines(x=c(thisData[thisLoc,posData]/1000000,thisData[thisLoc,posData]/1000000),
            y=c(as.numeric(thisData[thisLoc,"width_U95"])/1000,as.numeric(thisData[thisLoc,"width_L95"])/1000),lwd=1,lty=1)
      points(thisData[thisLoc,posData]/1000000,
             as.numeric(thisData[thisLoc,"width"])/1000,
             col='black',bg="black", lwd=1, pch=21,cex=1)
      if ((as.numeric(thisData[thisLoc,"width"])/1000)>widthMax) {
        points(thisData[thisLoc,posData]/1000000,15.2,
               col='black',bg="black", lwd=1, pch=21,cex=1)
      }
    }
    if (thisData[thisLoc,"width"]==0) {
      points(thisData[thisLoc,posData]/1000000,widthMax,
             col='black',bg=rgb(20,20,20,20,maxColorValue=255), lwd=1, pch=21,cex=1)
    }
  }
  # reference
  if (nrow(thisData_ref)>0) {
    if (LG!=6) { # to ensure reference is not placed down on ros_assembly plots
      lines(x=c(max(maxPlot)-0.05,max(maxPlot)-0.05),
            y=c(as.numeric(thisData_ref[1,"width_U95"])/1000,as.numeric(thisData_ref[1,"width_L95"])/1000),lwd=1,lty=1)
      points((max(maxPlot)-0.05),as.numeric(thisData_ref[1,"width"])/1000,col='black',bg="red", lwd=1, pch=21,cex=1)
    }
   
  }
  if (textPt==T) {
    xPos <- thisData[,posData]/1000000
    yPos <-as.numeric(thisData[,"width"])/1000
    yPos[yPos==0] <- widthMax-0.8
    textxy(X=xPos,Y=yPos,labs=as.vector(thisData[,"numberPlot"]),cex=1,offset=1)
  }
}
genomeScanPlotCentre_KASP <- function(thisData,deltaP,plotPos1,plotPos2,startPlot,minPlot,maxPlot,centreMax,
                                       stepPlot,posData,LG,textPt,thisData_ref,plotColGenes) {
  # plotPos1=c(0.03,0.9,0.02,0.20);plotPos2=c(0.86,0.98,0.02,0.20); startPlot=T
  # thisData <- LG4_all; minPlot=0;maxPlot=0.9;centreMax=20;stepPlot<- 0.1
  # deltaP <- 0.6
  # thisData <- clineFits
  # plotPos1 <- c(0.03,0.9,0.84,1.000); plotPos2 <- c(0.86,0.98,0.84, 1); startPlot=T
  # minPlot=0;maxPlot=50;stepPlot=2
  # thisData <- LG4_all
  # deltaP <- 0.6
  # if (ref!=F) {
  #   thisData_ref <- thisData[thisData[,"LG"]==6,]
  #   thisData_ref <- thisData_ref[thisData_ref[,"locusName"]==ref,]
  # }
  thisData <- thisData[thisData[,"LG"]==LG,]
  cat("\n---------------------------")
  # Sulf 34987000-34802000
  if (LG==4) {
    thisData_sulf <- thisData[thisData[,posData]>34802000 & thisData[,posData]<34987000,]
    #cat("\n sulf range positions: ",c(min(thisData_sulf[,"position"]),max(thisData_sulf[,"position"])))
    cat("\nClines sulf region total, n =",nrow(thisData_sulf))
    cat("\nClines thisData total, n =",nrow(thisData))
  }
  if (LG==6) {
    thisData_RosEl <- thisData[thisData[,"scaffold"]=="ros_assembly",]
    #cat("\n sulf range positions: ",c(min(thisData_sulf[,"position"]),max(thisData_sulf[,"position"])))
    cat("\nClines RosEl region total, n =",nrow(thisData_RosEl))
    cat("\nClines thisData total, n =",nrow(thisData))
  }
  Xvals <- seq(minPlot,maxPlot,stepPlot)
  Xvals_smallTick <- seq(minPlot,maxPlot,stepPlot/10)
  # plot
  par(fig = plotPos1, mar=c(1,1,1,1),new=((startPlot==TRUE)==FALSE))
  plot(NA, xlab="", xaxt="n",
       ylab="", main="", las=2, bty="n",cex=0.5,yaxt="n",xaxt="n",type='n',
       xlim=c(minPlot,maxPlot), ylim=c(-2,centreMax+1),
       cex.axis=1.4,cex.lab=1.4,col=rgb(100,100,100,120,maxColorValue=255), pch=16)
  # axis(2, at=c(-2,0,2,4,6,8,10,12,14,16,18,20), labels=NA,col.axis="black", las=1,cex.axis=1,lwd=2,pos=min(Xvals), tck = -.05)
  # axis(2, at=c(4,8,12,16,20), labels=c(4,8,12,16,20),col.axis="black", las=1,cex.axis=1,lwd=0,pos=min(Xvals))
  # axis(2, at=0, labels="n.f.",col.axis="black", las=1,cex.axis=1,lwd=0,pos=min(Xvals))
  # axis(1, at=Xvals,labels=NA,col.axis="black", las=1,cex.axis=1.2,lwd=2,pos=-2, tck = -.05)
  # axis(1, at=Xvals_smallTick,labels=NA,col.axis="black", las=1,cex.axis=1.2,lwd=1.4,pos=-2, tck = -.02)
  # axis(1, at=Xvals,labels=Xvals,col.axis="black", las=1,cex.axis=1,lwd=0,pos=-2.5)
  axis(2, at=seq(-2,centreMax,2), labels=NA,col.axis="black", las=1,cex.axis=1,lwd=2,pos=min(Xvals), tck = -.05)
  axis(2, at=seq(2,centreMax,2), labels=seq(2,centreMax,2),col.axis="black", las=1,cex.axis=1,lwd=0,pos=min(Xvals))
  axis(2, at=0, labels="n.f.",col.axis="black", las=1,cex.axis=1,lwd=0,pos=min(Xvals))
  axis(1, at=Xvals,labels=NA,col.axis="black", las=1,cex.axis=1.2,lwd=2,pos=-2, tck = -.05)
  axis(1, at=Xvals_smallTick,labels=NA,col.axis="black", las=1,cex.axis=1.2,lwd=1.4,pos=-2, tck = -.02)
  axis(1, at=Xvals,labels=Xvals,col.axis="black", las=1,cex.axis=1,lwd=0,pos=-2.3)
  mtext("centre",side=2,cex=1.4,line=0.8)
  if (plotColGenes==T) {
    # ROS1
    rect(xleft=(541738/1000000),xright=(544102/1000000),ybottom=-1.5,ytop=centreMax, col='pink',lty = "solid", lwd=2, border='pink')
    # ROS2
    rect(xleft=(549952/1000000),xright=(567141/1000000),ybottom=-1.5,ytop=centreMax, col='pink',lty = "solid", lwd=2, border='pink')
    # ROS3
    rect(xleft=(573791/1000000),xright=(576279/1000000),ybottom=-1.5,ytop=centreMax, col='pink',lty = "solid", lwd=2, border='pink')
    # Eluta
    rect(xleft=(705682/1000000),xright=(716427/1000000),ybottom=-1.5,ytop=centreMax, col='pink',lty = "solid", lwd=2, border='pink')
  }
  
  # WidthQuantiles <- quantile(thisData[,"width_est2"], na.rm=T, c(0.5, 0.05))
  # lines(x=c(0,max(range(thisData[,"position"]))),y=c(as.numeric(WidthQuantiles[1]),as.numeric(WidthQuantiles[1]))/1000,col = "black", lwd=1.5,lty = 3)
  # lines(x=c(0,max(range(thisData[,"position"]))),y=c(as.numeric(WidthQuantiles[2]),as.numeric(WidthQuantiles[2]))/1000,col = "black", lwd=1.5,lty = 3)
  # 
  # add points
  # ref line background
  if (nrow(thisData_ref)>0) {
    lines(x=c(min(Xvals),max(Xvals)),y=c(as.numeric(thisData_ref[1,"centre"])/1000,as.numeric(thisData_ref[1,"centre"])/1000)
          ,col = "darkgray", lwd=3,lty = 2)
  }
  centreRef <- as.numeric(thisData_ref[1,"centre"])/1000
  for (thisLoc in 1:nrow(thisData)) {
    # thisLoc <- 7
    if (thisData[thisLoc,"centre"]!=0) {
      lines(x=c(thisData[thisLoc,posData]/1000000,thisData[thisLoc,posData]/1000000),
            y=c(as.numeric(thisData[thisLoc,"centre_U95"])/1000,as.numeric(thisData[thisLoc,"centre_L95"])/1000),lwd=1,lty=1)
      points(thisData[thisLoc,posData]/1000000,
             as.numeric(thisData[thisLoc,"centre"])/1000,
             col='black',bg="black", lwd=1, pch=21,cex=1)
    }
    if (thisData[thisLoc,"centre"]==0) {
      points(thisData[thisLoc,posData]/1000000,0,
             col='black',bg=rgb(20,20,20,20,maxColorValue=255), lwd=1, pch=21,cex=1)
    }
  }
  if (nrow(thisData_ref)>0) {
    lines(x=c(max(maxPlot)-0.05,max(maxPlot)-0.05),
          y=c(as.numeric(thisData_ref[1,"centre_U95"])/1000,as.numeric(thisData_ref[1,"centre_L95"])/1000),lwd=1,lty=1)
    points(max(maxPlot)-0.05,
           as.numeric(thisData_ref[1,"centre"])/1000,
           col='black',bg="red", lwd=1, pch=21,cex=1)
  }
  if (textPt==T) {
    xPos <- thisData[,posData]/1000000
    yPos <-as.numeric(thisData[,"centre"])/1000
    yPos[yPos==0] <- 0.5
    textxy(X=xPos,Y=yPos,labs=as.vector(thisData[,"numberPlot"]),cex=1,offset=1)
  }
  # #Position second plot
  # #par(fig=c(0.85,1,0.84, 1), new=T, bty="n", mar=c(5,5,4,1))
  # par(fig=plotPos2, new=T, bty="n", mar=c(1,1,1,1))
  # yhist = hist(as.numeric(thisData[,"centre_est2"])/1000,plot=F,breaks=seq(0,22,0.5))
  # # line frequency distribution for background whole genome first
  # minClineBin <- min(yhist$mids[yhist$counts!=0])
  # yhist_counts <- yhist$counts[yhist$counts!=0]
  # yhist_mids <- yhist$mids[yhist$counts!=0]
  # plot(yhist_counts, yhist_mids, col="black",lwd=3,xlim=c(-1,max(yhist$counts)),ylim=c(0,22),
  #      type="l", xlab="",  ylab="",  bty="n",  xaxs = "i", yaxs = "i", xaxt="n", yaxt="n")
  # secondLab <- paste("delta P>",deltaP,sep="")
  # mtext(paste(paste("n=",nrow(thisData),sep=""),secondLab,sep="\n"),cex=0.8,side=1,line=-4)
}
