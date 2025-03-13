# Likelihood surfaces
# Only using the accepted values 
# Problems can arise when the max value is much lower than any values kept after thinning.
#   i.e. by chance, all the values near the max have been removed, so that there is only vals > -10 from max lnL

LL_credibleRegions <- function(SANN_output,locusName,SANN_clinefit_retvals_accept,modelSummary) {
  # SANN_output <- SANNfit_thisLocus
  par(mfrow=c(2,3))  
  par(mar=c(5, 5, 5, 2.1))  
  par(oma=c(1,1,1,1)) 
  # Accept and reject values
  SANN_clinefit_retvals <- SANN_output$retvals
  SANN_clinefit_retvals_reject <- SANN_clinefit_retvals[SANN_clinefit_retvals[,"accept"]==0,]
  SANN_clinefit_retvals_accept <- SANN_clinefit_retvals[SANN_clinefit_retvals[,"accept"]==1,]
  # plot(SANN_clinefit_retvals_accept[,"p1"],SANN_clinefit_retvals_accept[,"val"])
  # plot(SANN_clinefit_retvals_accept[,"p2"],SANN_clinefit_retvals_accept[,"val"])
  
  # SANN_clinefit_retvals[which(SANN_clinefit_retvals[,"val"]==max(SANN_clinefit_retvals[,"val"])),]
  
  # centre
  # within10LogL <- SANN_clinefit_retvals[SANN_clinefit_retvals[,"val"]>(SANN_output$minimum[1]-10),]
  # within2lnL <- SANN_clinefit_retvals[SANN_clinefit_retvals_accept[,"val"]>(SANN_output$minimum[1]-2),]
  within10LogL <- SANN_clinefit_retvals_accept[SANN_clinefit_retvals_accept[,"val"]>(max(SANN_clinefit_retvals[,"val"])-10),]
  within2lnL <- SANN_clinefit_retvals_accept[SANN_clinefit_retvals_accept[,"val"]>(max(SANN_clinefit_retvals[,"val"])-2),"p1"]
  lowerVals <- min(within2lnL)
  upperVals <- max(within2lnL)
  modelSummary[modelSummary[,"locus"]==locusName,"centre_L95"] <- lowerVals
  modelSummary[modelSummary[,"locus"]==locusName,"centre_U95"] <- upperVals
  centreVal <- as.numeric(as.vector(modelSummary[modelSummary[,"locus"]==locusName,"centre_accept"]))
  plot(SANN_clinefit_retvals[SANN_clinefit_retvals[,"val"]>(SANN_output$minimum[1]-10),"p1"],
       SANN_clinefit_retvals[SANN_clinefit_retvals[,"val"]>(SANN_output$minimum[1]-10),"val"],
       main=paste(locusName, " - centre",sep=" "),xlab="centre (m)",ylab="lnL",cex=0.9,cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
  valsSubHead2 <- c(paste("lnLmax=",
    round(as.numeric(as.vector(modelSummary[modelSummary[,"locus"]==locusName,"centre_accept"])),0),sep=""),
                    paste("(",round(lowerVals,0),sep=""),"-",
                    paste(round(upperVals,0),")",sep=""))
  mtext(paste(valsSubHead2,collapse=" "),cex=0.8)
  # lines(x=c(lowerVals,lowerVals), y=c(min(within10LogL[,"val"]),max(within10LogL[,"val"])), 
  #       type="l", pch=22, col="blue",lty=3,lwd=2.5)  
  # lines(x=c(upperVals,upperVals), y=c(min(within10LogL[,"val"]),max(within10LogL[,"val"])), 
  #       type="l", pch=22, col="blue",lty=3,lwd=2.5)  
  # lines(x=c(centreVal,centreVal), y=c(min(within10LogL[,"val"]),max(within10LogL[,"val"])), 
  #       type="l", pch=22, col="blue",lty=1,lwd=3)  
  
  lines(x=c(lowerVals,lowerVals), y=c(min(within10LogL[,"val"]),min(within10LogL[,"val"])-2), 
        type="l", pch=22, col="blue",lty=2,lwd=2.5)  
  lines(x=c(upperVals,upperVals), y=c(min(within10LogL[,"val"]),min(within10LogL[,"val"])-2), 
        type="l", pch=22, col="blue",lty=2,lwd=2.5)  
  lines(x=c(centreVal,centreVal), y=c(min(within10LogL[,"val"]),min(within10LogL[,"val"])-2), 
        type="l", pch=22, col="red",lty=1,lwd=5)  
  
  # width Likelihood surface
  within10LogL <- SANN_clinefit_retvals_accept[SANN_clinefit_retvals_accept[,"val"]>(SANN_output$minimum[1]-10),]
  within2lnL <- SANN_clinefit_retvals_accept[SANN_clinefit_retvals_accept[,"val"]>(SANN_output$minimum[1]-2),"p2"]
  lowerVals <- min(within2lnL)
  upperVals <- max(within2lnL)
  modelSummary[modelSummary[,"locus"]==locusName,"width_L95"] <- lowerVals
  modelSummary[modelSummary[,"locus"]==locusName,"width_U95"] <- upperVals

  widthVal <- as.numeric(as.vector(modelSummary[modelSummary[,"locus"]==locusName,"width_accept"]))
  plot(SANN_clinefit_retvals[SANN_clinefit_retvals[,"val"]>(SANN_output$minimum[1]-10),"p2"],
       SANN_clinefit_retvals[SANN_clinefit_retvals[,"val"]>(SANN_output$minimum[1]-10),"val"],
       main=paste(locusName, " - width",sep=" "),xlab="width (m)",ylab="lnL",cex=0.9,cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
  valsSubHead2 <- c(paste("lnLmax=",
                          round(as.numeric(as.vector(modelSummary[modelSummary[,"locus"]==locusName,"width_accept"])),0),sep=""),
                    paste("(",round(lowerVals,0),sep=""),"-",
                    paste(round(upperVals,0),")",sep=""))
  mtext(paste(valsSubHead2,collapse=" "),cex=0.8)
  lines(x=c(lowerVals,lowerVals), y=c(min(within10LogL[,"val"]),min(within10LogL[,"val"])-2), 
        type="l", pch=22, col="blue",lty=2,lwd=2.5)  
  lines(x=c(upperVals,upperVals), y=c(min(within10LogL[,"val"]),min(within10LogL[,"val"])-2), 
        type="l", pch=22, col="blue",lty=2,lwd=2.5)  
  lines(x=c(widthVal,widthVal), y=c(min(within10LogL[,"val"]),min(within10LogL[,"val"])-2), 
        type="l", pch=22, col="red",lty=1,lwd=5)  

  # pL Likelihood surface
  within10LogL <- SANN_clinefit_retvals_accept[SANN_clinefit_retvals_accept[,"val"]>(SANN_output$minimum[1]-10),]
  within2lnL <- SANN_clinefit_retvals_accept[SANN_clinefit_retvals_accept[,"val"]>(SANN_output$minimum[1]-2),"p3"]
  lowerVals <- min(within2lnL)
  upperVals <- max(within2lnL)
  modelSummary[modelSummary[,"locus"]==locusName,"pL_L95"] <- lowerVals
  modelSummary[modelSummary[,"locus"]==locusName,"pL_U95"] <- upperVals
  pLVal <- as.numeric(as.vector(modelSummary[modelSummary[,"locus"]==locusName,"pL_accept"]))
  plot(SANN_clinefit_retvals[SANN_clinefit_retvals[,"val"]>(SANN_output$minimum[1]-10),"p3"],
       SANN_clinefit_retvals[SANN_clinefit_retvals[,"val"]>(SANN_output$minimum[1]-10),"val"],
       main=paste(locusName, " - p0",sep=" "),xlab="p0",ylab="lnL",cex=0.9,cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
  valsSubHead2 <- c(paste("lnLmax=",
                          round(as.numeric(as.vector(modelSummary[modelSummary[,"locus"]==locusName,"pL_accept"])),3),sep=""),
                    paste("(",round(lowerVals,3),sep=""),"-",
                    paste(round(upperVals,3),")",sep=""))
  mtext(paste(valsSubHead2,collapse=" "),cex=0.8)
  lines(x=c(lowerVals,lowerVals), y=c(min(within10LogL[,"val"]),min(within10LogL[,"val"])-2), 
        type="l", pch=22, col="blue",lty=2,lwd=2.5)  
  lines(x=c(upperVals,upperVals), y=c(min(within10LogL[,"val"]),min(within10LogL[,"val"])-2), 
        type="l", pch=22, col="blue",lty=2,lwd=2.5)  
  lines(x=c(pLVal,pLVal), y=c(min(within10LogL[,"val"]),min(within10LogL[,"val"])-2), 
        type="l", pch=22, col="red",lty=1,lwd=5)  
  
  # pR Likelihood surface
  within10LogL <- SANN_clinefit_retvals_accept[SANN_clinefit_retvals_accept[,"val"]>(SANN_output$minimum[1]-10),]
  within2lnL <- SANN_clinefit_retvals_accept[SANN_clinefit_retvals_accept[,"val"]>(SANN_output$minimum[1]-2),"p4"]
  lowerVals <- min(within2lnL)
  upperVals <- max(within2lnL)
  modelSummary[modelSummary[,"locus"]==locusName,"pR_L95"] <- lowerVals
  modelSummary[modelSummary[,"locus"]==locusName,"pR_U95"] <- upperVals
  pRVal <- as.numeric(as.vector(modelSummary[modelSummary[,"locus"]==locusName,"pR_accept"]))
  plot(SANN_clinefit_retvals[SANN_clinefit_retvals[,"val"]>(SANN_output$minimum[1]-10),"p4"],
       SANN_clinefit_retvals[SANN_clinefit_retvals[,"val"]>(SANN_output$minimum[1]-10),"val"],
       main=paste(locusName, " - p1",sep=" "),xlab="p1",ylab="lnL",cex=0.9,cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
  valsSubHead2 <- c(paste("lnLmax=",
                          round(as.numeric(as.vector(modelSummary[modelSummary[,"locus"]==locusName,"pR_accept"])),3),sep=""),
                    paste("(",round(lowerVals,3),sep=""),"-",
                    paste(round(upperVals,3),")",sep=""))
  mtext(paste(valsSubHead2,collapse=" "),cex=0.8)
  lines(x=c(lowerVals,lowerVals), y=c(min(within10LogL[,"val"]),min(within10LogL[,"val"])-2), 
        type="l", pch=22, col="blue",lty=2,lwd=2.5)  
  lines(x=c(upperVals,upperVals), y=c(min(within10LogL[,"val"]),min(within10LogL[,"val"])-2), 
        type="l", pch=22, col="blue",lty=2,lwd=2.5)  
  lines(x=c(pRVal,pRVal), y=c(min(within10LogL[,"val"]),min(within10LogL[,"val"])-2), 
        type="l", pch=22, col="red",lty=1,lwd=5)  
  
  # Fst Likelihood surface
  within10LogL <- SANN_clinefit_retvals_accept[SANN_clinefit_retvals_accept[,"val"]>(SANN_output$minimum[1]-10),]
  within2lnL <- 1/(SANN_clinefit_retvals_accept[SANN_clinefit_retvals_accept[,"val"]>(SANN_output$minimum[1]-2),"p5"]+1)
  lowerVals <- min(within2lnL)
  upperVals <- max(within2lnL)
  modelSummary[modelSummary[,"locus"]==locusName,"Fst_L95"] <- lowerVals
  modelSummary[modelSummary[,"locus"]==locusName,"Fst_U95"] <- upperVals
  fstVal <- as.numeric(as.vector(modelSummary[modelSummary[,"locus"]==locusName,"Fst_accept"]))
  # Note, unlike parameters above - here a conversion from alpha to Fst required: Fst=1/(alpha+1)
  fst_p5 <- 1/(SANN_clinefit_retvals[SANN_clinefit_retvals[,"val"]>(SANN_output$minimum[1]-10),"p5"]+1)
  plot(fst_p5,SANN_clinefit_retvals[SANN_clinefit_retvals[,"val"]>(SANN_output$minimum[1]-10),"val"],
       main=paste(locusName, " - Fst",sep=" "),xlab="Fst",ylab="lnL",cex=0.9,cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
  valsSubHead2 <- c(paste("lnLmax=",
                          round(as.numeric(as.vector(modelSummary[modelSummary[,"locus"]==locusName,"Fst_accept"])),3),sep=""),
                    paste("(",round(lowerVals,3),sep=""),"-",
                    paste(round(upperVals,3),")",sep=""))
  mtext(paste(valsSubHead2,collapse=" "),cex=0.8)
  lines(x=c(lowerVals,lowerVals), y=c(min(within10LogL[,"val"]),min(within10LogL[,"val"])-2), 
        type="l", pch=22, col="blue",lty=2,lwd=2.5)  
  lines(x=c(upperVals,upperVals), y=c(min(within10LogL[,"val"]),min(within10LogL[,"val"])-2), 
        type="l", pch=22, col="blue",lty=2,lwd=2.5)  
  lines(x=c(fstVal,fstVal), y=c(min(within10LogL[,"val"]),min(within10LogL[,"val"])-2), 
        type="l", pch=22, col="red",lty=1,lwd=5)
  
  return(modelSummary)
}
