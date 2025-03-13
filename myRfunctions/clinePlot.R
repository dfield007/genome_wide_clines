# Individual Cline plot
# David Field
clinePlot <- function(clineModel,locusName,clineFit,curveFit,cexMain=1.1,cexAxisPlot= 1.85,
            axisLineWdth=2.5,cexPoints=1.8,lwdCurve=4,refLocus=NA,refPoints=F,HeaderDetails=F,FigLabel=NA,HaplotypeData) {
  # clineFit <- clineFits_test; curveFit=T; locusName="s1187_290152"; clineFit = modelSummary[thisRow,]
  # HaplotypeData <- HZ_Planoles_200mDemes_FLAr; # HeaderDetails <- T; FigLabel<- "(a)"
  # cexMain = 1.1; cexAxisPlot=1.85; axisLineWdth=2.5; cexPoints=1.8;lwdCurve=4;refLocus=NA;refPoints=F; 
  cline_model_sig_symm_plot <- function(XCordinate) {
    ((pL + ((pR-pL)/(1 + exp((-4*(XCordinate-centre))/width)))))
  }
  cline_model_sig_symm_plot_ref <- function(XCordinate_ref) {
    ((pL_ref + ((pR_ref-pL_ref)/(1 + exp((-4*(XCordinate_ref-centre_ref))/width_ref)))))
  }
  
  if (as.numeric(as.vector(clineFit[1,"N_total"]))==0) {
    next
  }
  if (as.numeric(as.vector(clineFit[1,"N_total"]))!=0) {
    
    if (HeaderDetails==T) {
      par(mfrow=c(1,1),mar=c(3.7,4.5,4.8,2),oma=c(1,1,1,1)) # (bottom, left, top, right)
    }
    if (HeaderDetails==F) {
      par(mfrow=c(1,1),mar=c(4.5,2,1,1),oma=c(2,2.5,2,1)) # (bottom, left, top, right)
    }
    deltaP <- round(as.numeric(as.vector(clineFit[1,"delta_p"])),2)
    TotalSamples <- as.numeric(as.vector(clineFit[1,"N_total"]))
    XCordinates <- as.numeric(as.vector(unlist(strsplit(as.vector(clineFit[1,"Xvals"]),";"))))
    p_adj <- as.numeric(as.vector(unlist(strsplit(as.vector(clineFit[1,"p"]),";"))))
    Ngenes <- as.numeric(as.vector(unlist(strsplit(as.vector(clineFit[1,"Ngenes"]),";"))))
    Nind <- Ngenes/2
    N_demes_mean <- mean(Nind)
    N_demes_sd <- sd(Nind)
    chartLabel <- locusName
    textSecond <- paste(paste("demes=",length(XCordinates),sep=""),
                        paste("N=",TotalSamples,sep=""),sep=", ")
    pointCol <- as.vector(clineFit[1,"pointCol"])
    curveCol <- as.vector(clineFit[1,"curveCol"])
    rep <- as.vector(clineFit[1,"rep"])
    iters <- as.vector(clineFit[1,"iters"])
  }
  if (curveFit==F & as.numeric(as.vector(clineFit[1,"N_total"]!=0))) {
    textSecond <- paste(textSecond,paste("dP=",deltaP,sep=""),sep=", ")
    plot(XCordinates,p_adj, yaxt="n",xaxt="n",type='n',bty="n",
         xlab="",ylab="",main=paste(chartLabel,textSecond,sep="\n"),ylim=c(-0.14,1.05),xlim=c(-640,25000),
         cex.main=cexMain,cex.axis=cexAxisPlot,cex.lab=cexAxisPlot)
    axis(2, at=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), labels=NA,col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=axisLineWdth,pos=-500, tck = -.01)
    axis(2, at=c(0,0.5,1.0), labels=c(0,0.5,1.0),col.axis="black", las=2.2,cex.axis=cexAxisPlot,lwd=0,pos=-820)
    axis(2, at=c(0,0.5,1.0), labels=NA,col.axis="black", las=2.2,cex.axis=cexAxisPlot,lwd=axisLineWdth*1.2,pos=-500,tck = -.018)
    axis(1, at=seq(-500,25000,500), labels=NA,col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=axisLineWdth,pos=-0.01, tck = -.01)
    axis(1, at=seq(0,25000,2000), labels=((seq(0,25000,2000)/1000)),col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=0,pos=-0.045)
    axis(1, at=seq(0,25000,1000), labels=NA,col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=axisLineWdth*1.2,pos=-0.01, tck = -.016)
    points(cbind(XCordinates,p_adj),cex=cexPoints,pch=21,col='black',bg=pointCol)
    mtext("geographic distance (km)",side=1,line=0,cex=cexAxisPlot)
    mtext("allele frequency (p)",side=2,line=2.9,cex=cexAxisPlot)
    lines(x=c(25000,25000),y=c(-0.01,1.02),lwd=axisLineWdth*1.2)
    lines(x=c(-500,25000),y=c(1.02,1.02),lwd=axisLineWdth*1.2)
    lines(x=c(-500,-500),y=c(1.00,1.02),lwd=axisLineWdth*1.2)
    lines(x=c(-500,-500),y=c(-0.01,0),lwd=axisLineWdth*1.2)
  }
  if (curveFit==T & as.numeric(as.vector(clineFit[1,"N_total"]!=0))) {
    textSecond <- paste(paste(textSecond,paste("dP=",deltaP,sep=""),sep=", "),clineModel,sep=", ")
    if (round(as.numeric(as.vector(clineFit[1,"centre_accept"])),0)!=0) {
      centre <- round(as.numeric(as.vector(clineFit[1,"centre_accept"])),0)
      centreL <- round(as.numeric(as.vector(clineFit[1,"centre_L95"])),0)
      centreU <- round(as.numeric(as.vector(clineFit[1,"centre_U95"])),0)
      width <- round(as.numeric(as.vector(clineFit[1,"width_accept"])),0)
      widthL <- round(as.numeric(as.vector(clineFit[1,"width_L95"])),0)
      widthU <- round(as.numeric(as.vector(clineFit[1,"width_U95"])),0)
      centreLabel <- paste(paste("c=",centre,sep=""),paste(paste("(",round(centreL,0),sep=""),
                                                           paste("-",round(centreU,0),sep="") ,")",sep=""),sep=" ")
      widthLabel <- paste(paste("w=",width,sep=""),paste(paste("(",round(widthL,0),sep=""),
                                                         paste("-",round(widthU,0),sep="") ,")",sep=""),sep=" ")
      textThird <- paste(centreLabel,widthLabel,sep=", ")
      pL <- round(as.numeric(as.vector(clineFit[1,"pL_accept"])),2)
      pLL <- round(as.numeric(as.vector(clineFit[1,"pL_L95"])),2)
      pLU <- round(as.numeric(as.vector(clineFit[1,"pL_U95"])),2)
      pLLabel <- paste(paste("p0=",pL,sep=""),paste(paste("(",round(pLL,2),sep=""),
                                                    paste("-",round(pLU,2),sep="") ,")",sep=""),sep=" ")
      pR <- round(as.numeric(as.vector(clineFit[1,"pR_accept"])),2)
      pRL <- round(as.numeric(as.vector(clineFit[1,"pR_L95"])),2)
      pRU <- round(as.numeric(as.vector(clineFit[1,"pR_U95"])),2)
      pRLabel <- paste(paste("p1=",pR,sep=""),paste(paste("(",round(pRL,2),sep=""),
                                                    paste("-",round(pRU,2),sep="") ,")",sep=""),sep=" ")
      fst <- round(as.numeric(as.vector(clineFit[1,"Fst_accept"])),3)
      fstL <- round(as.numeric(as.vector(clineFit[1,"Fst_L95"])),3)
      fstU <- round(as.numeric(as.vector(clineFit[1,"Fst_U95"])),3)
      fstLabel <- paste(paste("Fst=",fst,sep=""),paste(paste("(",round(fstL,3),sep=""),
                                                       paste("-",round(fstL,3),sep="") ,")",sep=""),sep=" ")
      textFourth <- paste(pLLabel,pRLabel,sep=", ")
      textFourth <- paste(textFourth,fstLabel,sep=", ")
      
      secondLevelHeader <- paste(textSecond,textThird,sep="\n")
      secondLevelHeader <- paste(secondLevelHeader,textFourth,sep="\n")
      # run details
      chartLabel <- paste(chartLabel,paste(paste("iter=",iters,sep=""),paste("rep ",rep,sep=""),sep=","))
      if (HeaderDetails==T) { 
        plot(XCordinates,p_adj, yaxt="n",xaxt="n",type='n',bty="n",
             xlab="",ylab="",main=paste(chartLabel,secondLevelHeader,sep="\n"),ylim=c(-0.14,1.05),xlim=c(-600,25000),
             cex.main=cexMain,cex.axis=cexAxisPlot,cex.lab=cexAxisPlot)
      }
      if (HeaderDetails==F & is.na(FigLabel)) {
        plot(XCordinates,p_adj, yaxt="n",xaxt="n",type='n',bty="n",
             xlab="",ylab="",main="",ylim=c(-0.14,1.05),xlim=c(-600,25000),
             cex.main=cexMain,cex.axis=cexAxisPlot,cex.lab=cexAxisPlot)
      }
      if (HeaderDetails==F & !is.na(FigLabel)) {
        plot(XCordinates,p_adj, yaxt="n",xaxt="n",type='n',bty="n",
             xlab="",ylab="",main="",ylim=c(-0.14,1.05),xlim=c(-600,25000),
             cex.main=cexMain,cex.axis=cexAxisPlot,cex.lab=cexAxisPlot)
        mtext(FigLabel, side=3, line=1, adj=0, cex=cexMain)
      }
      axis(2, at=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), labels=NA,col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=axisLineWdth,pos=-500, tck = -.01)
      axis(2, at=c(0,0.5,1.0), labels=c(0,0.5,1.0),col.axis="black", las=2.2,cex.axis=cexAxisPlot,lwd=0,pos=-820)
      axis(2, at=c(0,0.5,1.0), labels=NA,col.axis="black", las=2.2,cex.axis=cexAxisPlot,lwd=axisLineWdth*1.2,pos=-500,tck = -.018)
      axis(1, at=seq(-500,25000,500), labels=NA,col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=axisLineWdth,pos=-0.01, tck = -.01)
      axis(1, at=seq(0,25000,2000), labels=((seq(0,25000,2000)/1000)),col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=0,pos=-0.045)
      axis(1, at=seq(0,25000,1000), labels=NA,col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=axisLineWdth*1.2,pos=-0.01, tck = -.016)
      curve(cline_model_sig_symm_plot, 0, 24500, n=2001, add=TRUE, lwd=lwdCurve,lty=1, col=curveCol)
      if (!is.na(refLocus) & locusName != refLocus) {
        thisDrop_ref <- which(clineFit[,"locus"]==refLocus)
        if (as.numeric(as.vector(clineFit[thisDrop,"N_total"]))!=0) {
          TotalSamples_ref <- as.numeric(as.vector(clineFit[thisDrop_ref,"N_total"]))
          XCordinate_ref <- as.numeric(as.vector(unlist(strsplit(as.vector(clineFit[thisDrop_ref,"Xvals"]),";"))))
          p_adj_ref <- as.numeric(as.vector(unlist(strsplit(as.vector(clineFit[thisDrop_ref,"p"]),";"))))
          Ngenes_ref <- as.numeric(as.vector(unlist(strsplit(as.vector(clineFit[thisDrop_ref,"Ngenes"]),";"))))
          Nind_ref <- Ngenes/2
          centre_ref <- round(as.numeric(as.vector(clineFit[thisDrop_ref,"centre_accept"])),0)
          width_ref <- round(as.numeric(as.vector(clineFit[thisDrop_ref,"width_accept"])),0)
          pL_ref <- round(as.numeric(as.vector(clineFit[thisDrop_ref,"pL_accept"])),2)
          pR_ref <- round(as.numeric(as.vector(clineFit[thisDrop_ref,"pR_accept"])),2)
          curve(cline_model_sig_symm_plot_ref, 0, 24500, n=2001, add=TRUE, lwd=lwdCurve*1.3,lty=1, col="red")
          if (refPoints==T) {
            points(cbind(XCordinate_ref,p_adj_ref),cex=cexPoints*0.8,pch=22,col='black',bg="red")
          }
        }
      }
      
      if (!is.na(HaplotypeData)) {
        #str(HaplotypeData)
        #colnames(HaplotypeData[["Fla_2loc"]][["DemeMatrix2D_combine_noZero"]])
        thisHapData <- HaplotypeData[["Fla_2loc"]][["DemeMatrix2D_combine_noZero"]]
        XCordinate_hap <- as.numeric(as.vector(thisHapData[,"distAlongTransect"]))
        r_freq <- as.numeric(as.vector(thisHapData[,"r"]))
        points(cbind(XCordinate_hap,r_freq),cex=cexPoints,pch=21,col='black',bg="blue")
      }
      points(cbind(XCordinates,p_adj),cex=cexPoints,pch=21,col='black',bg=pointCol)
      mtext("geographic distance (km)",side=1,line=0,cex=cexAxisPlot)
      mtext("allele frequency (p)",side=2,line=2.9,cex=cexAxisPlot)
      lines(x=c(25000,25000),y=c(-0.01,1.02),lwd=axisLineWdth*1.2)
      lines(x=c(-500,25000),y=c(1.02,1.02),lwd=axisLineWdth*1.2)
      lines(x=c(-500,-500),y=c(1.00,1.02),lwd=axisLineWdth*1.2)
      lines(x=c(-500,-500),y=c(-0.01,0),lwd=axisLineWdth*1.2)
    }
  }
  
  
}
