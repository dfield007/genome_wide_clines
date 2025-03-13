# Individual or multiple loci on one Cline plot
# David Field
# 1/4/2024
clinePlot <- function(clineModel,modelSummaryInput,curveFit,cexMain=1.1,cexAxisPlot= 1.85,
            axisLineWdth=2.5,cexPoints=1.8,lwdCurve=4,refLocus=NA,refPoints=F,HeaderDetails=F,FigLabel=NA,multipleCurves=F,
            multiplePoints=F,HaplotypeData=NA,poolSeqClines,poolSeqStats,colourTable) {
  # David Field: 1/4/2024
  # this plot will draw any loci present in the modelSummaryInput table. If only one row,  
  
  # clineModel="symm_5p"; curveFit=T; modelSummaryInput = clineFit_scan_s261_720757[4,]; curveFit <- T; multipleCurves =T
  # HaplotypeData <- HZ_Planoles_150mDemes_FLAr; # HeaderDetails <- T; FigLabel<-NA; multiplePoints=T
  # cexMain = 1.1; cexAxisPlot=1.85; axisLineWdth=2.5; cexPoints=1.8;lwdCurve=3;refLocus=NA;refPoints=F; 
  # poolSeqClines <- allChr_Clines_Chr5clines ; poolSeqStats<- poolSeqStatsKeep; <- colourChoice 
  # refLocus <- F; HaplotypeData <- NA
  cline_model_sig_symm_plot <- function(XCordinate) {
    ((pL + ((pR-pL)/(1 + exp((-4*(XCordinate-centre))/width)))))
  }
  cline_model_sig_symm_plot_ref <- function(XCordinate_ref) {
    ((pL_ref + ((pR_ref-pL_ref)/(1 + exp((-4*(XCordinate_ref-centre_ref))/width_ref)))))
  }
  
  for (thisRow in 1:nrow(modelSummaryInput)) {
    # thisRow <- 1
    thisLocus <- modelSummaryInput[thisRow,"locus"]
    if (as.numeric(as.vector(modelSummaryInput[thisRow,"N_total"]))==0) {
      next
    }
    if (thisRow==1) {
      if (HeaderDetails==T) {
        par(mfrow=c(1,1),mar=c(3.7,4.5,4.8,2),oma=c(1,1,1,1)) # (bottom, left, top, right)
      }
      if (HeaderDetails==F) {
        par(mfrow=c(1,1),mar=c(4.5,2,1,1),oma=c(2,2.5,2,1)) # (bottom, left, top, right)
      }
    }
    if (as.numeric(as.vector(modelSummaryInput[thisRow,"N_total"]))!=0) {
      deltaP <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"delta_p"])),2)
      TotalSamples <- as.numeric(as.vector(modelSummaryInput[thisRow,"N_total"]))
      XCordinates <- as.numeric(as.vector(unlist(strsplit(as.vector(modelSummaryInput[thisRow,"Xvals"]),";"))))
      p_adj <- as.numeric(as.vector(unlist(strsplit(as.vector(modelSummaryInput[thisRow,"p"]),";"))))
      Ngenes <- as.numeric(as.vector(unlist(strsplit(as.vector(modelSummaryInput[thisRow,"Ngenes"]),";"))))
      Nind <- sum(Ngenes)/2
      N_demes_mean <- mean(Nind)
      N_demes_sd <- sd(Nind)
      chartLabel <- thisLocus
      textSecond <- paste(paste("demes=",length(XCordinates),sep=""),
                          paste("N=",Nind,sep=""),sep=", ")
      pointCol <- as.vector(modelSummaryInput[thisRow,"pointCol"])
      curveCol <- as.vector(modelSummaryInput[thisRow,"curveCol"])
      rep <- as.vector(modelSummaryInput[thisRow,"rep"])
      iters <- as.vector(modelSummaryInput[thisRow,"iters"])
    }
    
    if (curveFit==F & as.numeric(as.vector(modelSummaryInput[thisRow,"N_total"]))!=0) {
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
    if (curveFit==T & as.numeric(as.vector(modelSummaryInput[thisRow,"N_total"]))!=0) {
      textSecond <- paste(paste(textSecond,paste("dP=",deltaP,sep=""),sep=", "),clineModel,sep=", ")
      if (round(as.numeric(as.vector(modelSummaryInput[thisRow,"centre_accept"])),0)!=0) {
        centre <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"centre_accept"])),0)
        centreL <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"centre_L95"])),0)
        centreU <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"centre_U95"])),0)
        width <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"width_accept"])),0)
        widthL <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"width_L95"])),0)
        widthU <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"width_U95"])),0)
        centreLabel <- paste(paste("c=",round(centre/1000,2),sep=""),paste(paste("(",round(centreL/1000,2),sep=""),
                                                                           paste("-",round(centreU/1000,2),sep="") ,")",sep=""),sep=" ")
        widthLabel <- paste(paste("w=",round(width/1000,2),sep=""),paste(paste("(",round(widthL/1000,2),sep=""),
                                                                         paste("-",round(widthU/1000,2),sep="") ,")",sep=""),sep=" ")
        textThird <- paste(centreLabel,widthLabel,sep=", ")
        pL <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"pL_accept"])),2)
        pLL <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"pL_L95"])),2)
        pLU <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"pL_U95"])),2)
        pLLabel <- paste(paste("p0=",pL,sep=""),paste(paste("(",round(pLL,2),sep=""),
                                                      paste("-",round(pLU,2),sep="") ,")",sep=""),sep=" ")
        pR <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"pR_accept"])),2)
        pRL <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"pR_L95"])),2)
        pRU <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"pR_U95"])),2)
        pRLabel <- paste(paste("p1=",pR,sep=""),paste(paste("(",round(pRL,2),sep=""),
                                                      paste("-",round(pRU,2),sep="") ,")",sep=""),sep=" ")
        fst <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"Fst_accept"])),3)
        fstL <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"Fst_L95"])),3)
        fstU <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"Fst_U95"])),3)
        fstLabel <- paste(paste("Fst=",fst,sep=""),paste(paste("(",round(fstL,3),sep=""),
                                                         paste("-",round(fstL,3),sep="") ,")",sep=""),sep=" ")
        textFourth <- paste(pLLabel,pRLabel,sep=", ")
        textFourth <- paste(textFourth,fstLabel,sep=", ")
        
        secondLevelHeader <- paste(textSecond,textThird,sep="\n")
        secondLevelHeader <- paste(secondLevelHeader,textFourth,sep="\n")
        # run details
        chartLabel <- paste(thisLocus,paste(paste(", iter=",iters,sep=""),paste("rep ",rep,sep=""),sep=","))
        if (HeaderDetails==T & thisRow==1) { 
          plot(XCordinates,p_adj, yaxt="n",xaxt="n",type='n',bty="n",
               xlab="",ylab="",main=paste(chartLabel,secondLevelHeader,sep="\n"),ylim=c(-0.14,1.05),xlim=c(-600,25000),
               cex.main=cexMain,cex.axis=cexAxisPlot,cex.lab=cexAxisPlot)
          axis(2, at=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), labels=NA,col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=axisLineWdth,pos=-500, tck = -.01)
          axis(2, at=c(0,0.5,1.0), labels=c(0,0.5,1.0),col.axis="black", las=2.2,cex.axis=cexAxisPlot,lwd=0,pos=-820)
          axis(2, at=c(0,0.5,1.0), labels=NA,col.axis="black", las=2.2,cex.axis=cexAxisPlot,lwd=axisLineWdth*1.2,pos=-500,tck = -.018)
          axis(1, at=seq(-500,25000,500), labels=NA,col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=axisLineWdth,pos=-0.01, tck = -.01)
          axis(1, at=seq(0,25000,2000), labels=((seq(0,25000,2000)/1000)),col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=0,pos=-0.08)
          axis(1, at=seq(0,25000,1000), labels=NA,col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=axisLineWdth*1.2,pos=-0.01, tck = -.016)
          
          if (length(poolSeqClines)>0) {
            for (thisRowPoolSeq in 1:nrow(poolSeqClines)) {
              # thisRowPoolSeq <- 1
              poolSeqClineThisLocus <- cbind(poolSeqStats[2:7,"pool"],poolSeqStats[2:7,"DistAlongTransect"],as.vector(unlist(poolSeqClines[thisRowPoolSeq,9:14])))
              colnames(poolSeqClineThisLocus) <- c("pool","DistAlongTransect","p")
              poolSeqClineThisLocus <- as.data.frame(poolSeqClineThisLocus)
              poolSeqClineThisLocus[,"DistAlongTransect"] <- as.numeric(poolSeqClineThisLocus[,"DistAlongTransect"])
              poolSeqClineThisLocus[,"p"] <- as.numeric(poolSeqClineThisLocus[,"p"])
              poolSeqClineThisLocus <- poolSeqClineThisLocus[order(poolSeqClineThisLocus[,"DistAlongTransect"]),]
              
              # background loci not linked to any colour loci
              if (poolSeqClines[thisRowPoolSeq,"ColourGene"]==-9) {
                # lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=colourTable[colourTable[,"Region"]=="background","lineColour"])
                back_col <- colourTable[colourTable[,"Region"]=="background","lineColour"]
                back_col_transparent <- rgb(128, 128, 128, maxColorValue = 255, alpha = 0.3)
                lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=5,col=rgb(180, 180, 180, maxColorValue = 255, alpha = 1))
                points(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],pch=16,lwd=5,col=rgb(180, 180, 180, maxColorValue = 255, alpha = 1))
                
                }
              # ROS EL
              
              if (poolSeqClines[thisRow,"ColourGene"]=="ROS1" | poolSeqClines[thisRow,"ColourGene"]=="ROS2" | poolSeqClines[thisRow,"ColourGene"]=="ROS3") {
                lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=colourTable[colourTable[,"Region"]=="ROS","lineColour"])
              }
              if (poolSeqClines[thisRow,"ColourGene"]=="ELUTA") {
                
                
                lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=colourTable[colourTable[,"Region"]=="ELUTA","lineColour"])
              }
              if ((poolSeqClines[thisRow,"ColourGene"]=="ROSEA" | poolSeqClines[thisRow,"ColourGene"]=="ELUTA") & poolSeqClines[thisRow,"distanceToGene"]!="inGene") {
                
                
                lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=colourTable[colourTable[,"Region"]=="ROSelLinked","lineColour"])
              }
              # RUBIA
              if (poolSeqClines[thisRow,"ColourGene"]=="RUBIA" & poolSeqClines[thisRow,"distanceToGene"]!="inGene") {
                rubCol_linked <- colourTable[colourTable[,"Region"]=="RUBIA","lineColour"]
                rubCol_linked_trans <- rgb(228, 26, 28, maxColorValue = 255, alpha = 0.3)
                lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=0.8,col=rubCol_linked_trans)
              }
              if (poolSeqClines[thisRow,"ColourGene"]=="RUBIA" & poolSeqClines[thisRow,"distanceToGene"]=="inGene") {
                rubCol <- colourTable[colourTable[,"Region"]=="RUBIA","lineColour"]
                rubCol_trans <- rgb(228, 26, 28, maxColorValue = 255, alpha = 0.7)
                lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.3,col=rubCol_trans)
              }
              # SULF
              if (poolSeqClines[thisRow,"ColourGene"]=="SULF") {
                lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=colourTable[colourTable[,"Region"]=="SULFlinked","lineColour"])
              }
              # CREMOSA
              # no loci inGene for CREMOSA, therefore only linked and use one colour
              if (poolSeqClines[thisRow,"ColourGene"]=="CREMOSA" & poolSeqClines[thisRow,"distanceToGene"]!="inGene") {
                lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=colourTable[colourTable[,"Region"]=="CREMOSAlinked","lineColour"])
              }
              if (poolSeqClines[thisRow,"ColourGene"]=="CREMOSA" & poolSeqClines[thisRow,"distanceToGene"]=="inGene") {
                lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=colourTable[colourTable[,"Region"]=="CREMOSA","lineColour"])
              }
              
              # AURINA
              if (poolSeqClines[thisRow,"ColourGene"]=="AURINA") {
                lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=colourTable[colourTable[,"Region"]=="AURINA","lineColour"])
              }
              # AURINA
              if (poolSeqClines[thisRow,"ColourGene"]=="AURINA" & poolSeqClines[thisRow,"distanceToGene"]!="inGene") {
                lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=colourTable[colourTable[,"Region"]=="AURINAlinked","lineColour"])
              }
              # FLAVIA
              if (poolSeqClines[thisRow,"ColourGene"]=="FLAVIA" & poolSeqClines[thisRow,"distanceToGene"]!="inGene") {
                lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=0.8,col=colourTable[colourTable[,"Region"]=="FLAVIAlinked","lineColour"])
              }
              if (poolSeqClines[thisRow,"ColourGene"]=="FLAVIA" & poolSeqClines[thisRow,"distanceToGene"]=="inGene") {
                lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.4,col=colourTable[colourTable[,"Region"]=="FLAVIA","lineColour"])
              }
            }
          }
        }
        if (HeaderDetails==F & is.na(FigLabel)) {
          plot(XCordinates,p_adj, yaxt="n",xaxt="n",type='n',bty="n",
               xlab="",ylab="",main="",ylim=c(-0.14,1.05),xlim=c(-600,25000),
               cex.main=cexMain,cex.axis=cexAxisPlot,cex.lab=cexAxisPlot)
          axis(2, at=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), labels=NA,col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=axisLineWdth,pos=-500, tck = -.01)
          axis(2, at=c(0,0.5,1.0), labels=c(0,0.5,1.0),col.axis="black", las=2.2,cex.axis=cexAxisPlot,lwd=0,pos=-820)
          axis(2, at=c(0,0.5,1.0), labels=NA,col.axis="black", las=2.2,cex.axis=cexAxisPlot,lwd=axisLineWdth*1.2,pos=-500,tck = -.018)
          axis(1, at=seq(-500,25000,500), labels=NA,col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=axisLineWdth,pos=-0.01, tck = -.01)
          axis(1, at=seq(0,25000,2000), labels=((seq(0,25000,2000)/1000)),col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=0,pos=-0.08)
          axis(1, at=seq(0,25000,1000), labels=NA,col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=axisLineWdth*1.2,pos=-0.01, tck = -.016)
          
        }
        if (HeaderDetails==F & !is.na(FigLabel)) {
          plot(XCordinates,p_adj, yaxt="n",xaxt="n",type='n',bty="n",
               xlab="",ylab="",main="",ylim=c(-0.14,1.05),xlim=c(-600,25000),
               cex.main=cexMain,cex.axis=cexAxisPlot,cex.lab=cexAxisPlot)
          mtext(FigLabel, side=3, line=1, adj=0, cex=cexMain)
          axis(2, at=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), labels=NA,col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=axisLineWdth,pos=-500, tck = -.01)
          axis(2, at=c(0,0.5,1.0), labels=c(0,0.5,1.0),col.axis="black", las=2.2,cex.axis=cexAxisPlot,lwd=0,pos=-820)
          axis(2, at=c(0,0.5,1.0), labels=NA,col.axis="black", las=2.2,cex.axis=cexAxisPlot,lwd=axisLineWdth*1.2,pos=-500,tck = -.018)
          axis(1, at=seq(-500,25000,500), labels=NA,col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=axisLineWdth,pos=-0.01, tck = -.01)
          axis(1, at=seq(0,25000,2000), labels=((seq(0,25000,2000)/1000)),col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=0,pos=-0.08)
          axis(1, at=seq(0,25000,1000), labels=NA,col.axis="black", las=1,cex.axis=cexAxisPlot,lwd=axisLineWdth*1.2,pos=-0.01, tck = -.016)
          
        }
         
        if (multipleCurves==F) {
          curve(cline_model_sig_symm_plot, 0, 24500, n=2001, add=TRUE, lwd=lwdCurve,lty=1, col=curveCol)
        }
        if (multipleCurves==T) {
          if (nrow(modelSummaryInput)>1) {
            for (thisRow in 1:nrow(modelSummaryInput)) {
              # thisRow <-2
              XCordinates <- as.numeric(as.vector(unlist(strsplit(as.vector(modelSummaryInput[thisRow,"Xvals"]),";"))))
              p_adj <- as.numeric(as.vector(unlist(strsplit(as.vector(modelSummaryInput[thisRow,"p"]),";"))))
              pointCol <- as.vector(modelSummaryInput[thisRow,"pointCol"])
              curveCol <- as.vector(modelSummaryInput[thisRow,"curveCol"])
              pL <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"pL_accept"])),2)
              pR <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"pR_accept"])),2)
              centre <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"centre_accept"])),0)
              width <- round(as.numeric(as.vector(modelSummaryInput[thisRow,"width_accept"])),0)
              if (modelSummaryInput[thisRow,"flipCline"]==1) {
                p_adj <- 1-p_adj
                pL_temp <-1-pL
                pR_temp <-1-pR
                rm(pR); rm(pL)
                pR <- pR_temp
                pL <- pL_temp
              }
              curve(cline_model_sig_symm_plot, 0, 24500, n=2001, add=TRUE, lwd=lwdCurve,lty=1, col=curveCol)
            }
           }
          }
        if (refLocus!=F & thisLocus != refLocus) {
          thisDrop_ref <- which(modelSummaryInput[,"locus"]==refLocus)
          if (as.numeric(as.vector(modelSummaryInput[1,"N_total"]))!=0) {
            TotalSamples_ref <- as.numeric(as.vector(modelSummaryInput[thisDrop_ref,"N_total"]))
            XCordinate_ref <- as.numeric(as.vector(unlist(strsplit(as.vector(modelSummaryInput[thisDrop_ref,"Xvals"]),";"))))
            p_adj_ref <- as.numeric(as.vector(unlist(strsplit(as.vector(modelSummaryInput[thisDrop_ref,"p"]),";"))))
            Ngenes_ref <- as.numeric(as.vector(unlist(strsplit(as.vector(modelSummaryInput[thisDrop_ref,"Ngenes"]),";"))))
            Nind_ref <- Ngenes/2
            centre_ref <- round(as.numeric(as.vector(modelSummaryInput[thisDrop_ref,"centre_accept"])),0)
            width_ref <- round(as.numeric(as.vector(modelSummaryInput[thisDrop_ref,"width_accept"])),0)
            pL_ref <- round(as.numeric(as.vector(modelSummaryInput[thisDrop_ref,"pL_accept"])),2)
            pR_ref <- round(as.numeric(as.vector(modelSummaryInput[thisDrop_ref,"pR_accept"])),2)
            curve(cline_model_sig_symm_plot_ref, 0, 24500, n=2001, add=TRUE, lwd=lwdCurve*1.3,lty=1, col="red")
            if (refPoints==T) {
              points(cbind(XCordinate_ref,p_adj_ref),cex=cexPoints*0.8,pch=22,col='black',bg="red")
            }
          }
        }
       
        
        if (multiplePoints==T) {
          if (nrow(modelSummaryInput)>1) {
            for (thisRow in 1:nrow(modelSummaryInput)) {
              XCordinates <- as.numeric(as.vector(unlist(strsplit(as.vector(modelSummaryInput[thisRow,"Xvals"]),";"))))
              p_adj <- as.numeric(as.vector(unlist(strsplit(as.vector(modelSummaryInput[thisRow,"p"]),";"))))
              pointCol <- as.vector(modelSummaryInput[thisRow,"pointCol"])
              if (modelSummaryInput[thisRow,"flipCline"]==1) {
                p_adj <- 1-p_adj
                pR_temp <-1-pL
                pL_temp <-1-pR
                pR <- pR_temp
                pL <- pL_temp
              }
              points(cbind(XCordinates,p_adj),cex=cexPoints,pch=21,col='black',bg=pointCol)
            }
          }
        }
        if (multiplePoints==F) {
              XCordinates <- as.numeric(as.vector(unlist(strsplit(as.vector(modelSummaryInput[1,"Xvals"]),";"))))
              p_adj <- as.numeric(as.vector(unlist(strsplit(as.vector(modelSummaryInput[1,"p"]),";"))))
              pointCol <- as.vector(modelSummaryInput[1,"pointCol"])
              pointCol_transparent <- adjustcolor(pointCol, alpha.f = 0.8)
              points(cbind(XCordinates,p_adj),cex=cexPoints,pch=21,col=rgb(128, 128, 128, maxColorValue = 255, alpha = 0.5),bg=rgb(128, 128, 128, maxColorValue = 255, alpha = 0.5))
        }
        if (!is.na(HaplotypeData)) {
          #str(HaplotypeData)
          #colnames(HaplotypeData[["Fla_2loc"]][["DemeMatrix2D_combine_noZero"]])
          thisHapData <- HaplotypeData[["Fla_2loc"]][["DemeMatrix2D_combine_noZero"]]
          XCordinate_hap <- as.numeric(as.vector(thisHapData[,"distAlongTransect"]))
          r_freq <- as.numeric(as.vector(thisHapData[,"r"]))
          points(cbind(XCordinate_hap,r_freq),cex=cexPoints,pch=21,col='black',bg="blue")
          
          # fit a normal distribution curve here....
          #fit <- fitdistr(r_freq, "normal")
          # linear_model1 <- lm(y~x, data=sample_data) 
          # linear_model2 <- lm(r_freq~poly(XCordinate_hap,2,raw=TRUE)) 
          # linear_model3 <- lm(r_freq~poly(XCordinate_hap,3,raw=TRUE)) 
          # linear_model4 <- lm(r_freq~poly(XCordinate_hap,4,raw=TRUE)) 
          # linear_model5 <- lm(r_freq~poly(XCordinate_hap,5,raw=TRUE)) 
          # lines(XCordinate_hap, predict(linear_model2), col='red') 
          # lines(XCordinate_hap, predict(linear_model3), col='purple') 
          # lines(XCordinate_hap, predict(linear_model4), col='blue') 
          # 
          # lines(XCordinate_hap, predict(linear_model3, data.frame(x=x_axis)), col='purple') 
          # lines(XCordinate_hap, predict(linear_model4, data.frame(x=x_axis)), col='blue') 
          # points(bezierCurve(XCordinate_hap,r_freq,50), type="l", col="pink")
          # 
          thisData <- as.data.frame(cbind(XCordinate_hap,r_freq))
          
          loessFit <- loess(r_freq~XCordinate_hap, thisData, span = 0.9)
          loessFit <- data.frame(x=loessFit$x,y=loessFit$fitted)
          loessFit <- loessFit[order(loessFit$XCordinate_hap),]
          # turn curve off
          #lines(x=loessFit,col='blue',lwd=lwdCurve)
          
          #lowessFit <-lowess((thisData),f = .9,iter=1)
          #lines(lowessFit,col='red')
          
        }
        #points(cbind(XCordinates,p_adj),cex=cexPoints,pch=21,col='black',bg=pointCol)
        mtext("geographic distance (km)",side=1,line=0,cex=cexAxisPlot)
        mtext("allele frequency (p)",side=2,line=2.9,cex=cexAxisPlot)
        lines(x=c(25000,25000),y=c(-0.01,1.02),lwd=axisLineWdth*1.2)
        lines(x=c(-500,25000),y=c(1.02,1.02),lwd=axisLineWdth*1.2)
        lines(x=c(-500,-500),y=c(1.00,1.02),lwd=axisLineWdth*1.2)
        lines(x=c(-500,-500),y=c(-0.01,0),lwd=axisLineWdth*1.2)
      }
    }
  }
}
