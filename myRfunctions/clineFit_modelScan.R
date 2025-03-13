# Cline fit optimal transect
# For repeating cline fitting over a range of transect orientations
# David Field 14/2/2024

clineFit_modelScan <- function(inData,modelSummary,thisTrial,thisPop,dpThresh,
                              removePhenoMiss,clineModel="symm_5p",adjustXY,pointsCol,gridAdd,
                              plotAllele1D,plotDem,addPie,SNPlociListUpdated,plotAxisLab=F,demeNmin,demeNmax,
                              returnFreq,startParams,fileFormat="png",plotMarkDownAllDemes=F,
                              plotMarkDownZoomed=F,plotEastNorthing=T,tracePlots,plotClines=T,multipleCurves=F,
                              multiplePoints=F) {
   # inData <- GenoEcol_Planoles; thisTrial <- "OptimalTransect"; clineModel="symm_5p"; modelSummary <- lnLLsummary
   # demeNmax =50; demeNmin <- 10; thisPop <- "Planoles"; dpThresh <- 0.6 ;removePhenoMiss=F;adjustXY=F
   # pointsCol=F; gridAdd=F; plotAllele1D=T;plotDem=c(3,3);addPie=T;plotAxisLab=F;
   # returnFreq=1;startParams=c(10000,3000,0.1,0.9,15);fileFormat="png"; tracePlots=T;plotClines=T
    # modelSummary <- clineFit_scan_s91_39699
  # lociNames <- c("s1187_290152", # LG1
  #                "s48_16807","s154_504353", # LG 2 - 41.985
  #                "s316_93292", "s316_218136","s316_219369","s316_257789", # LG - 43.175
  #                "s60_358", # LG2 - 43.175
  #                "s1140_224946", # LG2 - 44.907
  #                "s91_39699","s91_78805", # LG 4 - sulf
  #                "s261_720757", # LG 5 - FLS
  #                "ros_assembly_543443","ros_assembly_715015", # ros & el
  #                "RosInferred", # phenotype
  #                "scpDNA_seq_58467") # cp genome
  # thisLocus <- "s91_39699"
  # thisLocus <- "s1187_290152"
  # thisLocus <- "ros_assembly_543443"
  # thisLocus <- "s261_720757"
  
   # New 2024
  
  for (thisRow in 1:nrow(modelSummary)) {
    # thisRow <- 6
    # note, some parameters coming in from modelSummary table in new version 2024
    thisLocus <- as.vector(modelSummary[thisRow,"locus"])
    locusName <- thisLocus
    demeSize <- as.numeric(as.vector(modelSummary[thisRow,"demeSize"]))
    rep <- as.vector(modelSummary[thisRow,"rep"])
    iterations <- as.numeric(as.vector(modelSummary[thisRow,"iters"]))
    burnin <- as.numeric(as.vector(modelSummary[thisRow,"burnin"]))
    thisGradient <- as.numeric(as.vector(modelSummary[thisRow,"gradient"]))
    thisIntercept <- as.numeric(as.vector(modelSummary[thisRow,"intercept"]))
    if (thisRow==1) {
      DemeData <- DemeSetup2D(HZdata=inData,lociNames=thisLocus,demeSize=demeSize,minSample=demeNmin,
                              removePhenoMiss=F,transectGradient = thisGradient, transectIntercept = thisIntercept,  
                              adjustXY=F, pointsCol = F, gridAdd=F,plotAllele1D = T, plotDem=c(3,3), 
                              addPie=T,SNPlociListUpdated,verbose=F,plotAxisLab=F,plotPerp=T,plotType="pdf",
                              poolSeqPlot=F,plotMap=F,noPlots=F,plotMarkDownAllDemes=F,plotEastNorthing=T,
                              plotMarkDownZoomed=F,plotTerrain=F,plotKASPrep=T,FigName=paste0("_ModelRow_",thisRow),plotCentreLine=T,
                              poolSeqClines=allChr_Clines_Chr1clines,poolSeqStats=poolLoc,colourTable=colourChoice)
      thisDataLoc <- DemeData[[which(names(DemeData)==thisLocus)]]
    }
    if (thisRow>1) {
      if (thisLocus != modelSummary[thisRow-1,"locus"] | thisGradient!=as.vector(modelSummary[thisRow-1,"gradient"])) {
        rm(DemeData)
        DemeData <- DemeSetup2D(HZdata=inData,lociNames=thisLocus,demeSize=demeSize,minSample=demeNmin,
                                removePhenoMiss=F,transectGradient = thisGradient, transectIntercept = thisIntercept,  
                                adjustXY=F, pointsCol = F, gridAdd=F,plotAllele1D = T, plotDem=c(3,3), 
                                addPie=T,SNPlociListUpdated,verbose=F,plotAxisLab=F,plotPerp=T,plotType="pdf",
                                poolSeqPlot=T,plotMap=F,noPlots=F,plotMarkDownAllDemes=F,plotEastNorthing=T,
                                plotMarkDownZoomed=F,plotTerrain=F,plotKASPrep=T,FigName=paste0("_ModelRow_",thisRow),plotCentreLine=T,
                                poolSeqClines=allChr_Clines_Chr1clines,poolSeqStats=poolLoc,colourTable=colourChoice)
        thisDataLoc <- DemeData[[which(names(DemeData)==thisLocus)]]
      }
    }
    
    # NEW approach 2020 and updated again for 2024
    XCordinates <- thisDataLoc[["DemeMatrix2D_DistanceRefLine"]][,2]
    thisLoc_DemeMatrix2D_NoZero <- thisDataLoc[["DemeMatrix2D_combine_noZero"]]
    
    # mininum sample size
    sampleSize <- (thisLoc_DemeMatrix2D_NoZero[,"pCount"]+thisLoc_DemeMatrix2D_NoZero[,"qCount"])/2
    theseDemesGood <- which(sampleSize>=demeNmin)
    TotalSamples <- sum(sampleSize[sampleSize>=demeNmin])
    
    # SANN cline fit
    # gather values for demes
    Xvals <- XCordinates 
    pCount_raw <- thisLoc_DemeMatrix2D_NoZero$pCount
    qCount_raw <- thisLoc_DemeMatrix2D_NoZero$qCount
    Ngenes_raw <- pCount_raw + qCount_raw
    pCount <- pCount_raw
    qCount <- qCount_raw
    Ngenes <- Ngenes_raw
    obsVals <- cbind(pCount,qCount,Ngenes)
    
    # skip loci with weak allele frequency differences
    striatum_p <- mean(thisLoc_DemeMatrix2D_NoZero[1:3,"p"])
    pseudom_p <- mean(thisLoc_DemeMatrix2D_NoZero[(nrow(thisLoc_DemeMatrix2D_NoZero)-3):(nrow(thisLoc_DemeMatrix2D_NoZero)),"p"])
    deltaP <- abs(pseudom_p-striatum_p)
    #if (deltaP<0.6) {
    #  next
    #}
    
    thisLoc_Demes_good <- thisLoc_DemeMatrix2D_NoZero[theseDemesGood,]
    
    # allele frequencies
    alleleFreqs <- thisLoc_Demes_good[,"p"]
    
    # plot raw data
    par(mfrow=c(1,1),mar=c(5,5,4,2),oma=c(2,2,2,2)) # (bottom, left, top, right)
    #chartLabel <- markerLabelPieChart(locusName,SNPlociListUpdated)
    textSec <- paste(paste("demes=",length(theseDemesGood),sep=""),
                     paste("n=",TotalSamples,sep=""),sep=", ")
    plot(XCordinates,alleleFreqs, main=paste(locusName,textSec,sep="\n"), 
         cex.main=0.9,cex.axis=1.1,cex.lab=1.1,ylim=c(0,1), xlim=c(min(XCordinates)-500,max(XCordinates)+500), 
         type="n",pch=21,xlab = "distance (m)",ylab = "allele frequency",lwd=2)
    points(cbind(XCordinates,alleleFreqs),cex=1,pch=21,col='black',bg='grey')
    # gather values for demes
    pCount_raw <- as.numeric(as.vector(thisLoc_Demes_good[,"pCount"]))
    qCount_raw <- as.numeric(as.vector(thisLoc_Demes_good[,"qCount"]))
    Ngenes_raw <- pCount_raw + qCount_raw
    demeN <- (pCount_raw + qCount_raw)/2
    pCount <- pCount_raw
    qCount <- qCount_raw
    Ngenes <- Ngenes_raw
    p_raw <- pCount_raw/Ngenes_raw
    
    # reduce size of exceptionally large N for some demes (subsample) 
    # Note - otherwise computational problem with pochhammer function
    # demeMaxSize <- 90
    if (any(demeN > demeNmax)) {
      largeNDemes <- which(demeN > demeMaxSize)
      for (thisDeme in 1:length(largeNDemes)) {
        # thisDeme <- 2
        thisN <- Ngenes[largeNDemes[thisDeme]]
        thisP <- pCount[largeNDemes[thisDeme]]
        thisQ <- qCount[largeNDemes[thisDeme]]
        subsample <- sample(c(1,0),demeNmax,prob=c((thisP/thisN),(thisQ/thisN)),replace=T)
        pCount[largeNDemes[thisDeme]] <- sum(subsample==1)
        qCount[largeNDemes[thisDeme]] <- sum(subsample==0)
        Ngenes[largeNDemes[thisDeme]] <- sum(subsample==1) + sum(subsample==0)
      }
    }
    
    obsVals <- cbind(pCount,qCount,Ngenes)
    
    # # new adjusted allele frequency
    p_adj <- pCount/Ngenes
    # demeNmax of 50 looks better - highly correlated with original values
    # plot(p_raw,p_adj)

    modelSummary[thisRow,"numDemes"] <- length(pCount)
    modelSummary[thisRow,"N_total"] <- sum(Ngenes)/2
    modelSummary[thisRow,"Xvals"] <- paste(XCordinates,collapse=";")
    modelSummary[thisRow,"Ngenes"] <- paste(Ngenes,collapse=";")
    modelSummary[thisRow,"p"] <- paste(p_adj,collapse=";")
    modelSummary[thisRow,"delta_p"] <- deltaP
    
    # old routine from slowClines_RosEl 2018 Aug (which works)
    cat(c("..Cline fitting.."))
    SANNfit_thisLocus <- SANN_clinefit(fn = logLik_5par_SANN,obsVals, positions=XCordinates, 
                                       initialParam = startParams,scale=0.3,retvals=TRUE,nmax = iterations,
                                       verbose=TRUE,acceptscale = 1.05, rejectscale = (1/1.05), 
                                       retfreq=returnFreq)
    # SANNfit_thisLocus$retvals[which(SANNfit_thisLocus$retvals[,"val"]==max(SANNfit_thisLocus$retvals[,"val"])),]
    
    SANNfit_thisLocus_backup <- SANNfit_thisLocus
    #SANNfit_thisLocus <- SANNfit_thisLocus_backup
    # Burnin remove
    if (burnin>0) {
      SANNfit_thisLocus$retvals <- SANNfit_thisLocus$retvals[(round((burnin/returnFreq),0)):nrow(SANNfit_thisLocus$retvals),]
    }
    # around 50% accept rate good
    # thisRow <- 6
    acceptRate <- round(sum(SANNfit_thisLocus$retvals[,"accept"])/nrow(SANNfit_thisLocus$retvals),3)
    modelSummary[thisRow,"acceptRate"] <- acceptRate
    #modelSummary[thisRow,"locusName"] <- locusName
    # Best estimate
    modelSummary[thisRow,"lnLmax"] <- SANNfit_thisLocus$minimum[1]
    modelSummary[thisRow,"centre"] <- SANNfit_thisLocus$estimate[1]
    modelSummary[thisRow,"width"] <- SANNfit_thisLocus$estimate[2]
    modelSummary[thisRow,"pL"] <- SANNfit_thisLocus$estimate[3]
    modelSummary[thisRow,"pR"] <- SANNfit_thisLocus$estimate[4]
    modelSummary[thisRow,"p0"] <- SANNfit_thisLocus$estimate[3]
    modelSummary[thisRow,"p1"] <- SANNfit_thisLocus$estimate[4]
    modelSummary[thisRow,"Fst"] <- 1/(SANNfit_thisLocus$estimate[5]+1)
    
    # Best estimate from remaining values returned
    # Note, if thinning >1 and burnin >0, you may loose the max lnL in the returned values
    SANNfit_thisLocus_retvals <- as.data.frame(SANNfit_thisLocus$retvals)
    SANNfit_thisLocus_retvals_reject <- SANNfit_thisLocus_retvals[SANNfit_thisLocus_retvals[,"accept"]==0,]
    SANNfit_thisLocus_retvals_accept <- SANNfit_thisLocus_retvals[SANNfit_thisLocus_retvals[,"accept"]==1,]
    modelSummary[thisRow,"centre_accept"] <- SANNfit_thisLocus_retvals_accept[SANNfit_thisLocus_retvals_accept[,"val"]==max(SANNfit_thisLocus_retvals_accept[,"val"]),"p1"]
    modelSummary[thisRow,"width_accept"] <- SANNfit_thisLocus_retvals_accept[SANNfit_thisLocus_retvals_accept[,"val"]==max(SANNfit_thisLocus_retvals_accept[,"val"]),"p2"]
    modelSummary[thisRow,"pL_accept"] <- SANNfit_thisLocus_retvals_accept[SANNfit_thisLocus_retvals_accept[,"val"]==max(SANNfit_thisLocus_retvals_accept[,"val"]),"p3"]
    modelSummary[thisRow,"pR_accept"] <- SANNfit_thisLocus_retvals_accept[SANNfit_thisLocus_retvals_accept[,"val"]==max(SANNfit_thisLocus_retvals_accept[,"val"]),"p4"]
    modelSummary[thisRow,"Fst_accept"] <- (1/(SANNfit_thisLocus_retvals_accept[SANNfit_thisLocus_retvals_accept[,"val"]==max(SANNfit_thisLocus_retvals_accept[,"val"]),"p5"]+1))
    #alleleFreq_Demes <- exportAlleleFreqDemes(DemeData)
    #write.csv(as.data.frame(alleleFreq_Demes),"alleleFreq_Demes.csv",row.names=F)
    
    # trace plots
    tracePlot(SANNfit_thisLocus,locusName,acceptRate)
    
    # Likelihood surfaces and support limits
    modelSum_thisRow <- LL_credibleRegions(SANNfit_thisLocus,locusName=thisLocus,SANNfit_thisLocus_retvals_accept,modelSummary[thisRow,])
    
    # Output trace plots and likelihood surfaces to file
    if (tracePlots==T) {
      # Trace plots
      if (fileFormat=="eps") {
        lastPart_trace <- paste0("rep_",rep,"_gradient_",round(thisGradient,3),"_intercept_",thisIntercept)
        fileTraceName <- paste0(locusName,"_","TracePlots_",thisTrial,"_",lastPart_trace,".eps")
        cairo_ps(fileTraceName,width = 24, height = 9.8, pointsize = 14,antialias="default")
        tracePlot(SANNfit_thisLocus,locusName,acceptRate)
        dev.off()
      }
      if (fileFormat=="png") {
        lastPart_trace <- paste0("rep_",rep,"_gradient_",round(thisGradient,3),"_intercept_",thisIntercept)
        fileTraceName <- paste0(locusName,"_","TracePlots_",thisTrial,"_",lastPart_trace,".png")
        png(fileTraceName,width = 800, height = 600, pointsize = 14)
        tracePlot(SANNfit_thisLocus,locusName,acceptRate)
        dev.off()
      }
      # Likelihood surfaces and support limits
      if (fileFormat=="eps") {
        lastPart <- paste0("rep_",rep,"_gradient_",round(thisGradient,3),"_intercept_",thisIntercept)
        fileLLfitsName <- paste0(locusName,"_","LLfits_",thisTrial,"_",lastPart,".eps")
        cairo_ps(fileLLfitsName,width = 24, height = 9.8, pointsize = 14,antialias="default")
        LL_credible <- LL_credibleRegions(SANNfit_thisLocus,locusName,SANNfit_thisLocus_retvals_accept,modelSummary[thisRow,])
        dev.off()
      }
      if (fileFormat=="png") {
        lastPart <- paste0("rep_",rep,"_gradient_",round(thisGradient,3),"_intercept_",thisIntercept)
        fileLLfitsName <- paste0(locusName,"_","LLfits_",thisTrial,"_",lastPart,".png")
        png(fileLLfitsName,width = 800, height = 600, pointsize = 14,antialias="default")
        LL_credibleRegions(SANNfit_thisLocus,locusName,SANNfit_thisLocus_retvals_accept,modelSummary[thisRow,])
        dev.off()
      }
      # update cline fits to file
      # fileName <- paste(paste("summaryModelsAllloci",thisTrial,sep="_"),".csv",sep="")
      # write.csv(as.data.frame(summaryModelsAllloci),fileName,row.names=F)
    }
    if (plotClines==T) {
      # plot cline
      clinePlot(clineModel="symm_5p",modelSummaryInput=modelSum_thisRow,curveFit=T,cexMain=1.1,cexAxisPlot= 1.85,
                axisLineWdth=2.5,cexPoints=1.8,lwdCurve=4,refLocus=F,refPoints=F,HeaderDetails=T,FigLabel="",multipleCurves=F,
                multiplePoints=F,HaplotypeData=NA)      
      # as pdf
      lastPart_cline <- paste0("rep_",rep,"_gradient_",round(thisGradient,3),"_intercept_",thisIntercept)
      fileClinePlotName <- paste0(locusName,"_","ClinePlot_",thisTrial,"_",lastPart_cline,".pdf")
      pdf(fileClinePlotName,width = 15, height = 15,pointsize=12)
      clinePlot(clineModel="symm_5p",modelSummaryInput=modelSum_thisRow,curveFit=T,cexMain=1.1,cexAxisPlot= 2.1,
                axisLineWdth=3.8,cexPoints=3,lwdCurve=4.3,refLocus=F,refPoints=F,HeaderDetails=T,FigLabel="",HaplotypeData=NA)
      dev.off()
      # as png
      lastPart_cline <- paste0("rep_",rep,"_gradient_",round(thisGradient,3),"_intercept_",thisIntercept)
      fileClinePlotName <- paste0(locusName,"_","ClinePlot_",thisTrial,"_",lastPart_cline,".png")
      png(fileClinePlotName,width = 1000, height = 1000, pointsize = 14,antialias="default")
      clinePlot(clineModel="symm_5p",modelSummaryInput=modelSum_thisRow,curveFit=T,cexMain=1.1,cexAxisPlot= 2.1,
                axisLineWdth=3.8,cexPoints=3,lwdCurve=4.3,refLocus=F,refPoints=F,HeaderDetails=T,FigLabel="",HaplotypeData=NA)
      dev.off()
      obsVals_thisLocus <- NULL
      curveFit <- NULL
    }
    
   modelSum_thisRow[,"N_total"] <- sum(Ngenes)/2
   modelSum_thisRow[,"maxLL"] <- modelSum_thisRow[,"lnLmax"]
   modelSum_thisRow[,"p0"] <- modelSum_thisRow[,"pL"]
   modelSum_thisRow[,"p1"] <- modelSum_thisRow[,"pR"]
   modelSummary[thisRow,]<- modelSum_thisRow
   write.csv(as.data.frame(modelSummary),paste0(thisLocus,"_","modelSummary_",thisTrial,".csv"),row.names=F)
   # write.csv(as.data.frame(modelSummary),"SULF_BradleyetALfreqs.csv",row.names=F)
   
  }
  # clineFit_scan_s1187_290152 <- modelSummary
  return(modelSummary)
}
