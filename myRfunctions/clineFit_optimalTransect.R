# Cline fit optimal transect
# For repeating cline fitting over a range of transect orientations

# David Field 14/2/2024
clineFit_optimalTransect <- function(inData,modelSummary,thisTrial,thisPop,dpThresh,demeSize,minSample,
                              removePhenoMiss,clineModel="symm_5p",transectGradient, transectIntercept, adjustXY,pointsCol,gridAdd,
                              plotAllele1D,plotDem,addPie,SNPlociListUpdated,plotAxisLab=F,minDemeSize,demeMaxSize,
                              returnFreq,burnin,iterations,startParams,fileFormat="png",plotMarkDownAllDemes=F,
                              plotMarkDownZoomed=F,plotEastNorthing=T) {
   inData <- GenoEcol_Planoles 
   #thisLocus <- modelSummary[,"locus"]
   #rep <- modelSummary[,"rep"]
   thisTrial <- "rep1"; clineModel="symm_5p"; rep=1
   modelSummary <- lnLLsummary
  
   thisPop <- "Planoles"; dpThresh <- 0.6
   demeSize=150;minSample=10;removePhenoMiss=F;adjustXY=F
   pointsCol=F; gridAdd=F; plotAllele1D=T;plotDem=c(3,3);addPie=T;plotAxisLab=F;
   minDemeSize=10; demeMaxSize=30; plotPerp=T
   returnFreq=1;burnin=1000;iterations=8000;startParams=c(10000,3000,0.1,0.9,15);fileFormat="png"
  
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
  # 
  # lociNames <- c("s1187_290152","ros_assembly_543443")
  # lociNames <- c("ros_assembly_543443")
  # locusName <- c("s1187_290152")
  # demeSizeSelect <- 100
  # transectGradient = thisGradient;transectIntercept = thisIntercept 
   
   # new 2020
  # DemeData <- DemeSetup2D(HZdata=inData,lociNames=locusSet,demeSize=demeSize, minSample=10, 
  #                         removePhenoMiss=F, transectGradient = thisGradient, 
  #                         transectIntercept = thisIntercept, adjustXY=F, pointsCol=F,gridAdd=T, 
  #                         plotAllele1D=T, plotDem=c(3,3), addPie=T, SNPlociListUpdated, verbose = F, 
  #                         plotAxisLab=T, plotPerp=T, poolSeqPlot=T, plotType="pdf", plotMap=F,plotMarkDownAllDemes=T,
  #                         plotMarkDownZoomed=T,plotEastNorthing=T) 
  
  # New 2024
  
  for (thisRow in 1:nrow(modelSummary)) {
    # thisRow <- 1
    thisLocus <- modelSummary[,"locus"]
    rep <- modelSummary[,"rep"]
    thisGradient <- as.vector(modelSummary[thisRow,"gradient"])
    thisIntercept <- as.vector(modelSummary[thisRow,"intercept"])
    if (thisRow==1) {
      DemeData <- DemeSetup2D(HZdata=inData,lociNames=thisLocus,demeSize=demeSize,minSample=minDemeSize,
                              removePhenoMiss=F,transectGradient = thisGradient, transectIntercept = thisIntercept,  
                              adjustXY=F, pointsCol = F, gridAdd=F,plotAllele1D = T, plotDem=c(3,3), 
                              addPie=T,SNPlociListUpdated,verbose=F,plotAxisLab=F,plotPerp=T,plotType="pdf",
                              poolSeqPlot=T,plotMap=F,noPlots=F,plotMarkDownAllDemes=F,plotEastNorthing=T,
                              plotMarkDownZoomed=F,plotTerrain=F,plotKASPrep=T,FigName="Chr1",plotCentreLine=T,
                              poolSeqClines=allChr_Clines_Chr1clines,poolSeqStats=poolLoc,colourTable=colourChoice)
      thisDataLoc <- DemeData[[which(names(DemeData)==thisLocus)]]
    }
    if (thisRow>1) {
      if (thisLocus != modelSummary[thisRow-1,"locus"] | thisGradient!=as.vector(modelSummary[thisRow-1,"gradient"])) {
        rm(DemeData)
        DemeData <- DemeSetup2D(HZdata=inData,lociNames=thisLocus,demeSize=demeSize,minSample=minDemeSize,
                                removePhenoMiss=F,transectGradient = thisGradient, transectIntercept = thisIntercept,  
                                adjustXY=F, pointsCol = F, gridAdd=F,plotAllele1D = T, plotDem=c(3,3), 
                                addPie=T,SNPlociListUpdated,verbose=F,plotAxisLab=F,plotPerp=T,plotType="pdf",
                                poolSeqPlot=T,plotMap=F,noPlots=F,plotMarkDownAllDemes=F,plotEastNorthing=T,
                                plotMarkDownZoomed=F,plotTerrain=F,plotKASPrep=T,FigName="Chr1",plotCentreLine=T,
                                poolSeqClines=allChr_Clines_Chr1clines,poolSeqStats=poolLoc,colourTable=colourChoice)
        thisDataLoc <- DemeData[[which(names(DemeData)==thisLocus)]]
      }
    }
    
    # NEW approach 2020 and updated again for 2024
    #names(thisDataLoc)
    XCordinates <- thisDataLoc[["DemeMatrix2D_DistanceRefLine"]][,2]
    thisLoc_DemeMatrix2D_NoZero <- thisDataLoc[["DemeMatrix2D_combine_noZero"]]
    #colnames(thisLoc_DemeMatrix2D_NoZero)
    #nrow(thisLoc_DemeMatrix2D_NoZero)==nrow(XCordinates)
    # colnames(thisLoc_DemeMatrix2D_NoZero)
    # mininum sample size
    sampleSize <- (thisLoc_DemeMatrix2D_NoZero[,"pCount"]+thisLoc_DemeMatrix2D_NoZero[,"qCount"])/2
    theseDemesGood <- which(sampleSize>=minDemeSize)
    TotalSamples <- sum(sampleSize[sampleSize>=minDemeSize])
    
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
    striatum_p <- mean(thisLoc_DemeMatrix2D_NoZero[1:6,"p"])
    pseudom_p <- mean(thisLoc_DemeMatrix2D_NoZero[(nrow(thisLoc_DemeMatrix2D_NoZero)-6):(nrow(thisLoc_DemeMatrix2D_NoZero)),"p"])
    deltaP <- abs(pseudom_p-striatum_p)
    if (deltaP<0.6) {
      next
    }
    
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
    pCount <- pCount_raw
    qCount <- qCount_raw
    Ngenes <- Ngenes_raw
    p_raw <- pCount_raw/Ngenes_raw
    
    # reduce size of exceptionally large N for some demes (subsample) 
    # Note - otherwise computational problem with pochhammer function
    # demeMaxSize <- 40
    if (any(Ngenes > demeMaxSize)) {
      largeNDemes <- which(Ngenes > demeMaxSize)
      for (thisDeme in 1:length(largeNDemes)) {
        # thisDeme <- 2
        thisN <- Ngenes[largeNDemes[thisDeme]]
        thisP <- pCount[largeNDemes[thisDeme]]
        thisQ <- qCount[largeNDemes[thisDeme]]
        subsample <- sample(c(1,0),demeMaxSize,prob=c((thisP/thisN),(thisQ/thisN)),replace=T)
        pCount[largeNDemes[thisDeme]] <- sum(subsample==1)
        qCount[largeNDemes[thisDeme]] <- sum(subsample==0)
        Ngenes[largeNDemes[thisDeme]] <- sum(subsample==1) + sum(subsample==0)
      }
    }
    
    obsVals <- cbind(pCount,qCount,Ngenes)
    
    # # new adjusted allele frequency
    p_adj <- pCount/Ngenes
    plot(p_raw,p_adj)
    # Send deme data to summary table
    # summaryModelsAllloci[thisDrop,"N_total"] <- TotalSamples
    # summaryModelsAllloci[thisDrop,"N_demes_mean"] <- mean(Ngenes)
    # summaryModelsAllloci[thisDrop,"N_demes_sd"] <- sd(Ngenes)
    # summaryModelsAllloci[thisDrop,"Xvals"] <- paste(XCordinates,collapse=";")
    # summaryModelsAllloci[thisDrop,"Ngenes"] <- paste(Ngenes,collapse=";")
    # summaryModelsAllloci[thisDrop,"p"] <- paste(p_adj,collapse=";")
    # 
    
    modelSummary[thisRow,"locus"] <- thisLocus
    modelSummary[thisRow,"demeSize"] <- demeSize
    modelSummary[thisRow,"numDemes"] <- length(pCount)
    modelSummary[thisRow,"rep"] <- 1
    modelSummary[thisRow,"iters"] <- iterations
    modelSummary[thisRow,"burnin"] <- burnin
    modelSummary[thisRow,"N_total"] <- TotalSamples
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
    
    # Burnin remove
    if (burnin>0) {
      SANNfit_thisLocus$retvals <- SANNfit_thisLocus$retvals[(round((burnin/returnFreq),0)):nrow(SANNfit_thisLocus$retvals),]
    }
    # around 50% accept rate good
    acceptRate <- round(sum(SANNfit_thisLocus$retvals[,"accept"])/nrow(SANNfit_thisLocus$retvals),3)
    summaryModelsAllloci[thisDrop,"locusName"] <- locusName
    # Best estimate
    modelSummary[thisRow,"lnLmax"] <- SANNfit_thisLocus$minimum[1]
    modelSummary[thisRow,"centre"] <- SANNfit_thisLocus$estimate[1]
    modelSummary[thisRow,"width"] <- SANNfit_thisLocus$estimate[2]
    modelSummary[thisRow,"pL"] <- SANNfit_thisLocus$estimate[3]
    modelSummary[thisRow,"pR"] <- SANNfit_thisLocus$estimate[4]
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
    
    
    
  }
  
  
  
  
  
  # str(HZ_50mDemes)
  
  
  # names(DemeData_200m)
  # save allele frequencies to file 
  # for each locus (need this for RosEl manuscript re-plotting)
  alleleFreq_Demes <- exportAlleleFreqDemes(DemeData)
  write.csv(as.data.frame(alleleFreq_Demes),"alleleFreq_Demes.csv",row.names=F)
  
  fullLociNames <- colnames(inData)[(which(colnames(inData)=="s1187_290152")):(which(colnames(inData)=="s735_701957"))]
  summaryModelsAllloci <- matrix(0,length(fullLociNames),40)
  colnames(summaryModelsAllloci) <- c("locusName","LG","cM","pop","demeSize","minN","maxN","iters","model","rep",
                                      "N_total","N_demes_mean","N_demes_sd","deltaP",
                                      "lnLmax","centre","width","pL","pR","Fst",
                                      "centre_accept","width_accept","pL_accept","pR_accept","Fst_accept",
                                      "centre_L95","centre_U95","width_L95","width_U95","pL_L95",
                                      "pL_U95","pR_L95","pR_U95","Fst_L95","Fst_U95","Xvals","Ngenes","p","pointCol","curveCol")
  summaryModelsAllloci[,"locusName"] <- fullLociNames
  summaryModelsAllloci[,"demeSize"] <- demeSize
  summaryModelsAllloci[,"minN"] <- minSample
  summaryModelsAllloci[,"maxN"] <- demeMaxSize
  summaryModelsAllloci[,"iters"] <- iterations
  summaryModelsAllloci[,"model"] <- clineModel
  summaryModelsAllloci[,"rep"] <- rep
  summaryModelsAllloci[,"pop"] <- thisPop
  summaryModelsAllloci[,"pointCol"] <- "gray"
  summaryModelsAllloci[,"curveCol"] <- "black"
  for (thisLocus in 1:length(locusSet)) {
    # thisLocus <- 1
    locusName <- locusSet[thisLocus]
    thisDrop <- which(summaryModelsAllloci[,"locusName"]==locusName)
    summaryModelsAllloci[thisDrop,"LG"] <- as.vector(SNPlociListUpdated[which(SNPlociListUpdated[,"LocusName"]==locusName),"LG"])
    summaryModelsAllloci[thisDrop,"cM"] <- as.vector(SNPlociListUpdated[which(SNPlociListUpdated[,"LocusName"]==locusName),"cM"])
    # extract this locus deme data
    if (locusName%in%names(DemeData)) {
      thisDataLoc <- DemeData[[which(names(DemeData)==locusName)]]
      # position of demes along transect
      # XCordinates <- thisDataLoc[["DemeMatrix2D_pqCountNoZero"]][,"distAlongTransect"]
      # all data for each deme (non zero samples)
      thisLoc_DemeMatrix2D_pqCountNoZero <- thisDataLoc[["DemeMatrix2D_combine_noZero"]]
      
      # mininum sample size
      sampleSize <- as.vector(apply(thisLoc_DemeMatrix2D_pqCountNoZero[,6:7],1,sum)/2)
      #cat("\n sampleSize: ",sampleSize)
      theseDemesGood <- which(sampleSize>=minDemeSize)
      #cat("\n theseDemesGood: ",theseDemesGood)
      TotalSamples <- sum(sampleSize[sampleSize>=minDemeSize])
      
      # reduce to demes with sufficient sample size (minDemeSize)
      thisLoc_Demes_good <- thisLoc_DemeMatrix2D_pqCountNoZero[theseDemesGood,]
      # X cordinates 
      XCordinates <- thisLoc_Demes_good[,"distAlongTransect"]
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
      pCount <- pCount_raw
      qCount <- qCount_raw
      Ngenes <- Ngenes_raw
      # reduce size of exceptionally large N for some demes (subsample) 
      # Note - otherwise computational problem with pochhammer function
      if (any(Ngenes > demeMaxSize)) {
        largeNDemes <- which(Ngenes > demeMaxSize)
        for (thisDeme in 1:length(largeNDemes)) {
          # thisDeme <- 2
          thisN <- Ngenes[largeNDemes[thisDeme]]
          thisP <- pCount[largeNDemes[thisDeme]]
          thisQ <- qCount[largeNDemes[thisDeme]]
          subsample <- sample(c(1,0),demeMaxSize,prob=c((thisP/thisN),(thisQ/thisN)),replace=T)
          pCount[largeNDemes[thisDeme]] <- sum(subsample==1)
          qCount[largeNDemes[thisDeme]] <- sum(subsample==0)
          Ngenes[largeNDemes[thisDeme]] <- sum(subsample==1) + sum(subsample==0)
        }
      }
      # Delta p
      striatum_p <- mean(thisLoc_DemeMatrix2D_pqCountNoZero[1:6,"p"])
      pseudom_p <- mean(thisLoc_DemeMatrix2D_pqCountNoZero[(nrow(thisLoc_DemeMatrix2D_pqCountNoZero)-6):(nrow(thisLoc_DemeMatrix2D_pqCountNoZero)),"p"])
      deltaP <- abs(pseudom_p-striatum_p)
      summaryModelsAllloci[thisDrop,"deltaP"] <- deltaP
      obsVals_thisLocus <- cbind(pCount,qCount,Ngenes)
      # new adjusted allele frequency
      p_adj <- pCount/Ngenes
      
      # Send deme data to summary table
      summaryModelsAllloci[thisDrop,"N_total"] <- TotalSamples
      summaryModelsAllloci[thisDrop,"N_demes_mean"] <- mean(Ngenes)
      summaryModelsAllloci[thisDrop,"N_demes_sd"] <- sd(Ngenes)
      summaryModelsAllloci[thisDrop,"Xvals"] <- paste(XCordinates,collapse=";")
      summaryModelsAllloci[thisDrop,"Ngenes"] <- paste(Ngenes,collapse=";")
      summaryModelsAllloci[thisDrop,"p"] <- paste(p_adj,collapse=";")
      # Main SANN routine
      cat(c("\n Cline fitting SANN"))
      cat(c("\n ",locusName))
      cat(c("."))
      cat(c("deltaP",round(deltaP,3)))
      cat(c("; "))
      
      # skip loci without strong allele frequency differences
      if (deltaP < dpThresh) {
        # plot cline
        clinePlot(clineModel="no fit",locusName,clineSummary=summaryModelsAllloci,curveFit=F,cexMain=1.3,cexAxisPlot= 1.85,
                  axisLineWdth=2.5,cexPoints=1.8,lwdCurve=4)      
        # as pdf
        fileClineFitName <- paste(paste(paste(locusName,"ClineFit",sep="_"),thisTrial,sep="_"),".pdf",sep="")
        pdf(fileClineFitName,width = 15, height = 15,pointsize=12)
        clinePlot(clineModel="not fit",locusName,summaryModelsAllloci,curveFit=F,cexMain=2.5,cexAxisPlot= 2.7,
                  axisLineWdth=3.8,cexPoints=3,lwdCurve=4.3)
        dev.off()
        # # as tiff
        # fileClineFitName <- paste(paste(paste(locusName,"ClineFit",sep="_"),thisTrial,sep="_"),".tiff",sep="")
        # tiff(paste(folder_output,fileClineFitName,sep=""),width = 1500, height = 1500,pointsize=7,res=400)
        # clinePlot(clineModel="",locusName,summaryModelsAllloci,curveFit=T,cexMain=1,cexAxisPlot= 1.1,
        #           axisLineWdth=1.5,cexPoints=1.2,lwdCurve=2)      
        # dev.off()
        # as png
        fileClineFitName <- paste(paste(paste(locusName,"ClineFit",sep="_"),thisTrial,sep="_"),".png",sep="")
        png(fileClineFitName,width = 1000, height = 1000, pointsize = 14,antialias="default")
        clinePlot(clineModel="no fit",locusName,summaryModelsAllloci,curveFit=F,cexMain=1.4,cexAxisPlot= 2,
                  axisLineWdth=2,cexPoints=1.8,lwdCurve=4)
        dev.off()
        obsVals_thisLocus <- NULL
        curveFit <- NULL
        # update cline fits to file
        fileName <- paste(paste("summaryModelsAllloci",thisTrial,sep="_"),".csv",sep="")
        write.csv(as.data.frame(summaryModelsAllloci),fileName,row.names=F)
        obsVals_thisLocus <- NULL
        curveFit <- NULL
        cat(c(".skip locus for fitting"))
        next
      }
      
      # NEW approach 2020
      XCordinates <- thisDataLoc[["DemeMatrix2D_DistanceRefLine"]][,2]
      thisLoc_DemeMatrix2D_NoZero <- thisDataLoc[["DemeMatrix2D_combine_noZero"]]
      #nrow(thisLoc_DemeMatrix2D_NoZero)==nrow(XCordinates)
      # colnames(thisLoc_DemeMatrix2D_NoZero)
      # mininum sample size
      sampleSize <- as.vector(apply(thisLoc_DemeMatrix2D_NoZero[,6:7],1,sum)/2)
      theseDemesGood <- which(sampleSize>=minDemeSize)
      TotalSamples <- sum(sampleSize[sampleSize>=minDemeSize])
      
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
      
      # skip loci without strong allele frequency differences
      striatum_p <- mean(thisLoc_DemeMatrix2D_NoZero[1:6,"p"])
      pseudom_p <- mean(thisLoc_DemeMatrix2D_NoZero[(nrow(thisLoc_DemeMatrix2D_NoZero)-6):(nrow(thisLoc_DemeMatrix2D_NoZero)),"p"])
      deltaP <- abs(pseudom_p-striatum_p)
      if (deltaP<0.6) {
        next
      }
      
      # reduce size of exceptionally large N for some demes (subsample) 
      # Note - otherwise computational problem with pochhammer function
      if (any(Ngenes_raw > demeMaxSize)) {
        largeNDemes <- which(Ngenes_raw > demeMaxSize)
        for (thisDeme in 1:length(largeNDemes)) {
          # thisDeme <- 2
          thisN <- Ngenes_raw[largeNDemes[thisDeme]]
          thisP <- pCount_raw[largeNDemes[thisDeme]]
          thisQ <- qCount_raw[largeNDemes[thisDeme]]
          subsample <- sample(c(1,0),demeMaxSize,prob=c((thisP/thisN),(thisQ/thisN)),replace=T)
          pCount[largeNDemes[thisDeme]] <- sum(subsample==1)
          qCount[largeNDemes[thisDeme]] <- sum(subsample==0)
          Ngenes[largeNDemes[thisDeme]] <- sum(subsample==1) + sum(subsample==0)
        }
      }
      # Main SANN routine
      # SANNfit_thisLocus <- SANN_clinefit(fn = logLik_5par_SANN, positions=XCordinates, obsVals,initialParam = startParams,deltap = NULL, scale=0.4,retvals=TRUE,nmax = 5000 ,verbose=TRUE,acceptscale = 1.01, rejectscale = 0.99, retfreq=returnFreq)
      
      # old routine from slowClines_RosEl 2018 Aug (which works)
      cat(c("..Cline fitting.."))
      SANNfit_thisLocus <- SANN_clinefit(fn = logLik_5par_SANN,obsVals, positions=XCordinates, 
                                         initialParam = startParams,scale=0.3,retvals=TRUE,nmax = iterations,
                                         verbose=TRUE,acceptscale = 1.05, rejectscale = (1/1.05), 
                                         retfreq=returnFreq)
      # SANNfit_thisLocus$retvals[which(SANNfit_thisLocus$retvals[,"val"]==max(SANNfit_thisLocus$retvals[,"val"])),]
      
      SANNfit_thisLocus_backup <- SANNfit_thisLocus
      # Burnin remove
      if (burnin>0) {
        SANNfit_thisLocus$retvals <- SANNfit_thisLocus$retvals[(round((burnin/returnFreq),0)):nrow(SANNfit_thisLocus$retvals),]
      }
      acceptRate <- round(sum(SANNfit_thisLocus$retvals[,"accept"])/nrow(SANNfit_thisLocus$retvals),3)
      summaryModelsAllloci[thisDrop,"locusName"] <- locusName
      # Best estimate
      summaryModelsAllloci[thisDrop,"lnLmax"] <- SANNfit_thisLocus$minimum[1]
      summaryModelsAllloci[thisDrop,"centre"] <- SANNfit_thisLocus$estimate[1]
      summaryModelsAllloci[thisDrop,"width"] <- SANNfit_thisLocus$estimate[2]
      summaryModelsAllloci[thisDrop,"pL"] <- SANNfit_thisLocus$estimate[3]
      summaryModelsAllloci[thisDrop,"pR"] <- SANNfit_thisLocus$estimate[4]
      summaryModelsAllloci[thisDrop,"Fst"] <- 1/(SANNfit_thisLocus$estimate[5]+1)
      # Best estimate from remaining values returned
      # Note, if thinning >1 and burnin >0, you may loose the max lnL in the returned values
      SANNfit_thisLocus_retvals <- as.data.frame(SANNfit_thisLocus$retvals)
      SANNfit_thisLocus_retvals_reject <- SANNfit_thisLocus_retvals[SANNfit_thisLocus_retvals[,"accept"]==0,]
      SANNfit_thisLocus_retvals_accept <- SANNfit_thisLocus_retvals[SANNfit_thisLocus_retvals[,"accept"]==1,]
      summaryModelsAllloci[thisDrop,"centre_accept"] <- SANNfit_thisLocus_retvals_accept[SANNfit_thisLocus_retvals_accept[,"val"]==max(SANNfit_thisLocus_retvals_accept[,"val"]),"p1"]
      summaryModelsAllloci[thisDrop,"width_accept"] <- SANNfit_thisLocus_retvals_accept[SANNfit_thisLocus_retvals_accept[,"val"]==max(SANNfit_thisLocus_retvals_accept[,"val"]),"p2"]
      summaryModelsAllloci[thisDrop,"pL_accept"] <- SANNfit_thisLocus_retvals_accept[SANNfit_thisLocus_retvals_accept[,"val"]==max(SANNfit_thisLocus_retvals_accept[,"val"]),"p3"]
      summaryModelsAllloci[thisDrop,"pR_accept"] <- SANNfit_thisLocus_retvals_accept[SANNfit_thisLocus_retvals_accept[,"val"]==max(SANNfit_thisLocus_retvals_accept[,"val"]),"p4"]
      summaryModelsAllloci[thisDrop,"Fst_accept"] <- (1/(SANNfit_thisLocus_retvals_accept[SANNfit_thisLocus_retvals_accept[,"val"]==max(SANNfit_thisLocus_retvals_accept[,"val"]),"p5"]+1))
      # Trace plots
      tracePlot(SANNfit_thisLocus,locusName,acceptRate)
      if (fileFormat=="eps") {
        fileTraceName <- paste(paste(paste(locusName,"Trace",sep="_"),thisTrial,sep="_"),".eps",sep="")
        cairo_ps(fileTraceName,width = 24, height = 9.8, pointsize = 14,antialias="default")
        tracePlot(SANNfit_thisLocus,locusName,acceptRate)
        dev.off()
      }
      if (fileFormat=="png") {
        fileTraceName <- paste(paste(paste(locusName,"Trace",sep="_"),thisTrial,sep="_"),".png",sep="")
        png(fileTraceName,width = 800, height = 600, pointsize = 14)
        tracePlot(SANNfit_thisLocus,locusName,acceptRate)
        dev.off()
      }
      # Likelihood surfaces and support limits
      summaryModelsAllloci <- LL_credibleRegions(SANNfit_thisLocus,locusName,SANNfit_thisLocus_retvals_accept,summaryModelsAllloci)
      if (fileFormat=="eps") {
        fileLLfitsName <- paste(paste(paste(locusName,"LLfits",sep="_"),thisTrial,sep="_"),".eps",sep="")
        cairo_ps(fileLLfitsName,width = 24, height = 9.8, pointsize = 14,antialias="default")
        LL_credibleRegions(SANNfit_thisLocus,locusName,SANNfit_thisLocus_retvals_accept,summaryModelsAllloci)
        dev.off()
      }
      if (fileFormat=="png") {
        fileLLfitsName <- paste(paste(paste(locusName,"LLfits",sep="_"),thisTrial,sep="_"),".png",sep="")
        png(fileLLfitsName,width = 800, height = 600, pointsize = 14,antialias="default")
        LL_credibleRegions(SANNfit_thisLocus,locusName,SANNfit_thisLocus_retvals_accept,summaryModelsAllloci)
        dev.off()
      }
      # update cline fits to file
      fileName <- paste(paste("summaryModelsAllloci",thisTrial,sep="_"),".csv",sep="")
      write.csv(as.data.frame(summaryModelsAllloci),fileName,row.names=F)
      # plot cline
      clinePlot(clineModel="symm_5p",locusName,summaryModelsAllloci,curveFit=T,cexMain=1.3,cexAxisPlot= 1.85,
                axisLineWdth=2.5,cexPoints=1.8,lwdCurve=4)      
      # as pdf
      fileClineFitName <- paste(paste(paste(locusName,"ClineFit",sep="_"),thisTrial,sep="_"),".pdf",sep="")
      pdf(fileClineFitName,width = 15, height = 15,pointsize=12)
      clinePlot(clineModel="symm_5p",locusName,summaryModelsAllloci,curveFit=T,cexMain=2.5,cexAxisPlot= 2.7,
                axisLineWdth=3.8,cexPoints=3,lwdCurve=4.3)
      dev.off()
      # # as tiff
      # fileClineFitName <- paste(paste(paste(locusName,"ClineFit",sep="_"),thisTrial,sep="_"),".tiff",sep="")
      # tiff(paste(folder_output,fileClineFitName,sep=""),width = 1500, height = 1500,pointsize=7,res=400)
      # clinePlot(clineModel="symm_5p",locusName,summaryModelsAllloci,curveFit=T,cexMain=1,cexAxisPlot= 1.1,
      #           axisLineWdth=1.5,cexPoints=1.2,lwdCurve=2)      
      # dev.off()
      # as png
      fileClineFitName <- paste(paste(paste(locusName,"ClineFit",sep="_"),thisTrial,sep="_"),".png",sep="")
      png(fileClineFitName,width = 1000, height = 1000, pointsize = 14,antialias="default")
      clinePlot(clineModel="symm_5p",locusName,summaryModelsAllloci,curveFit=T,cexMain=1.4,cexAxisPlot= 2,
                axisLineWdth=2,cexPoints=1.8,lwdCurve=4)
      dev.off()
      obsVals_thisLocus <- NULL
      curveFit <- NULL
      
    }
    if (locusName%in%names(HZ_200mDemes)==F) {
      summaryModelsAllloci[thisDrop,"N_total"] <- "noData"
    }
  }
  summaryModelsAllloci[summaryModelsAllloci[,"N_total"]=="noData","N_total"] <- 0
  return(summaryModelsAllloci)
}
