# DemeSetup2D_genos_PoolSeq
# Author: David. L. Field
# 1/4/2024

# This is a function for separating individuals into demes in two-dimensions and calculating various pop gen statistics.
# This function designed for genotypes and allele frequencies across a set of SNP loci.
# Produces a list for each deme, containing a datatable with estimates various values such as allele frequencies, 
# genotype counts etc.. in demes in two-dimensions. 
# Note- this version works with PoolSeq data preparations for the WGS clines paper
# Note - there is another version which works for slowClines
# lociNames = which loci to include
# demeSize = in metres
# minSample = minimum sample size for calculating pop gen estimates
# removePhenoMiss = T/F, whether to remove individuals with no colour phenotypes
# transectGradient = for clines analyses
# transectIntercept = for cline analyses
# adjustXY = T/F, whether to reset the mininum X and Y spatial values to 0
# plotProgress = T/F, whether to show allele frequencies in plots for each locus
# pointsCol = T/F, colour points for individuals based on colour phenotypes
# gridAdd = T/F, add a grid based on deme size in the background of 2D map plots
# plotAxisLab = 
# plotPerp = show perpendicular intersects between deme and transect (demostrate flattening from 2D to 1D)
# poolSeqPlot = shows location of poolSeq (Tavares et al 2018 PNAS) along the transect
# plotType = 'pdf'
# plotMap = 

DemeSetup2D_PoolSeq <- function(HZdata,lociNames,demeSize,minSample,removePhenoMiss,transectGradient,transectIntercept,
                        adjustXY,pointsCol,gridAdd,plotAllele1D,plotDem,addPie,SNPlociListUpdated,verbose=F,
                        plotAxisLab,plotPerp,poolSeqPlot,plotType,plotMap,noPlots=T,plotMarkDownAllDemes,
                        plotMarkDownZoomed,plotEastNorthing,plotTerrain,plotKASPrep,FigName,plotCentreLine,poolSeqClines,
                        poolSeqStats,plotKASPcurve=T,colourTable) {
  
  # nrow(HZdata)
  # HZdata <- GenoEcol_Planoles
  # FigName=paste0("_ModelRow_",thisRow)
  # lociNames=thisLocus
  # lociNames <- c("s1187_290152")
  # lociNames <- c("s1187_290152","s316_93292", "s316_257789","ros_assembly_543443")
  # plotAxisLab <- T ; demeSize <- 200; minSample <- 10; removePhenoMiss <- FALSE; gridAdd <- F; poolSeqPlot=F; adjustXY<- F
  # plotAllele1D <- T; location <- "Planoles"; plotPerp = T; plotType ="pdf"; verbose=T; plotMap=T; addPie=T
  # plotEastNorthing <- T; plotMarkDownAllDemes <- F; plotMarkDownZoomed <- F
  # in Easting
  # transectGradient <- -0.24; transectIntercept <- 4788100 # Planoles
  # transectGradient <- 3.125; transectIntercept <- 3412000 # Cadi
  # in metres
  # transectGradient <- -0.20; transectIntercept <- 4900 # Planoles # this is the final one
  # transectGradient <- -0.24; transectIntercept <- 6000 # Planoles
  # 
  # transectGradient = -0.2; transectIntercept = 4700
  # transectGradient <- 3.125; transectIntercept <- 3412000 # Cadi
  # adjustXY <- FALSE; pointsCol <- FALSE; pieRadius <- 60; plotDem <- c(2,2); addPie <- T; verbose=T; poolSeqPlot=F
  # noPlots <- F; plotMarkDownAllDemes = F; plotMarkDownZoomed = F; plotEastNorthing=T; plotTerrain <- F; plotMap <- F
  # poolSeqClines <- allChr_Clines_Chr5clines; poolSeqStats <- poolSeqStats; FigName="Chr4"; plotKASPrep=T
  # colourTable <- colourChoice
  transectIntercept <- as.vector(transectIntercept)
  transectGradient <- as.vector(transectGradient)
  cat("\nProcessing demes... Total number of individuals at start = ",nrow(HZdata))
  flush.console()
  markerLabelPieChart <- function(thisLocusName,SNPlociListUpdated) {
    # SampleSize <- 1514
    markerLabelMake <- NULL
    SampleLabel <- NULL
    colnames(SNPlociListUpdated)
    linkageGroup <- SNPlociListUpdated[SNPlociListUpdated[,"LocusName"]==thisLocusName,"LG"]
    cM <- SNPlociListUpdated[SNPlociListUpdated[,"LocusName"]==thisLocusName,"cM"]
    fullLabel <- paste(thisLocusName,
                       paste(paste("LG",linkageGroup,sep=", "),sep=", "))
    return(fullLabel)
  }
  
  if (adjustXY==TRUE) {
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
  # remove plants with missing position
  HZdata <- HZdata[!is.na(HZdata[,"Easting"]),]
  HZdata <- HZdata[HZdata[,"Easting"]!=-9,]
  minX <- as.numeric(min(HZdata[,"Easting"]))
  maxX <- as.numeric(max(HZdata[,"Easting"]))
  minY <- as.numeric(min(HZdata[,"Northing"]))
  maxY <- as.numeric(max(HZdata[,"Northing"]))
  
  # Using midpoints of square grid
  Xboundaries <- sort(seq(minX,maxX,demeSize))
  Yboundaries <- sort(seq(minY,maxY,demeSize))
  Xboundaries_midpoint <- Xboundaries+(demeSize/2)
  Yboundaries_midpoint <- Yboundaries+(demeSize/2)
 
  # Deme data stored in a list for each locus
  DemeData <- list()
  if (verbose==T) {
    cat(c(", Locus.."))
    cat(c(", minX = ",round(minX,3)))
    cat(c(", minY = ",round(minY,3)))
  }
  # calculate deme stats for each locus
  for (thisLocus in 1:length(lociNames)) {
    # thisLocus <- 1
    if (verbose==T) {
      cat(c("..",thisLocus))
      flush.console()
    }
    thisLocusName <- lociNames[thisLocus]
    # thisLocusName <- "ros_assembly_717045"
    fullData_thisLocus <- HZdata[,colnames(HZdata)==thisLocusName]
    # backup <- fullData_thisLocus
    summaryfullData_thisLocus <- table(fullData_thisLocus)
    # summaryfullData_thisLocus_backup <- table(backup)

    #sumGenotypes <- summaryfullData_thisLocus
    if (length(summaryfullData_thisLocus)==1) {
      next
    }
    DemeMatrix2D <- matrix(0,(length(Yboundaries)),(length(Xboundaries)))
    colnames(DemeMatrix2D) <- Xboundaries_midpoint
    rownames(DemeMatrix2D) <- rev(Yboundaries_midpoint)
    DemeMatrix2D[,] <- -9
    # allele frequencies
    DemeMatrix2D_p <- DemeMatrix2D
    DemeMatrix2D_q <- DemeMatrix2D
    DemeMatrix2D_pq <- DemeMatrix2D
    # allele counts
    DemeMatrix2D_Ngenes <- DemeMatrix2D
    DemeMatrix2D_pCount <- DemeMatrix2D
    DemeMatrix2D_qCount <- DemeMatrix2D
    # genotype counts
    DemeMatrix2D_ppCount <- DemeMatrix2D
    DemeMatrix2D_qqCount <- DemeMatrix2D
    DemeMatrix2D_pqCount <- DemeMatrix2D
    # genotype frequencies
    DemeMatrix2D_pp <- DemeMatrix2D
    DemeMatrix2D_qq <- DemeMatrix2D
    DemeMatrix2D_pq <- DemeMatrix2D
    
    # Individuals 
    DemeMatrix2D_inds <- DemeMatrix2D
    # Genotypes
    DemeMatrix2D_genos <- DemeMatrix2D
    
    # Keep list of individuals
    DemeIndividuals <- matrix(0,ncol(DemeMatrix2D_p)*nrow(DemeMatrix2D_p),7)
    colnames(DemeIndividuals) <- c("Xmin","Xmax","Ymin","Ymax","numInds","Inds","Genos")
    demeCounter <- 0
    # 
    if (plotMap==T) {
      par(mfrow=c(1,1))
      # cores and flanks - 421000,426000,4685500,4687800
      #plotHZbasic(HZdata,409000,437000,4684500,4692000,TRUE,FALSE,TitleStart=paste("Planoles ",thisLocusName,sep=""),plotAge=F)
      plotHZpedigree(HZdata,409000,437000,4684500,4692000,TRUE,FALSE,TitleStart=paste("Planoles ",thisLocusName,sep=""))
      #myCols <- viridis(30,1,0,1,option="D") # viridis
      myCols <- rainbow(20) # 
    }
    for (thisX in 1:ncol(DemeMatrix2D_p)) {
      for (thisY in 1:nrow(DemeMatrix2D_p)) {
        demeCounter <- demeCounter + 1
        # thisX <- 1
        # thisY <- 1
        # old approach
        #theseX <- as.numeric(unlist(strsplit(colnames(DemeMatrix2D)[thisX],"_")))
        #theseY <- as.numeric(unlist(strsplit(rownames(DemeMatrix2D)[thisY],"_")))
        # new approach
        theseX <- c(as.numeric(colnames(DemeMatrix2D)[thisX])-(demeSize/2),as.numeric(colnames(DemeMatrix2D)[thisX])+(demeSize/2))
        theseY <- c(as.numeric(rownames(DemeMatrix2D)[thisY])-(demeSize/2),as.numeric(rownames(DemeMatrix2D)[thisY])+(demeSize/2))
        theseRows2D <- as.numeric(HZdata[,"Easting"]) >= theseX[1] &
          as.numeric(HZdata[,"Easting"]) <= theseX[2] &
          as.numeric(HZdata[,"Northing"]) <= theseY[2] &
          as.numeric(HZdata[,"Northing"]) >= theseY[1]
        if (all(theseRows2D==FALSE)) {
          next
        }
        #jLoc_kDeme[2] <- "NGY"
        #(all(jLoc_kDeme!="NGY"))
        #HZdataSetXYdata[,colnames(HZdataSetXYdata)==thisLocusName]
        if (any(theseRows2D==TRUE)) {
          #cat("\nthisX: ",thisX)
          #cat(",thisY: ",thisY)
          
          # data for this Locus and this Deme (set by boundary) of Easting and Northing 
          # i.e. the jth Locus and kth Deme
          jLoc_kDeme <- HZdata[theseRows2D,colnames(HZdata)==thisLocusName]
          inds_kDeme <- HZdata[theseRows2D,"PlantID_final"]
          
          # problem with code here (fixed?)
          if (all(jLoc_kDeme==-10 | jLoc_kDeme==-9)) {
            next
          }
          else {
            TotalCount <- (sum(jLoc_kDeme==0) + sum(jLoc_kDeme==1) + sum(jLoc_kDeme==2) + sum(jLoc_kDeme=="-9"))
            Ncount <- (sum(jLoc_kDeme==0) + sum(jLoc_kDeme==1) + sum(jLoc_kDeme==2))
            Ngenes <- 2*sum(sum(jLoc_kDeme==0) + sum(jLoc_kDeme==1) + sum(jLoc_kDeme==2))
            pCount <- (2*sum(jLoc_kDeme==0)) + sum(jLoc_kDeme==1)
            qCount <- (2*sum(jLoc_kDeme==2)) + sum(jLoc_kDeme==1)
            missCount <- sum(jLoc_kDeme=="-9")
            TotalGenotyped <- TotalCount-missCount
            # Genotypes
            ppCount <- sum(jLoc_kDeme==0)
            qqCount <- sum(jLoc_kDeme==2)
            pqCount <- sum(jLoc_kDeme==1)
            # check sums
            #pTotal+qTotal == totalGenes
            p_freq <- ((2*sum(jLoc_kDeme==0))+sum(jLoc_kDeme==1))/(2*sum(sum(jLoc_kDeme==0)+sum(jLoc_kDeme==1)+sum(jLoc_kDeme==2)))
            q_freq <- ((2*sum(jLoc_kDeme==2))+sum(jLoc_kDeme==1))/(2*sum(sum(jLoc_kDeme==0)+sum(jLoc_kDeme==1)+sum(jLoc_kDeme==2)))
            miss_freq <- sum(jLoc_kDeme=="-9")/(sum(jLoc_kDeme==0)+sum(jLoc_kDeme==1)+sum(jLoc_kDeme==2)+sum(jLoc_kDeme=="-9"))
            Homs1 <- sum(jLoc_kDeme==0)
            Homs2 <- sum(jLoc_kDeme==2)
            Hets <- sum(jLoc_kDeme==1)
            Homs_total <- sum(jLoc_kDeme==0) + sum(jLoc_kDeme==2)
            Hom_obs <- sum(sum(jLoc_kDeme==0) + sum(jLoc_kDeme==2))/(sum(jLoc_kDeme==0) + sum(jLoc_kDeme==1) + sum(jLoc_kDeme==2))
            Het_obs <- sum(jLoc_kDeme==1) /(sum(jLoc_kDeme==0) + sum(jLoc_kDeme==1) + sum(jLoc_kDeme==2))
            Het_exp <- 2*p_freq*q_freq
            F_coeff <- 1-Het_exp
            # send to data matrix
            DemeMatrix2D_p[thisY,thisX] <- p_freq
            DemeMatrix2D_q[thisY,thisX] <- q_freq
            DemeMatrix2D_Ngenes[thisY,thisX] <- Ngenes
            DemeMatrix2D_pCount[thisY,thisX] <- pCount
            DemeMatrix2D_qCount[thisY,thisX] <- qCount
            # genotype counts
            DemeMatrix2D_ppCount[thisY,thisX] <-  ppCount
            DemeMatrix2D_qqCount[thisY,thisX] <-  qqCount
            DemeMatrix2D_pqCount[thisY,thisX] <-  Hets
            # genotype freqs
            DemeMatrix2D_pp[thisY,thisX] <-  ppCount/TotalGenotyped
            DemeMatrix2D_qq[thisY,thisX] <-  qqCount/TotalGenotyped
            DemeMatrix2D_pq[thisY,thisX] <-  Hets/TotalGenotyped
            # record individuals
            DemeMatrix2D_inds[thisY,thisX] <- paste(as.vector(inds_kDeme),collapse="/")
            # Genotypes
            DemeMatrix2D_genos[thisY,thisX] <- paste(jLoc_kDeme,collapse=",")
          }
        }
      }
    }
    # write.csv(DemeIndividuals,"DemeIndividuals.csv")
    
    # individuals
    DemeIndividuals_melted <- reshape2::melt(t(DemeIndividuals))
    # allele freq tables
    DemeMatrix2D_p_melted <- reshape2::melt(t(DemeMatrix2D_p))
    colnames(DemeMatrix2D_p_melted) <- c("X","Y","p")
    DemeMatrix2D_q_melted <- reshape2::melt(t(DemeMatrix2D_q))
    colnames(DemeMatrix2D_q_melted) <- c("X","Y","q")
    # allele count tables
    DemeMatrix2D_pCount_melted <- reshape2::melt(t(DemeMatrix2D_pCount))
    colnames(DemeMatrix2D_pCount_melted) <- c("X","Y","pCount")
    DemeMatrix2D_qCount_melted <- reshape2::melt(t(DemeMatrix2D_qCount))
    colnames(DemeMatrix2D_qCount_melted) <- c("X","Y","qCount")
    # genotype count tables
    DemeMatrix2D_ppCount_melted <- reshape2::melt(t(DemeMatrix2D_ppCount))
    colnames(DemeMatrix2D_ppCount_melted) <- c("X","Y","ppCount")
    DemeMatrix2D_qqCount_melted <- reshape2::melt(t(DemeMatrix2D_qqCount))
    colnames(DemeMatrix2D_qqCount_melted) <- c("X","Y","qqCount")
    DemeMatrix2D_pqCount_melted <- reshape2::melt(t(DemeMatrix2D_pqCount))
    colnames(DemeMatrix2D_pqCount_melted) <- c("X","Y","pqCount")
    # genotype freq tables
    DemeMatrix2D_pp_melted <- reshape2::melt(t(DemeMatrix2D_pp))
    colnames(DemeMatrix2D_pp_melted) <- c("X","Y","pp")
    DemeMatrix2D_qq_melted <- reshape2::melt(t(DemeMatrix2D_qq))
    colnames(DemeMatrix2D_qq_melted) <- c("X","Y","qq")
    DemeMatrix2D_pq_melted <- reshape2::melt(t(DemeMatrix2D_pq))
    colnames(DemeMatrix2D_pq_melted) <- c("X","Y","pq")
    
    # Inds
    DemeMatrix2D_inds_melted <- reshape2::melt(t(DemeMatrix2D_inds))
    colnames(DemeMatrix2D_inds_melted) <- c("X","Y","inds")
    # Genotypes
    DemeMatrix2D_genos_melted <- reshape2::melt(t(DemeMatrix2D_genos))
    colnames(DemeMatrix2D_genos_melted) <- c("X","Y","genos")
    # nrow(DemeMatrix2D_pq_melted); nrow(DemeMatrix2D_inds_melted); nrow(DemeMatrix2D_genos_melted)
    # convert Easting to start from 0 (+1000). This starting at 1000 is to avoid tricky intercept problems for far West demes
    X_adj <- (DemeMatrix2D_pCount_melted[,"X"])-min(Xboundaries_midpoint)
    Y_adj <- (DemeMatrix2D_pCount_melted[,"Y"])-min(Yboundaries_midpoint)
    # combine
    DemeMatrix2D_combine <- cbind(DemeMatrix2D_pCount_melted[,1:2],X_adj,Y_adj,DemeMatrix2D_pCount_melted$pCount,
                                  DemeMatrix2D_qCount_melted$qCount,DemeMatrix2D_p_melted$p,DemeMatrix2D_q_melted$q,
                                  DemeMatrix2D_ppCount_melted$ppCount,DemeMatrix2D_qqCount_melted$qqCount,
                                  DemeMatrix2D_pqCount_melted$pqCount,DemeMatrix2D_pp_melted$pp,DemeMatrix2D_qq_melted$qq,
                                  DemeMatrix2D_pq_melted$pq,DemeMatrix2D_inds_melted$inds,DemeMatrix2D_genos_melted$genos)
    colnames(DemeMatrix2D_combine) <- c("X","Y","X_adj","Y_adj","pCount","qCount","p","q","ppCount","qqCount","pqCount","pp","qq","pq","inds","genos")
    # just in case the demes with no data vary among loci - we find the minX now.
    xMin <- min(DemeMatrix2D_combine[,"X_adj"])
    xMax <- max(DemeMatrix2D_combine[,"X_adj"])
    yMin <- min(DemeMatrix2D_combine[,"Y_adj"])
    yMax <- max(DemeMatrix2D_combine[,"Y_adj"])
    
    DemeMatrix2D_combine_noZero <- DemeMatrix2D_combine[DemeMatrix2D_combine[,"pCount"]!=-9,]
    # DemeMatrix2D_combine_noZero[,1:10]
    # sum(DemeMatrix2D_combine_noZero[,"pCount"] + DemeMatrix2D_combine_noZero[,"qCount"])/2
    # head(DemeMatrix2D_combine_noZero)
    # head(DemeMatrix2D_pqCount)
    # min(DemeMatrix2D_combine_noZero[,"X_adj"])
    # DemeMatrix2D_combine[DemeMatrix2D_combine[,"pp"]!=-9,c(1:10)]
    
    ##########################
    # Pie plot with transect #
    ##########################
    # noPlots=F

    if (noPlots==F) {
      if (length(lociNames)>=4) {
        # now set up multipanel
        par(mfrow=c(2,1))
        par(mar=c(3,3,3,3))
        par(oma=c(3,3,3,3))
      }
      if (length(lociNames)<4) {
        # now set up multipanel
        par(mfrow=c(2,1))
        par(mar=c(3,3,3,3))
        par(oma=c(3,3,3,3))
      } 
      #if (plotTerrain==T) {
      #  par(mfrow=c(1,1))
      #  par(mar=c(5.5,4.5,4,2) + 0.1) 
      #} 
      # note, can only calculate min and max once. Because if we exclude demes with <minSample, and update,
      # we would shift all the x and y positions of demes. Best to keep originals, because we dont know in advance
      # how the sample sizes vary across loci.
      # xMin <- min(DemeMatrix2D_combine_noZero[,"X_adj"])
      # xMax <- max(DemeMatrix2D_combine_noZero[,"X_adj"])
      # yMin <- min(DemeMatrix2D_combine_noZero[,"Y_adj"])
      # yMax <- max(DemeMatrix2D_combine_noZero[,"Y_adj"])
      chartLabel <- markerLabelPieChart(thisLocusName,SNPlociListUpdated)
      
      # plot(DEM_merged,main=locusToPlot,ylim=yRange,xlim=xRange,alpha=alphaTerm,
      #      axes=FALSE,yaxt="n",xaxt="n",ann=F,cex=2,
      #      col=colPallete(colPalleteNum))
      
      plot(DemeMatrix2D_combine_noZero[,"X_adj"],DemeMatrix2D_combine_noZero[,"Y_adj"],type='n',
           main="",cex.main=1,bty="n",yaxt="n",xaxt="n",
           xlab = "",ylab = "",col = "black", bg = "black",xaxs = "i", yaxs = "i",
           xlim=c(xMin-1000,24300),
           ylim=c(yMin-1000,10500))
      axis(1,at=seq(xMin,22000,1000),labels=NA,col.axis="black",las=2,cex.axis=1.1,lwd=1.6,pos=yMin, tck = -.02)
      axis(1,at=seq(xMin,22000,2000),labels=seq(xMin,22000,2000)/1000,col.axis="black",tck=-0.3,las=1,cex.axis=1.1,lwd=0,pos=yMin)
      axis(2,at=seq(yMin,yMax+1000,1000),labels=NA,col.axis="black",las=1,cex.axis=1.1,lwd=1.6,pos=xMin,tck = -.02)
      axis(2,at=seq(yMin,yMax+1000,2000),labels=seq(yMin,yMax+1000,2000)/1000,col.axis="black",tck=-.03,las=1,cex.axis=1.1,lwd=0,pos=xMin-0.5)
      xMinLine <- min(seq(xMin,22000,2000))
      xMaxLine <- max(seq(xMin,22000,2000))
      yMinLine <- min(seq(yMin,yMax+1000,1000))
      yMaxLine <- max(seq(yMin,yMax+1000,1000))
      lines(x=c(xMaxLine,xMaxLine),y=c(yMinLine,yMaxLine),lwd=1.6)
      lines(x=c(xMinLine,xMaxLine),y=c(yMaxLine,yMaxLine),lwd=1.6)
      title(paste(chartLabel,paste(demeSize,"m demes",sep=""),sep=", "),line=2.25,cex.main=0.9)
      mtext("East (km) / Easting",side=1,cex=1,line=2.2)
      mtext("North (km) / Northing",side=2,cex=1,line=2.2)
      curve(((transectGradient*x)+(transectIntercept)), from=c(min(DemeMatrix2D_combine_noZero[,"X_adj"]-300),to=max(DemeMatrix2D_combine_noZero[,"X_adj"])), col="darkgray",add = TRUE,type = "l",lwd=2,lty=1)
      
      # Easting and Northing axis
      # head(DemeMatrix2D_combine_noZero)
      #min(DemeMatrix2D_combine_noZero[,"X"])
      #max(DemeMatrix2D_combine_noZero[,"X"])
      #min(DemeMatrix2D_combine_noZero[,"Y"])
      #max(DemeMatrix2D_combine_noZero[,"Y"])
      #min(DemeMatrix2D_combine_noZero[,"Y_adj"])
      #max(DemeMatrix2D_combine_noZero[,"Y_adj"])
      # DemeMatrix2D_combine_noZero[DemeMatrix2D_combine_noZero[,"Y_adj"]==min(DemeMatrix2D_combine_noZero[,"Y_adj"]),1:6]
      # DemeMatrix2D_combine_noZero[DemeMatrix2D_combine_noZero[,"Y_adj"]==max(DemeMatrix2D_combine_noZero[,"Y_adj"]),1:6]
      if (plotEastNorthing==T) {
        axis(1,at=seq(-435,24000,2000),labels=NA,outer=F,col.axis="black",las=2,cex.axis=1.1,lwd=1.6,pos=yMaxLine, tck = 0.02)
        axis(1,at=seq(-435,24000,2000),labels=seq(410000,435435,2000),col.axis="black",tck=-0.3,las=1,cex.axis=1.1,lwd=0,pos=yMaxLine+3000)
        axis(2,at=seq(-2246,10000,2000),labels=NA,col.axis="black",las=1,cex.axis=1.1,lwd=1.6,pos=xMaxLine,tck = 0.02)
        axis(2,at=seq(-2246,10000,2000),labels=seq(4682000,4694000,2000),col.axis="black",tck=-.03,las=1,cex.axis=1.1,lwd=0,pos=xMaxLine+3500)
      }
      if (addPie==T) {
        # plotting pie charts
        # min(DemeMatrix2D_combine_noZero[,"X_adj"])
        totalN <- (DemeMatrix2D_combine_noZero[,"pCount"]+DemeMatrix2D_combine_noZero[,"qCount"])/2
        #minSample <- 2
        keep <- which(totalN>=minSample)
        # this part causing the discrepancy. Once we remove demes with < minSample, the min X_adj is shifted. 
        # If we update the min X_adj then the two sets of values are no longer consistent.
        DemeMatrix2D_combine_noZero <- DemeMatrix2D_combine_noZero[keep,]
        allthisP <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[,"pCount"]))
        allthisQ <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[,"qCount"]))
        allSum <- allthisP+allthisQ
        allFreq <- allthisP/allSum
        if (mean(allFreq[1:4])>0.6) {
          allthisP_temp <- allthisQ
          allthisQ_temp <- allthisP
          allthisQ <- allthisP_temp
          allthisP <- allthisQ_temp
        }
        for (thisDeme2D in 1:nrow(DemeMatrix2D_combine_noZero)) {
          # thisDeme2D <- 1
          # head(DemeMatrix2D_combine_noZero)
          if (plotTerrain==F) {
            thisXraw <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[thisDeme2D,"X_adj"]))
            thisYraw <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[thisDeme2D,"Y_adj"]))
          }
          if (plotTerrain==T) {
            thisXraw <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[thisDeme2D,"X"]))
            thisYraw <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[thisDeme2D,"Y"]))
          }
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
          add.pie(z=c(thisP,thisQ), x=thisXraw, y=thisYraw, radius=120,col=c((rgb(255,255,0,255,maxColorValue=255)),(rgb(255,0,0,255,maxColorValue=255))),labels="")
          #}
        }
      }
      if (addPie==F) {
        points(DemeMatrix2D_combine_noZero[,"X_adj"],DemeMatrix2D_combine_noZero[,"Y_adj"])
      }
    }
   
    #######################################
    # Calculating distance along transect #
    #######################################
    # noPlots <- T
    # now find the intersection to transect and how far along the transect each deme sits
    distTransect <- matrix(0,nrow(DemeMatrix2D_combine_noZero),2)
    colnames(distTransect) <- c("DistanceRefaxis","DistanceRefaxis_adj")
    
    # could be a problem with DistanceRefaxis_adj, is it consistent across loci? We cant have different adjustments
    # creep in at each locus because the number of samples differs a little. This may shift which demes have sufficient
    # numbers among loci and influence the adjustment.
    # calculate distance along transect
    xCross_start_set <- NULL
    for (thisRow in 1:nrow(distTransect)) {
      # thisRow <- 1
      # 1. find perpendicular linear equation
      # y - y1 = m(x - x1)
      # head(DemeMatrix2D_combine_noZero)
      DemeMatrix2D_combine_noZero[,"X_adj"]
      thisY <- DemeMatrix2D_combine_noZero[thisRow,"Y_adj"]
      thisX <- DemeMatrix2D_combine_noZero[thisRow,"X_adj"]
      # interceptPerp <- ((4.16*(-thisX))+thisY)
      # eqnPerp <- (4.16*x) + 4788100
      #  perpendicular gradient can be found with: -(1/x)
      # find point of contact with reference line
      # y values must be equal, therefore
      # mx + b = mx + b
      # original
      # y = -(0.24*420000) + 4788100
      # perpendicular
      # y = (4.16*420000) + 4788100
      # difference <- (-(0.24*420000) + 4788100) - ((4.16*420000) + 4788100)
      # 4788100 - difference
      # Position and distance along reference line. Using basic trig.
      # 1. Perpedicular gradient = -(1/x)
      gradientPerp <- (1/transectGradient)*-1
      # 2. Finding linear equation of a line: y-y1 = m(x-x1)
      # y = mx - mx1 + y1
      interceptPerp <- ((gradientPerp*(-thisX))+thisY)
      # 3. Finding where two lines intersect.
      # Set them to be equal. Knowing b and b1, and setting y equal we only have to solve for x
      # i.e. mx + b = m1x + b1
      # mx - m1x = b1 - b
      # therefore, x = (b-b1)/[-(m-m1)]
      xCross <- (interceptPerp-transectIntercept)/((gradientPerp-transectGradient)*-1)
      yCross <- (transectGradient*xCross)+transectIntercept
      if (thisRow == 1) {
        xCross_start <- xCross
        yCross_start <- yCross
        xCross_start_set <- xCross_start
        yCross_start_set <- yCross_start
      }
      # Lastly, determine how far is the intersection point along the reference line (compared to x=0 and y=intercept)
      #xvals <- c(-1000,xCross)
      #yvals <- c(((transectGradient*-1000)+transectIntercept),yCross)
      xvals <- c(xCross_start,xCross)
      yvals <- c(((transectGradient*xCross_start)+transectIntercept),yCross)
      distTransect[thisRow,2] <- as.vector(dist(cbind(xvals,yvals),method="euclidean"))
      #curve(((transectGradient*x)+(transectIntercept)), from=min(DemeMatrix2D_combine_noZero[,"X_adj"]),
      #      to=max(DemeMatrix2D_combine_noZero[,"X_adj"]), col="darkgray",add = TRUE,type = "l",lwd=3,lty=1)
      # An alternative is to use the solve function, gives the same result
      # linearEqn <- matrix(c(transectIntercept,transectGradient,interceptPerp,gradientPerp),2,2,byrow=T)
      # XYcross <- c(-solve(cbind(linearEqn[,2],-1)) %*% linearEqn[,1])
      # plotPerp<- T
      # noPlots <- F
      if (noPlots==F) {
        if (plotPerp==T) {
          lines(x=c(thisX,xCross),y=c(thisY,yCross),lwd=0.8,col='blue')
          #points(x=xCross,y=yCross,pch=21,col='red',cex=0.6,lwd=0.5)
        }
      }
    }
    # original in column 1
    distTransect[,1] <- distTransect[,2] 
    # rescale to min val (may need to keep an eye on this and make it further negative , <5000. Depends on where the perp hits)
    #distTransect[,2] <- distTransect[,2] - min(distTransect[,2])
    # head(DemeMatrix2D_combine_noZero)
    #######################
    # Pool Seq comparison #
    #######################
    # only do this once
    if (poolSeqPlot==T & thisLocus==1) {
        # setup matrix for poolSeq comparison
        # bring in KASP geno
          poolSeqStats[1,4:5] <- DemeMatrix2D_combine_noZero[1,1:2]
          poolSeqStats[1,6:7] <- DemeMatrix2D_combine_noZero[1,3:4]
          Deme2D_pools <- DemeMatrix2D
          Deme2D_pools[,] <- -9
        
        demeCounter <- 0
        theseRows2D <- NULL
        for (thisXcol in 1:ncol(Deme2D_pools)) {
          for (thisYrow in 1:nrow(Deme2D_pools)) {
            demeCounter <- demeCounter + 1
            # thisXcol <- 26
            # thisYrow <- 17
            thisDeme_mid_X <- as.numeric(colnames(Deme2D_pools)[thisXcol])
            thisDeme_mid_Y <- as.numeric(rownames(Deme2D_pools)[thisYrow])
            theseX <- c(as.numeric(colnames(Deme2D_pools)[thisXcol])-(demeSize/2),as.numeric(colnames(Deme2D_pools)[thisXcol])+(demeSize/2))
            theseY <- c(as.numeric(rownames(Deme2D_pools)[thisYrow])-(demeSize/2),as.numeric(rownames(Deme2D_pools)[thisYrow])+(demeSize/2))
            theseRows2D <- as.numeric(poolSeqStats[,"Easting"]) >= theseX[1] & as.numeric(poolSeqStats[,"Easting"]) <= theseX[2] &
              as.numeric(poolSeqStats[,"Northing"]) <= theseY[2] & as.numeric(poolSeqStats[,"Northing"]) >= theseY[1]
            #jLoc_kDeme[2] <- "NGY"
            #(all(jLoc_kDeme!="NGY"))
            #HZdataSetXYdata[,colnames(HZdataSetXYdata)==thisLocusName]
            #xCross_start_set <- poolSeqStats[1,"X_adj"]
            #yCross_start_set <- poolSeqStats[1,"Y_adj"]
            
            if (any(theseRows2D[2:7])==T) {
              #print ("Found one")
              #print (paste0("thisXcol: ",thisXcol))
              #print (paste0("thisYrow: ",thisYrow))
              poolSeqStats[which(theseRows2D==T),"XPos_deme"] <- thisDeme_mid_X
              poolSeqStats[which(theseRows2D==T),"YPos_deme"] <- thisDeme_mid_Y
              # Need to use the most upper left deme from KASP data as starting reference
              # the Easting & Northing was adjusted by min values as follows:
              # poolSeqStats[1,"XPos"]-min(Xboundaries_midpoint)
              # poolSeqStats[1,"YPos"]-min(Yboundaries_midpoint)
              poolSeqStats[which(theseRows2D==T),6] <- poolSeqStats[which(theseRows2D==T),4] - min(Xboundaries_midpoint)
              poolSeqStats[which(theseRows2D==T),7] <- poolSeqStats[which(theseRows2D==T),5] - min(Yboundaries_midpoint)
              
              # thisRow <- 2
              # 1. find perpendicular linear equation
              # y - y1 = m(x - x1)
              # interceptPerp <- ((4.16*(-thisX))+thisY)
              # eqnPerp <- (4.16*x) + 4788100
              #  perpendicular gradient can be found with: -(1/x)
              # find point of contact with reference line
              # y values must be equal, therefore
              # mx + b = mx + b
              # original
              # y = -(0.24*420000) + 4788100
              # perpendicular
              # y = (4.16*420000) + 4788100
              # difference <- (-(0.24*420000) + 4788100) - ((4.16*420000) + 4788100)
              # 4788100 - difference
              # Position and distance along reference line. Using basic trig.
              # 1. Perpendicular gradient = -(1/x)
              gradientPerp <- (1/transectGradient)*-1
              # 2. Finding linear equation of a line: y-y1 = m(x-x1)
              # y = mx - mx1 + y1
              interceptPerp <- ((gradientPerp*(-poolSeqStats[which(theseRows2D==T),"X_adj"]))+poolSeqStats[which(theseRows2D==T),"Y_adj"])
              # 3. Finding where two lines intersect.
              # Set them to be equal. Knowing b and b1, and setting y equal we only have to solve for x
              # i.e. mx + b = m1x + b1
              # mx - m1x = b1 - b
              # therefore, x = (b-b1)/[-(m-m1)]
              xCross <- (interceptPerp-transectIntercept)/((gradientPerp-transectGradient)*-1)
              yCross <- (transectGradient*xCross)+transectIntercept
              # yCross_start and xCross_start (first X pos of deme with data) defined with KASP sample set above
              # Lastly, determine how far is the intersection point along the reference line (compared to x=0 and y=intercept)
              #xvals <- c(-1000,xCross)
              #yvals <- c(((transectGradient*-1000)+transectIntercept),yCross)
              xvals <- c(xCross_start_set,xCross)
              yvals <- c(((transectGradient*xCross_start_set)+transectIntercept),yCross)
              
              # error could also be occuring here for 200m demes
              # cat("\n euclidean distance test 1: ",as.vector(dist(cbind(xvals,yvals),method="euclidean")))
              # flush.console()
              # print("poolSeqStats table")
              # print(poolSeqStats)
               
              # cat("\n theseRows2D==T: ",as.vector(which(theseRows2D==T)))
              # flush.console()
              
               #poolSeqStats[which(theseRows2D==T),"DistAlongTransect"] <- as.vector(dist(cbind(xvals,yvals),method="euclidean"))
              #print (paste0("thisDistAlongTransect: ", poolSeqStats[which(theseRows2D==T),"DistAlongTransect"]))
              #curve(((transectGradient*x)+(transectIntercept)), from=min(DemeMatrix2D_combine_noZero[,"X_adj"]),
              #      to=max(DemeMatrix2D_combine_noZero[,"X_adj"]), col="darkgray",add = TRUE,type = "l",lwd=3,lty=1)
              # An alternative is to use the solve function, gives the same result
              # linearEqn <- matrix(c(transectIntercept,transectGradient,interceptPerp,gradientPerp),2,2,byrow=T)
              # XYcross <- c(-solve(cbind(linearEqn[,2],-1)) %*% linearEqn[,1])
              if (noPlots==F) { 
                if (plotPerp==T) {
                  lines(x=c(poolSeqStats[which(theseRows2D==T),"X_adj"],xCross),y=c(poolSeqStats[which(theseRows2D==T),"Y_adj"],yCross),lwd=1.1,col='black')
                  points(x=poolSeqStats[which(theseRows2D==T),"X_adj"],y=poolSeqStats[which(theseRows2D==T),"Y_adj"],pch=23,col='black',bg=c(rgb(255,255,255,230,maxColorValue=255)),lwd=1.5,cex=1.6)
                  #points(x=xCross,y=yCross,pch=21,col='black',bg='black',cex=0.6,lwd=0.5)
                }
              }
            }
          }
        }
      }
    # poolSeq within bounds of KASP genotypes? Yes
    # min(poolSeqStats[,"XPos_deme"])>minX; max(poolSeqStats[,"XPos_deme"])<maxX 
    # min(poolSeqStats[,"YPos_deme"])>minY; max(poolSeqStats[,"YPos_deme"])<maxY 
    
    ################################################
    # plot allele freqs by distance along transect #
    ################################################
    SampleSize <- (DemeMatrix2D_combine_noZero[,"pCount"]+ DemeMatrix2D_combine_noZero[,"qCount"])/2
    SampleSizeMin <- SampleSize > minSample
    DemeMatrix2D_combine_noZero <- cbind(DemeMatrix2D_combine_noZero[,1:4],round(distTransect[,2],1),DemeMatrix2D_combine_noZero[,5:ncol(DemeMatrix2D_combine_noZero)])
    colnames(DemeMatrix2D_combine_noZero)[5] <- "distAlongTransect"
      # re-orientate p (magenta) & q (yellow)
    if (mean(DemeMatrix2D_combine_noZero[1:4,"p"])
        > mean(DemeMatrix2D_combine_noZero[(nrow(DemeMatrix2D_combine_noZero)-4):nrow(DemeMatrix2D_combine_noZero),"p"])) {
        DemeMatrix_temp <- DemeMatrix2D_combine_noZero
        DemeMatrix_temp[,"q"] <- DemeMatrix2D_combine_noZero[,"p"]
        DemeMatrix_temp[,"p"] <- DemeMatrix2D_combine_noZero[,"q"]
        DemeMatrix_temp[,"qCount"] <- DemeMatrix2D_combine_noZero[,"pCount"]
        DemeMatrix_temp[,"pCount"] <- DemeMatrix2D_combine_noZero[,"qCount"]
        DemeMatrix_temp[,"ppCount"] <- DemeMatrix2D_combine_noZero[,"qqCount"]
        DemeMatrix_temp[,"qqCount"] <- DemeMatrix2D_combine_noZero[,"ppCount"]
        DemeMatrix_temp[,"pp"] <- DemeMatrix2D_combine_noZero[,"qq"]
        DemeMatrix_temp[,"qq"] <- DemeMatrix2D_combine_noZero[,"pp"]
        DemeMatrix2D_combine_noZero <- DemeMatrix_temp
      }
      
    SampleSizeMin <- SampleSize >= minSample
    totalSampleSizeMinDemes <- sum(SampleSize[SampleSizeMin])
    textSec <- paste(paste("demes=",sum(SampleSizeMin),sep=""),
                       paste("n=",totalSampleSizeMinDemes,sep=""),sep=", ")
      #xMin <- min(DemeMatrix2D_combine_noZero[,"X_adj"])
      #xMax <- max(DemeMatrix2D_combine_noZero[,"X_adj"])
      #yMin <- min(DemeMatrix2D_combine_noZero[,"Y_adj"])
      #yMax <- max(DemeMatrix2D_combine_noZero[,"Y_adj"])
    if (noPlots==F) {
      par(mfrow=c(1,1))
      par(mar=c(3,3,3,3))
      par(oma=c(3,3,3,3))
      plot(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
           DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],bg = "black",xaxs = "i", yaxs = "i",
             xlab="",ylab="",main="",xlim=c(xMin-1000,22000),
           ylim=c(-0.1,1.1),type='n',cex.main=0.9,bty="n",yaxt="n",xaxt="n")
      axis(1,at=seq(xMin,22000,1000),labels=NA,col.axis="black",las=2,cex.axis=1.1,lwd=1.6,pos=-0.01, tck = -.02)
      axis(1,at=seq(xMin,22000,2000),labels=seq(xMin,22000,2000)/1000,col.axis="black",tck=-0.3,las=1,cex.axis=1.1,lwd=0,pos=-0.02)
      axis(2,at=seq(-0.01,1.01,0.1),labels=NA,col.axis="black",las=1,cex.axis=1.1,lwd=1.6,pos=xMin,tck = -.02)
      axis(2,at=seq(0,1.0,0.2),labels=seq(0,1,0.2),col.axis="black",tck=-.03,las=1,cex.axis=1.1,lwd=0,pos=xMin)
      title(paste(chartLabel,textSec,sep="\n"),line=0.1,cex.main=0.9)
      lines(x=c(22000,22000),y=c(0,1),lwd=1.6)
      lines(x=c(xMinLine,22000),y=c(1,1),lwd=1.6)
      # mark centre of phenotypic cline
      #arrows(x0 = 12580,y0 = -0.22,x1 = 12580,y1 = -0.015,col='red',lwd=2)
      
      if (plotAxisLab==T) {
          mtext("Distance along transect (km)",side=1,cex=1,line=2.2)
          mtext("Allele frequency (p)",side=2,cex=1,line=2.5)
      }
    
      # Some poolSeq comparisons
      # allele frequences drawn as lines with no points per deme (too many loci to show points)
      if (poolSeqPlot==T) {
        if (length(poolSeqClines)>0) {
          for (thisRow in 1:nrow(poolSeqClines)) {
            # thisRow <- 1
            poolSeqClineThisLocus <- cbind(poolSeqStats[2:7,"pool"],poolSeqStats[2:7,"DistAlongTransect"],as.vector(unlist(poolSeqClines[thisRow,9:14])))
            colnames(poolSeqClineThisLocus) <- c("pool","DistAlongTransect","p")
            poolSeqClineThisLocus <- as.data.frame(poolSeqClineThisLocus)
            poolSeqClineThisLocus[,"DistAlongTransect"] <- as.numeric(poolSeqClineThisLocus[,"DistAlongTransect"])
            poolSeqClineThisLocus[,"p"] <- as.numeric(poolSeqClineThisLocus[,"p"])
            poolSeqClineThisLocus <- poolSeqClineThisLocus[order(poolSeqClineThisLocus[,"DistAlongTransect"]),]
            
            # background loci not linked to any colour loci
            if (poolSeqClines[thisRow,"ColourGene"]==-9) {
              back_col <- colourTable[colourTable[,"Region"]=="background","lineColour"]
              back_col_transparent <- adjustcolor(back_col, alpha.f = 0.3)
              lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=back_col_transparent)
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
              lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=colourTable[colourTable[,"Region"]=="RUBIAlinked","lineColour"])
            }
            if (poolSeqClines[thisRow,"ColourGene"]=="RUBIA" & poolSeqClines[thisRow,"distanceToGene"]=="inGene") {
              lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=colourTable[colourTable[,"Region"]=="RUBIA","lineColour"])
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
          if (plotKASPrep==T) {
            for (thisLoc in lociNames) {
              # Cremosa KASP
              if (thisLoc == "s1187_290152") {
                # SNP KASP points
                points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                       DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="CREMOSA","pointColour"],lwd=1.5)
              
                # cbind(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],DemeMatrix2D_combine_noZero[SampleSizeMin,"pCount"],DemeMatrix2D_combine_noZero[SampleSizeMin,"qCount"],DemeMatrix2D_combine_noZero[SampleSizeMin,"p"])
                }
              # Flavia KASP down
              if (thisLoc == "s316_93292" ) {
                # SNP KASP points
                points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                       DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="FLAVIA","pointColour"],lwd=1.5)
              }
              # Flavia KASP up
              if (thisLoc == "s316_257789") {
                # SNP KASP points
                points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                       DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="FLAVIA","pointColour"],lwd=1.5)
              }
              # SULF KASP
              if (thisLoc == "s91_39699") {
                # SNP KASP points
                points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                       DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="SULF","pointColour"],lwd=1.5)
              }
              # RUBIA KASP
              if (thisLoc == "s261_720757") {
                # SNP KASP points
                points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                       DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="RUBIA","pointColour"],lwd=1.5)
              }
              # ROS el KASP
              if (thisLoc == "ros_assembly_543443" & FigName!="Chr6") {
                # SNP KASP points
                points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                       DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="ROS","pointColour"],lwd=1.5)
              }
              
              # highlight key ROS1 locus in black
              #if (thisLoc == "ros_assembly_543443" & FigName=="Chr6") {
                
              if (thisLoc == "ros_assembly_543443" & FigName=="Chr6") {
                poolSeqClineThisLocus <- cbind(poolSeqStats[2:7,"pool"],poolSeqStats[2:7,"DistAlongTransect"],as.vector(unlist(poolSeqClines[poolSeqClines[,"position"]==52889078,c(9:14)])))
                colnames(poolSeqClineThisLocus) <- c("pool","DistAlongTransect","p")
                poolSeqClineThisLocus <- as.data.frame(poolSeqClineThisLocus)
                poolSeqClineThisLocus[,"DistAlongTransect"] <- as.numeric(poolSeqClineThisLocus[,"DistAlongTransect"])
                poolSeqClineThisLocus[,"p"] <- as.numeric(poolSeqClineThisLocus[,"p"])
                poolSeqClineThisLocus <- poolSeqClineThisLocus[order(poolSeqClineThisLocus[,"DistAlongTransect"]),]
                points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                       DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="ROS","pointColour"],lwd=1.5)
                lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=2.5,col='black')
                points(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],pch=23,col='black',bg=c(rgb(255,255,255,230,maxColorValue=255)),lwd=1.5,cex=1.6)
                
              }
            }
          }
        }
      }
    } 
    # pdf version of above double figure
    # now set up multipanel
    fileDemes <- paste(paste(thisLocusName,"_Transect_PoolSeq_vs_KASP",FigName,sep=""),".pdf",sep="_")
    pdf(fileDemes,width = 14.5, height = 8.2,pointsize=18)
    par(mfrow=c(1,1))
    par(mar=c(3,3,3,3))
    par(oma=c(3,3,3,3))    
    #chartLabel <- markerLabelPieChart(thisLocusName,SNPlociListUpdated)
    
    # plot(DEM_merged,main=locusToPlot,ylim=yRange,xlim=xRange,alpha=alphaTerm,
    #      axes=FALSE,yaxt="n",xaxt="n",ann=F,cex=2,
    #      col=colPallete(colPalleteNum))
    plot(DemeMatrix2D_combine_noZero[,"X_adj"],DemeMatrix2D_combine_noZero[,"Y_adj"],type='n',
         main="",cex.main=1,bty="n",yaxt="n",xaxt="n",
         xlab = "",ylab = "",col = "black", bg = "black",xaxs = "i", yaxs = "i",
         xlim=c(xMin-1000,24300),
         ylim=c(yMin-1000,10500))
    axis(1,at=seq(xMin,xMax,1000),labels=NA,col.axis="black",las=2,cex.axis=1.1,lwd=1.6,pos=yMin, tck = -.02)
    axis(1,at=seq(xMin,xMax,2000),labels=seq(xMin,xMax,2000)/1000,col.axis="black",tck=-0.3,las=1,cex.axis=1.1,lwd=0,pos=yMin)
    axis(2,at=seq(yMin,yMax+1000,1000),labels=NA,col.axis="black",las=1,cex.axis=1.1,lwd=1.6,pos=xMin,tck = -.02)
    axis(2,at=seq(yMin,yMax+1000,2000),labels=seq(yMin,yMax+1000,2000)/1000,col.axis="black",tck=-.03,las=1,cex.axis=1.1,lwd=0,pos=xMin-0.5)
    xMinLine <- min(seq(xMin,xMax,2000))
    xMaxLine <- max(seq(xMin,xMax,2000))
    yMinLine <- min(seq(yMin,yMax+1000,1000))
    yMaxLine <- max(seq(yMin,yMax+1000,1000))
    lines(x=c(xMaxLine,xMaxLine),y=c(yMinLine,yMaxLine),lwd=1.6)
    lines(x=c(xMinLine,xMaxLine),y=c(yMaxLine,yMaxLine),lwd=1.6)
    secondLevelHeader <- paste0("demeSize= ",demeSize,", gradient= ",round(transectGradient,3),", ","transect= ", transectIntercept)
    chartLabelFinal <- paste(thisLocusName,secondLevelHeader,sep="; ")
    # run details
    #title(paste(chartLabel,paste(demeSize,"m demes",sep=""),sep=", "),line=2.25,cex.main=0.9)
    title(chartLabelFinal,line=2.25,cex.main=0.9)
    mtext("East (km) / Easting",side=1,cex=1,line=2.2)
    mtext("North (km) / Northing",side=2,cex=1,line=2.2)
    curve(((transectGradient*x)+(transectIntercept)), from=c(min(DemeMatrix2D_combine_noZero[,"X_adj"]-300),to=max(DemeMatrix2D_combine_noZero[,"X_adj"])), col="darkgray",add = TRUE,type = "l",lwd=2,lty=1)
    
    # Easting and Northing axis
    # head(DemeMatrix2D_combine_noZero)
    #min(DemeMatrix2D_combine_noZero[,"X"])
    #max(DemeMatrix2D_combine_noZero[,"X"])
    #min(DemeMatrix2D_combine_noZero[,"Y"])
    #max(DemeMatrix2D_combine_noZero[,"Y"])
    #min(DemeMatrix2D_combine_noZero[,"Y_adj"])
    #max(DemeMatrix2D_combine_noZero[,"Y_adj"])
    # DemeMatrix2D_combine_noZero[DemeMatrix2D_combine_noZero[,"Y_adj"]==min(DemeMatrix2D_combine_noZero[,"Y_adj"]),1:6]
    # DemeMatrix2D_combine_noZero[DemeMatrix2D_combine_noZero[,"Y_adj"]==max(DemeMatrix2D_combine_noZero[,"Y_adj"]),1:6]
    if (plotEastNorthing==T) {
      axis(1,at=seq(-435,24000,2000),labels=NA,outer=F,col.axis="black",las=2,cex.axis=1.1,lwd=1.6,pos=yMaxLine, tck = 0.02)
      axis(1,at=seq(-435,24000,2000),labels=seq(410000,435435,2000),col.axis="black",tck=-0.3,las=1,cex.axis=1.1,lwd=0,pos=yMaxLine+3000)
      axis(2,at=seq(-2246,10000,2000),labels=NA,col.axis="black",las=1,cex.axis=1.1,lwd=1.6,pos=xMaxLine,tck = 0.02)
      axis(2,at=seq(-2246,10000,2000),labels=seq(4682000,4694000,2000),col.axis="black",tck=-.03,las=1,cex.axis=1.1,lwd=0,pos=xMaxLine+4000)
    }
    
    if (addPie==T) {
      # plotting pie charts
      # min(DemeMatrix2D_combine_noZero[,"X_adj"])
      totalN <- (DemeMatrix2D_combine_noZero[,"pCount"]+DemeMatrix2D_combine_noZero[,"qCount"])/2
      keep <- which(totalN>=minSample)
      # this part causing the discrepancy. Once we remove demes with < minSample, the min X_adj is shifted. 
      # If we update the min X_adj then the two sets of values are no longer consistent.
      DemeMatrix2D_combine_noZero <- DemeMatrix2D_combine_noZero[keep,]
      allthisP <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[,"pCount"]))
      allthisQ <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[,"qCount"]))
      allSum <- allthisP+allthisQ
      allFreq <- allthisP/allSum
      
      #######################################
      # Calculating distance along transect #
      #######################################
      # noPlots <- T
      # now find the intersection to transect and how far along the transect each deme sits
      distTransect <- matrix(0,nrow(DemeMatrix2D_combine_noZero),2)
      colnames(distTransect) <- c("DistanceRefaxis","DistanceRefaxis_adj")
      
      # could be a problem with DistanceRefaxis_adj, is it consistent across loci? We cant have different adjustments
      # creep in at each locus because the number of samples differs a little. This may shift which demes have sufficient
      # numbers among loci and influence the adjustment.
      # calculate distance along transect
      xCross_start_set <- NULL
      for (thisRow in 1:nrow(distTransect)) {
        # thisRow <- 1
        # 1. find perpendicular linear equation
        # y - y1 = m(x - x1)
        # head(DemeMatrix2D_combine_noZero)
        DemeMatrix2D_combine_noZero[,"X_adj"]
        thisY <- DemeMatrix2D_combine_noZero[thisRow,"Y_adj"]
        thisX <- DemeMatrix2D_combine_noZero[thisRow,"X_adj"]
        # interceptPerp <- ((4.16*(-thisX))+thisY)
        # eqnPerp <- (4.16*x) + 4788100
        #  perpendicular gradient can be found with: -(1/x)
        # find point of contact with reference line
        # y values must be equal, therefore
        # mx + b = mx + b
        # original
        # y = -(0.24*420000) + 4788100
        # perpendicular
        # y = (4.16*420000) + 4788100
        # difference <- (-(0.24*420000) + 4788100) - ((4.16*420000) + 4788100)
        # 4788100 - difference
        # Position and distance along reference line. Using basic trig.
        # 1. Perpedicular gradient = -(1/x)
        gradientPerp <- (1/transectGradient)*-1
        # 2. Finding linear equation of a line: y-y1 = m(x-x1)
        # y = mx - mx1 + y1
        interceptPerp <- ((gradientPerp*(-thisX))+thisY)
        # 3. Finding where two lines intersect.
        # Set them to be equal. Knowing b and b1, and setting y equal we only have to solve for x
        # i.e. mx + b = m1x + b1
        # mx - m1x = b1 - b
        # therefore, x = (b-b1)/[-(m-m1)]
        xCross <- (interceptPerp-transectIntercept)/((gradientPerp-transectGradient)*-1)
        yCross <- (transectGradient*xCross)+transectIntercept
        if (thisRow == 1) {
          xCross_start <- xCross
          yCross_start <- yCross
          xCross_start_set <- xCross_start
          yCross_start_set <- yCross_start
        }
        # Lastly, determine how far is the intersection point along the reference line (compared to x=0 and y=intercept)
        #xvals <- c(-1000,xCross)
        #yvals <- c(((transectGradient*-1000)+transectIntercept),yCross)
        xvals <- c(xCross_start,xCross)
        yvals <- c(((transectGradient*xCross_start)+transectIntercept),yCross)
        distTransect[thisRow,2] <- as.vector(dist(cbind(xvals,yvals),method="euclidean"))
        #curve(((transectGradient*x)+(transectIntercept)), from=min(DemeMatrix2D_combine_noZero[,"X_adj"]),
        #      to=max(DemeMatrix2D_combine_noZero[,"X_adj"]), col="darkgray",add = TRUE,type = "l",lwd=3,lty=1)
        # An alternative is to use the solve function, gives the same result
        # linearEqn <- matrix(c(transectIntercept,transectGradient,interceptPerp,gradientPerp),2,2,byrow=T)
        # XYcross <- c(-solve(cbind(linearEqn[,2],-1)) %*% linearEqn[,1])
        # plotPerp<- T
        # noPlots <- F
        if (noPlots==F) {
          if (plotPerp==T) {
            lines(x=c(thisX,xCross),y=c(thisY,yCross),lwd=0.6,col=c(rgb(0,0,255,180,maxColorValue=255)))
            #points(x=xCross,y=yCross,pch=21,col='red',cex=0.6,lwd=0.5)
          }
        }
      }
      # original in column 1
      distTransect[,1] <- distTransect[,2] 
      # rescale to min val (may need to keep an eye on this and make it further negative , <5000. Depends on where the perp hits)
      #distTransect[,2] <- distTransect[,2] - min(distTransect[,2])
      # head(DemeMatrix2D_combine_noZero)
      
      if (mean(allFreq[1:4])>0.6) {
        allthisP_temp <- allthisQ
        allthisQ_temp <- allthisP
        allthisQ <- allthisQ_temp
        allthisP <- allthisP_temp
      }
      for (thisDeme2D in 1:nrow(DemeMatrix2D_combine_noZero)) {
        # thisDeme2D <- 1
        # head(DemeMatrix2D_combine_noZero)
        if (plotTerrain==F) {
          thisXraw <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[thisDeme2D,"X_adj"]))
          thisYraw <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[thisDeme2D,"Y_adj"]))
        }
        if (plotTerrain==T) {
          thisXraw <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[thisDeme2D,"X"]))
          thisYraw <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[thisDeme2D,"Y"]))
        }
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
        add.pie(z=c(thisP,thisQ), x=thisXraw, y=thisYraw, radius=140,col=c((rgb(255,0,0,200,maxColorValue=255)),(rgb(255,255,0,200,maxColorValue=255))),labels="")
        #}
      }
    }
    if (addPie==F) {
      points(DemeMatrix2D_combine_noZero[,"X_adj"],DemeMatrix2D_combine_noZero[,"Y_adj"])
    }
    
    #######################
    # Pool Seq comparison #
    #######################
    # only do this once
    if (poolSeqPlot==T) {
      # setup matrix for poolSeq comparison
      # bring in KASP geno
      poolSeqStats[1,4:5] <- DemeMatrix2D_combine_noZero[1,1:2]
      poolSeqStats[1,6:7] <- DemeMatrix2D_combine_noZero[1,3:4]
      Deme2D_pools <- DemeMatrix2D
      Deme2D_pools[,] <- -9
      demeCounter <- 0
      theseRows2D <- NULL
      for (thisXcol in 1:ncol(Deme2D_pools)) {
        for (thisYrow in 1:nrow(Deme2D_pools)) {
          demeCounter <- demeCounter + 1
          # thisXcol <- 26
          # thisYrow <- 17
          thisDeme_mid_X <- as.numeric(colnames(Deme2D_pools)[thisXcol])
          thisDeme_mid_Y <- as.numeric(rownames(Deme2D_pools)[thisYrow])
          theseX <- c(as.numeric(colnames(Deme2D_pools)[thisXcol])-(demeSize/2),as.numeric(colnames(Deme2D_pools)[thisXcol])+(demeSize/2))
          theseY <- c(as.numeric(rownames(Deme2D_pools)[thisYrow])-(demeSize/2),as.numeric(rownames(Deme2D_pools)[thisYrow])+(demeSize/2))
          theseRows2D <- as.numeric(poolSeqStats[,"Easting"]) >= theseX[1] & as.numeric(poolSeqStats[,"Easting"]) <= theseX[2] &
            as.numeric(poolSeqStats[,"Northing"]) <= theseY[2] & as.numeric(poolSeqStats[,"Northing"]) >= theseY[1]
          #jLoc_kDeme[2] <- "NGY"
          #(all(jLoc_kDeme!="NGY"))
          #HZdataSetXYdata[,colnames(HZdataSetXYdata)==thisLocusName]
          #xCross_start_set <- poolSeqStats[1,"X_adj"]
          #yCross_start_set <- poolSeqStats[1,"Y_adj"]
          
          if (any(theseRows2D[2:7])==T) {
            #print ("Found one")
            #print (paste0("thisXcol: ",thisXcol))
            #print (paste0("thisYrow: ",thisYrow))
            poolSeqStats[which(theseRows2D==T),"XPos_deme"] <- thisDeme_mid_X
            poolSeqStats[which(theseRows2D==T),"YPos_deme"] <- thisDeme_mid_Y
            # Need to use the most upper left deme from KASP data as starting reference
            # the Easting & Northing was adjusted by min values as follows:
            # poolSeqStats[1,"XPos"]-min(Xboundaries_midpoint)
            # poolSeqStats[1,"YPos"]-min(Yboundaries_midpoint)
            poolSeqStats[which(theseRows2D==T),6] <- poolSeqStats[which(theseRows2D==T),4] - min(Xboundaries_midpoint)
            poolSeqStats[which(theseRows2D==T),7] <- poolSeqStats[which(theseRows2D==T),5] - min(Yboundaries_midpoint)
            
            # thisRow <- 2
            # 1. find perpendicular linear equation
            # y - y1 = m(x - x1)
            # interceptPerp <- ((4.16*(-thisX))+thisY)
            # eqnPerp <- (4.16*x) + 4788100
            #  perpendicular gradient can be found with: -(1/x)
            # find point of contact with reference line
            # y values must be equal, therefore
            # mx + b = mx + b
            # original
            # y = -(0.24*420000) + 4788100
            # perpendicular
            # y = (4.16*420000) + 4788100
            # difference <- (-(0.24*420000) + 4788100) - ((4.16*420000) + 4788100)
            # 4788100 - difference
            # Position and distance along reference line. Using basic trig.
            # 1. Perpedicular gradient = -(1/x)
            gradientPerp <- (1/transectGradient)*-1
            # 2. Finding linear equation of a line: y-y1 = m(x-x1)
            # y = mx - mx1 + y1
            interceptPerp <- ((gradientPerp*(-poolSeqStats[which(theseRows2D==T),"X_adj"]))+poolSeqStats[which(theseRows2D==T),"Y_adj"])
            # 3. Finding where two lines intersect.
            # Set them to be equal. Knowing b and b1, and setting y equal we only have to solve for x
            # i.e. mx + b = m1x + b1
            # mx - m1x = b1 - b
            # therefore, x = (b-b1)/[-(m-m1)]
            xCross <- (interceptPerp-transectIntercept)/((gradientPerp-transectGradient)*-1)
            yCross <- (transectGradient*xCross)+transectIntercept
            # yCross_start and xCross_start (first X pos of deme with data) defined with KASP sample set above
            # Lastly, determine how far is the intersection point along the reference line (compared to x=0 and y=intercept)
            #xvals <- c(-1000,xCross)
            #yvals <- c(((transectGradient*-1000)+transectIntercept),yCross)
            xvals <- c(xCross_start_set,xCross)
            yvals <- c(((transectGradient*xCross_start_set)+transectIntercept),yCross)
            
            
            # this is the line making errors at 200m demes
            
           #  cat("\n euclidean distance test 2: ",as.vector(dist(cbind(xvals,yvals),method="euclidean")))
          #   flush.console()
          #   print("poolSeqStats table")
          #   print(poolSeqStats[which(theseRows2D==T),"DistAlongTransect"])
            
            #### NOTE: turn this line below back on later
            poolSeqStats[which(theseRows2D==T),"DistAlongTransect"] <- as.vector(dist(cbind(xvals,yvals),method="euclidean"))
            
            #print (paste0("thisDistAlongTransect: ", poolSeqStats[which(theseRows2D==T),"DistAlongTransect"]))
            #curve(((transectGradient*x)+(transectIntercept)), from=min(DemeMatrix2D_combine_noZero[,"X_adj"]),
            #      to=max(DemeMatrix2D_combine_noZero[,"X_adj"]), col="darkgray",add = TRUE,type = "l",lwd=3,lty=1)
            # An alternative is to use the solve function, gives the same result
            # linearEqn <- matrix(c(transectIntercept,transectGradient,interceptPerp,gradientPerp),2,2,byrow=T)
            # XYcross <- c(-solve(cbind(linearEqn[,2],-1)) %*% linearEqn[,1])
            if (noPlots==F) { 
              if (plotPerp==T) {
                lines(x=c(poolSeqStats[which(theseRows2D==T),"X_adj"],xCross),y=c(poolSeqStats[which(theseRows2D==T),"Y_adj"],yCross),lwd=1.1,col='black')
                points(x=poolSeqStats[which(theseRows2D==T),"X_adj"],y=poolSeqStats[which(theseRows2D==T),"Y_adj"],pch=23,col='black',bg=c(rgb(255,255,255,230,maxColorValue=255)),lwd=1.5,cex=1.6)
                #points(x=xCross,y=yCross,pch=21,col='black',bg='black',cex=0.6,lwd=0.5)
              }
            }
          }
        }
      }
    }
    
    dev.off()
    
    ################################################
    # plot allele freqs by distance along transect #
    ################################################
    
    fileDemes <- paste(paste(thisLocusName,"_PoolSeq_vs_KASP_",FigName,sep=""),".pdf",sep="_")
    pdf(fileDemes,width = 14.5, height = 7,pointsize=18)
    par(mfrow=c(1,1))
    par(mar=c(3,3,3,3))
    par(oma=c(3,3,3,3))    
    SampleSize <- (DemeMatrix2D_combine_noZero[,"pCount"]+ DemeMatrix2D_combine_noZero[,"qCount"])/2
    SampleSizeMin <- SampleSize > minSample
    DemeMatrix2D_combine_noZero <- cbind(DemeMatrix2D_combine_noZero[,1:4],round(distTransect[,2],1),DemeMatrix2D_combine_noZero[,5:ncol(DemeMatrix2D_combine_noZero)])
    colnames(DemeMatrix2D_combine_noZero)[5] <- "distAlongTransect"
    # re-orientate p (magenta) & q (yellow)
    if (mean(DemeMatrix2D_combine_noZero[1:4,"p"])
        > mean(DemeMatrix2D_combine_noZero[(nrow(DemeMatrix2D_combine_noZero)-4):nrow(DemeMatrix2D_combine_noZero),"p"])) {
      DemeMatrix_temp <- DemeMatrix2D_combine_noZero
      DemeMatrix_temp[,"q"] <- DemeMatrix2D_combine_noZero[,"p"]
      DemeMatrix_temp[,"p"] <- DemeMatrix2D_combine_noZero[,"q"]
      DemeMatrix_temp[,"qCount"] <- DemeMatrix2D_combine_noZero[,"pCount"]
      DemeMatrix_temp[,"pCount"] <- DemeMatrix2D_combine_noZero[,"qCount"]
      DemeMatrix_temp[,"ppCount"] <- DemeMatrix2D_combine_noZero[,"qqCount"]
      DemeMatrix_temp[,"qqCount"] <- DemeMatrix2D_combine_noZero[,"ppCount"]
      DemeMatrix_temp[,"pp"] <- DemeMatrix2D_combine_noZero[,"qq"]
      DemeMatrix_temp[,"qq"] <- DemeMatrix2D_combine_noZero[,"pp"]
      DemeMatrix2D_combine_noZero <- DemeMatrix_temp
    }
    
    SampleSizeMin <- SampleSize >= minSample
    totalSampleSizeMinDemes <- sum(SampleSize[SampleSizeMin])
    textSec <- paste(paste("demes=",sum(SampleSizeMin),sep=""),
                     paste("n=",totalSampleSizeMinDemes,sep=""),sep=", ")
    #xMin <- min(DemeMatrix2D_combine_noZero[,"X_adj"])
    #xMax <- max(DemeMatrix2D_combine_noZero[,"X_adj"])
    #yMin <- min(DemeMatrix2D_combine_noZero[,"Y_adj"])
    #yMax <- max(DemeMatrix2D_combine_noZero[,"Y_adj"])
    
    if (noPlots==F) {
      plot(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
           DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],bg = "black",xaxs = "i", yaxs = "i",
           xlab="",ylab="",main="",xlim=c(xMin-1000,22200),
           ylim=c(-0.1,1.1),type='n',cex.main=0.9,bty="n",yaxt="n",xaxt="n")
      axis(1,at=seq(xMin,22000,1000),labels=NA,col.axis="black",las=2,cex.axis=1.1,lwd=1.6,pos=-0.01, tck = -.02)
      axis(1,at=seq(xMin,22000,2000),labels=seq(xMin,22000,2000)/1000,col.axis="black",tck=-0.3,las=1,cex.axis=1.1,lwd=0,pos=-0.02)
      axis(2,at=seq(-0.01,1.01,0.1),labels=NA,col.axis="black",las=1,cex.axis=1.1,lwd=1.6,pos=xMin,tck = -.02)
      axis(2,at=seq(0,1.0,0.2),labels=seq(0,1,0.2),col.axis="black",tck=-.03,las=1,cex.axis=1.1,lwd=0,pos=xMin)
      
      secondLevelHeader <- paste0("demeSize= ",demeSize,", gradient= ",round(transectGradient,3),", ","transect= ", transectIntercept)
      chartLabelFinal <- paste(thisLocusName,secondLevelHeader,sep="; ")
      # run details
      #title(paste(chartLabel,paste(demeSize,"m demes",sep=""),sep=", "),line=2.25,cex.main=0.9)
      title(chartLabelFinal,line=2.25,cex.main=0.9)
      
      lines(x=c(22000,22000),y=c(0,1),lwd=1.6)
      lines(x=c(xMinLine,22000),y=c(1,1),lwd=1.6)
      
      if (plotAxisLab==T) {
         mtext("Distance along transect (km)",side=1,cex=1,line=2.5)
         mtext("Allele frequency (p)",side=2,cex=1,line=2.5)
      }
        
      # Some poolSeq comparisons
      if (length(poolSeqClines)>0) {
        if (plotCentreLine ==T) {
          lines(x = c(12580,12580), y = c(-0.04,1.04),lty=5,lwd=4,col='black')
        }
        
        for (thisRow in 1:nrow(poolSeqClines)) {
          # thisRow <- 1
          poolSeqClineThisLocus <- cbind(poolSeqStats[2:7,"pool"],poolSeqStats[2:7,"DistAlongTransect"],as.vector(unlist(poolSeqClines[thisRow,9:14])))
          colnames(poolSeqClineThisLocus) <- c("pool","DistAlongTransect","p")
          poolSeqClineThisLocus <- as.data.frame(poolSeqClineThisLocus)
          poolSeqClineThisLocus[,"DistAlongTransect"] <- as.numeric(poolSeqClineThisLocus[,"DistAlongTransect"])
          poolSeqClineThisLocus[,"p"] <- as.numeric(poolSeqClineThisLocus[,"p"])
          poolSeqClineThisLocus <- poolSeqClineThisLocus[order(poolSeqClineThisLocus[,"DistAlongTransect"]),]
          
          # background loci not linked to any colour loci
          if (poolSeqClines[thisRow,"ColourGene"]==-9) {
            # lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=colourTable[colourTable[,"Region"]=="background","lineColour"])
            #back_col <- colourTable[colourTable[,"Region"]=="background","lineColour"]
            back_col <- "darkgrey"
            back_col_transparent <- adjustcolor(back_col, alpha.f = 0.8)
            lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=0.9,col=back_col_transparent)
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
            rubCol_linked_transparent <- adjustcolor(rubCol_linked, alpha.f = 0.5)
            lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=0.8,col=rubCol_linked_transparent)
          }
          if (poolSeqClines[thisRow,"ColourGene"]=="RUBIA" & poolSeqClines[thisRow,"distanceToGene"]=="inGene") {
            rubCol <- colourTable[colourTable[,"Region"]=="RUBIA","lineColour"]
            rubCol_transparent <- adjustcolor(rubCol_linked, alpha.f = 0.8)
            lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.3,col=rubCol)
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
        if (plotKASPrep==T) {
          for (thisLoc in lociNames) {
            # Cremosa KASP
            if (thisLoc == "s1187_290152") {
              # SNP KASP points
              points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                     DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="CREMOSA","pointColour"],lwd=1.5)
            }
            # Flavia KASP down
            if (thisLoc == "s316_93292" ) {
              # SNP KASP points
              points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                     DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="FLAVIA","pointColour"],lwd=1.5)
            }
            # Flavia KASP up
            if (thisLoc == "s316_257789") {
              # SNP KASP points
              points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                     DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="FLAVIA","pointColour"],lwd=1.5)
            }
            # SULF KASP
            if (thisLoc == "s91_39699") {
              # SNP KASP points
              points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                     DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="SULF","pointColour"],lwd=1.5)
            }
            # RUBIA KASP
            if (thisLoc == "s261_720757") {
              
              if (plotKASPcurve == T) {
                cline_model_sig_symm_plot <- function(XCordinate) {
                  ((pL + ((pR-pL)/(1 + exp((-4*(XCordinate-centre))/width)))))
                }
              }
                # SNP KASP points
                rubCol <- colourTable[colourTable[,"Region"]=="RUBIA","pointColour"]
                rubCol_transparent <- adjustcolor(rubCol, alpha.f = 0.3)
                points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                       DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=rubCol_transparent,lwd=1.1)
              }
              
            
            # ROS el KASP
            if (thisLoc == "ros_assembly_543443" & FigName!="Chr6") {
              # SNP KASP points
              points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                     DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=19,col='black',bg=colourTable[colourTable[,"Region"]=="ROS","pointColour"],lwd=1.5)
            }
            
            # highlight key ROS1 locus in black
            if (thisLoc == "ros_assembly_543443" & FigName=="Chr6") {
              poolSeqClineThisLocus <- cbind(poolSeqStats[2:7,"pool"],poolSeqStats[2:7,"DistAlongTransect"],as.vector(unlist(poolSeqClines[poolSeqClines[,"position"]==52889078,c(9:14)])))
              colnames(poolSeqClineThisLocus) <- c("pool","DistAlongTransect","p")
              poolSeqClineThisLocus <- as.data.frame(poolSeqClineThisLocus)
              poolSeqClineThisLocus[,"DistAlongTransect"] <- as.numeric(poolSeqClineThisLocus[,"DistAlongTransect"])
              poolSeqClineThisLocus[,"p"] <- as.numeric(poolSeqClineThisLocus[,"p"])
              poolSeqClineThisLocus <- poolSeqClineThisLocus[order(poolSeqClineThisLocus[,"DistAlongTransect"]),]
              points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                     DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="ROS","pointColour"],lwd=1.5)
              lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=2.5,col='black')
              points(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],pch=23,col='black',bg=c(rgb(255,255,255,230,maxColorValue=255)),lwd=1.5,cex=1.6)
              
            }
          }
        }
      }

    }
    
    dev.off()
    
    ##################################################
    # Allele Freq transect "Square" version of above #
    ##################################################
    
    fileDemes <- paste(paste(thisLocusName,"_PoolSeq_vs_KASP_Square",FigName,sep=""),".pdf",sep="_")

    pdf(fileDemes,width = 7, height = 7,pointsize=18)
    par(mfrow=c(1,1))
    par(mar=c(3,3,3,3))
    par(oma=c(3,3,3,3))        
    SampleSize <- (DemeMatrix2D_combine_noZero[,"pCount"]+ DemeMatrix2D_combine_noZero[,"qCount"])/2
    SampleSizeMin <- SampleSize > minSample
    DemeMatrix2D_combine_noZero <- cbind(DemeMatrix2D_combine_noZero[,1:4],round(distTransect[,2],1),DemeMatrix2D_combine_noZero[,5:ncol(DemeMatrix2D_combine_noZero)])
    colnames(DemeMatrix2D_combine_noZero)[5] <- "distAlongTransect"
    # re-orientate p (magenta) & q (yellow)
    if (mean(DemeMatrix2D_combine_noZero[1:4,"p"])
          > mean(DemeMatrix2D_combine_noZero[(nrow(DemeMatrix2D_combine_noZero)-4):nrow(DemeMatrix2D_combine_noZero),"p"])) {
        DemeMatrix_temp <- DemeMatrix2D_combine_noZero
        DemeMatrix_temp[,"q"] <- DemeMatrix2D_combine_noZero[,"p"]
        DemeMatrix_temp[,"p"] <- DemeMatrix2D_combine_noZero[,"q"]
        DemeMatrix_temp[,"qCount"] <- DemeMatrix2D_combine_noZero[,"pCount"]
        DemeMatrix_temp[,"pCount"] <- DemeMatrix2D_combine_noZero[,"qCount"]
        DemeMatrix_temp[,"ppCount"] <- DemeMatrix2D_combine_noZero[,"qqCount"]
        DemeMatrix_temp[,"qqCount"] <- DemeMatrix2D_combine_noZero[,"ppCount"]
        DemeMatrix_temp[,"pp"] <- DemeMatrix2D_combine_noZero[,"qq"]
        DemeMatrix_temp[,"qq"] <- DemeMatrix2D_combine_noZero[,"pp"]
        DemeMatrix2D_combine_noZero <- DemeMatrix_temp
      }
      
    SampleSizeMin <- SampleSize >= minSample
    totalSampleSizeMinDemes <- sum(SampleSize[SampleSizeMin])
    textSec <- paste(paste("demes=",sum(SampleSizeMin),sep=""),
                       paste("n=",totalSampleSizeMinDemes,sep=""),sep=", ")
      #xMin <- min(DemeMatrix2D_combine_noZero[,"X_adj"])
      #xMax <- max(DemeMatrix2D_combine_noZero[,"X_adj"])
      #yMin <- min(DemeMatrix2D_combine_noZero[,"Y_adj"])
      #yMax <- max(DemeMatrix2D_combine_noZero[,"Y_adj"])
    if (noPlots==F) {
      plot(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
           DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],bg = "black",xaxs = "i", yaxs = "i",
           xlab="",ylab="",main="",xlim=c(xMin-1000,22200),
           ylim=c(-0.1,1.1),type='n',cex.main=0.9,bty="n",yaxt="n",xaxt="n")
      axis(1,at=seq(xMin,22000,1000),labels=NA,col.axis="black",las=2,cex.axis=1.1,lwd=1.6,pos=-0.01, tck = -.02)
      axis(1,at=seq(xMin,22000,2000),labels=seq(xMin,22000,2000)/1000,col.axis="black",tck=-0.3,las=1,cex.axis=1.1,lwd=0,pos=-0.02)
      axis(2,at=seq(-0.01,1.01,0.1),labels=NA,col.axis="black",las=1,cex.axis=1.1,lwd=1.6,pos=xMin,tck = -.02)
      axis(2,at=seq(0,1.0,0.2),labels=seq(0,1,0.2),col.axis="black",tck=-.03,las=1,cex.axis=1.1,lwd=0,pos=xMin)
      secondLevelHeader <- paste0("demeSize= ",demeSize,", gradient= ",round(transectGradient,3),", ","transect= ", transectIntercept)
      chartLabelFinal <- paste(thisLocusName,secondLevelHeader,sep="; ")
      # run details
      #title(paste(chartLabel,paste(demeSize,"m demes",sep=""),sep=", "),line=2.25,cex.main=0.9)
      title(chartLabelFinal,line=0.1,cex.main=0.5)
      
      #title(paste(chartLabel,textSec,sep="\n"),line=0.1,cex.main=0.9)
      lines(x=c(22000,22000),y=c(0,1),lwd=1.6)
      lines(x=c(xMinLine,22000),y=c(1,1),lwd=1.6)
      
      if (plotAxisLab==T) {
        mtext("Distance along transect (km)",side=1,cex=1,line=2.5)
        mtext("Allele frequency (p)",side=2,cex=1,line=2.5)
      }
      
      # Some poolSeq comparisons
      clineLWD <- 0.7
      if (poolSeqPlot!=F) {
        if (length(poolSeqClines)>0) {
          for (thisRow in 1:nrow(poolSeqClines)) {
            # thisRow <- 1
            poolSeqClineThisLocus <- cbind(poolSeqStats[2:7,"pool"],poolSeqStats[2:7,"DistAlongTransect"],as.vector(unlist(poolSeqClines[thisRow,9:14])))
            colnames(poolSeqClineThisLocus) <- c("pool","DistAlongTransect","p")
            poolSeqClineThisLocus <- as.data.frame(poolSeqClineThisLocus)
            poolSeqClineThisLocus[,"DistAlongTransect"] <- as.numeric(poolSeqClineThisLocus[,"DistAlongTransect"])
            poolSeqClineThisLocus[,"p"] <- as.numeric(poolSeqClineThisLocus[,"p"])
            poolSeqClineThisLocus <- poolSeqClineThisLocus[order(poolSeqClineThisLocus[,"DistAlongTransect"]),]
            
            # background loci not linked to any colour loci
            if (poolSeqClines[thisRow,"ColourGene"]==-9) {
              lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=colourTable[colourTable[,"Region"]=="background","lineColour"])
            }
            # ROS EL
            if ((poolSeqClines[thisRow,"ColourGene"]=="ROSEA" | poolSeqClines[thisRow,"ColourGene"]=="ELUTA") & poolSeqClines[thisRow,"distanceToGene"]!="inGene") {
              lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=colourTable[colourTable[,"Region"]=="ROSelLinked","lineColour"])
            }
            if (poolSeqClines[thisRow,"ColourGene"]=="ROS1" | poolSeqClines[thisRow,"ColourGene"]=="ROS2" | poolSeqClines[thisRow,"ColourGene"]=="ROS3") {
              lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=colourTable[colourTable[,"Region"]=="ROS","lineColour"])
            }
            if (poolSeqClines[thisRow,"ColourGene"]=="ELUTA") {
              lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=colourTable[colourTable[,"Region"]=="ELUTA","lineColour"])
            }
            # RUBIA
            if (poolSeqClines[thisRow,"ColourGene"]=="RUBIA" & poolSeqClines[thisRow,"distanceToGene"]!="inGene") {
              lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=colourTable[colourTable[,"Region"]=="RUBIAlinked","lineColour"])
            }
            if (poolSeqClines[thisRow,"ColourGene"]=="RUBIA" & poolSeqClines[thisRow,"distanceToGene"]=="inGene") {
              lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=1.1,col=colourTable[colourTable[,"Region"]=="RUBIA","lineColour"])
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
          if (plotKASPrep==T) {
            for (thisLoc in lociNames) {
              # Cremosa KASP
              if (thisLoc == "s1187_290152") {
                # SNP KASP points
                points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                       DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="CREMOSA","pointColour"],lwd=1.5)
              }
              # Flavia KASP down
              if (thisLoc == "s316_93292" ) {
                # SNP KASP points
                points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                       DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="FLAVIA","pointColour"],lwd=1.5)
              }
              # Flavia KASP up
              if (thisLoc == "s316_257789") {
                # SNP KASP points
                points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                       DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="FLAVIA","pointColour"],lwd=1.5)
              }
              # SULF KASP
              if (thisLoc == "s91_39699") {
                # SNP KASP points
                points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                       DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="SULF","pointColour"],lwd=1.5)
              }
              # RUBIA KASP
              if (thisLoc == "s261_720757") {
                # SNP KASP points
                points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                       DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="RUBIA","pointColour"],lwd=1.5)
              }
              # ROS el KASP
              if (thisLoc == "ros_assembly_543443" & FigName!="Chr6") {
                # SNP KASP points
                points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                       DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="ROS","pointColour"],lwd=1.5)
              }
              
              # highlight key ROS1 locus in black
              if (thisLoc == "ros_assembly_543443" & FigName=="Chr6") {
                poolSeqClineThisLocus <- cbind(poolSeqStats[2:7,"pool"],poolSeqStats[2:7,"DistAlongTransect"],as.vector(unlist(poolSeqClines[poolSeqClines[,"position"]==52889078,c(9:14)])))
                colnames(poolSeqClineThisLocus) <- c("pool","DistAlongTransect","p")
                poolSeqClineThisLocus <- as.data.frame(poolSeqClineThisLocus)
                poolSeqClineThisLocus[,"DistAlongTransect"] <- as.numeric(poolSeqClineThisLocus[,"DistAlongTransect"])
                poolSeqClineThisLocus[,"p"] <- as.numeric(poolSeqClineThisLocus[,"p"])
                poolSeqClineThisLocus <- poolSeqClineThisLocus[order(poolSeqClineThisLocus[,"DistAlongTransect"]),]
                points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],
                       DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],pch=21,col='black',bg=colourTable[colourTable[,"Region"]=="ROS","pointColour"],lwd=1.5)
                lines(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],lty=1,lwd=2.5,col='black')
                points(x=poolSeqClineThisLocus[,"DistAlongTransect"],y=poolSeqClineThisLocus[,"p"],pch=23,col='black',bg=c(rgb(255,255,255,230,maxColorValue=255)),lwd=1.5,cex=1.6)
                
              }
            }
          }
        }
      }
      if (plotCentreLine ==T) {
        lines(x = c(12580,12580), y = c(-0.04,1.04),lty=5,lwd=4,col='red')
      }
    }
    dev.off()
    
    # Send data to list
    totalSampleSize2D <- (sum(DemeMatrix2D_combine_noZero[,"pCount"])+sum(DemeMatrix2D_combine_noZero[,"qCount"]))/2
    # send all 2D data to list
    DemeData[[thisLocusName]][["DemeMatrix2D_pCount"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_qCount"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_Ngenes"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_combine_noZero"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_ppCount"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_qqCount"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_pqCount"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_pp"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_qq"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_pq"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_Individuals"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_DistanceRefLine"]] <- list()
    DemeData[[thisLocusName]][["poolSeqStats"]] <- list()
    
    DemeData[[thisLocusName]][["DemeMatrix2D_pCount"]] <- DemeMatrix2D_pCount
    DemeData[[thisLocusName]][["DemeMatrix2D_qCount"]] <- DemeMatrix2D_qCount
    DemeData[[thisLocusName]][["DemeMatrix2D_Ngenes"]] <- DemeMatrix2D_Ngenes
    DemeData[[thisLocusName]][["DemeMatrix2D_ppCount"]] <- DemeMatrix2D_ppCount
    DemeData[[thisLocusName]][["DemeMatrix2D_qqCount"]] <- DemeMatrix2D_qqCount
    DemeData[[thisLocusName]][["DemeMatrix2D_pqCount"]] <- DemeMatrix2D_pqCount
    DemeData[[thisLocusName]][["DemeMatrix2D_pp"]] <- DemeMatrix2D_pp
    DemeData[[thisLocusName]][["DemeMatrix2D_qq"]] <- DemeMatrix2D_qq
    DemeData[[thisLocusName]][["DemeMatrix2D_pq"]] <- DemeMatrix2D_pq
    DemeData[[thisLocusName]][["DemeMatrix2D_combine_noZero"]] <- DemeMatrix2D_combine_noZero
    DemeData[[thisLocusName]][["DemeMatrix2D_Individuals"]] <- DemeIndividuals
    DemeData[[thisLocusName]][["DemeMatrix2D_DistanceRefLine"]] <- distTransect
    DemeData[[thisLocusName]][["poolSeqStats"]] <- poolSeqStats
    
    # clear
    rm(DemeMatrix2D_p,DemeMatrix2D_q,DemeMatrix2D_pq,DemeMatrix2D_Ngenes,DemeMatrix2D_pCount,DemeMatrix2D_combine_noZero,
       DemeMatrix2D_qCount,DemeMatrix2D_ppCount,DemeMatrix2D_qqCount,DemeMatrix2D_pqCount,DemeMatrix2D_pp,DemeMatrix2D_qq,
       DemeMatrix2D_inds,DemeMatrix2D_genos,DemeMatrix2D_combine)
  } # end locus loop
  return(DemeData)
  # 1/4/2024
}
