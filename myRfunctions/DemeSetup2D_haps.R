# DemeSetup2D_haps
# Author: David. L. Field
# 08/05/2024

# This is a function for separating individuals into demes in two-dimensions and calculating various pop gen statistics 
# Designed for calculating haplotypes frequencies (parental and recombinant) at ROS, SULF and FLA loci. 
# For genotypes at standard SNP loci see the other script (DemeSetup2D_genos).
# NOTE: does not work for standard SNP loci, only for specific coding for haplotypes in the final data set
# Haplotypes were determined in the pipeline script 'setup4_EcolGenoFinal.Rmd'. Columns (e.g. Fla_2loc) 
# Produces a list for each deme, containing a datatable with estimates various values such as haplotype frequencies, 
# in demes in two-dimensions. 
# HZdata = the data,
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

DemeSetup2D_haplos <- function(HZdata,lociNames,demeSize,minSample,removePhenoMiss,transectGradient,transectIntercept,
                        adjustXY,pointsCol,gridAdd,plotAllele1D,plotDem,addPie,SNPlociListUpdated,verbose=F,
                        plotAxisLab,plotPerp,poolSeqPlot,plotType,plotMap,noPlots=T) {
  # HZdata <- HZ_Planoles
  # nrow(HZdata)
  # HZdata[,"ros_assembly_543443"]
  # lociNames <- locusSet_FLA
  # lociNames <- c("Fla_2loc")
  # lociNames <- c("s1187_290152","s316_93292", "s316_257789","ros_assembly_543443")
  # plotAxisLab <- T ; demeSize <- 200; minSample <- 10; removePhenoMiss <- FALSE; gridAdd <- TRUE; poolSeqPlot=T
  # plotAllele1D <- T; location <- "Planoles"; plotPerp = T; plotType ="pdf"; verbose=T; plotMap=T; addPie=T
  # in Easting
  # transectGradient <- -0.24; transectIntercept <- 4788100 # Planoles
  # transectGradient <- 3.125; transectIntercept <- 3412000 # Cadi
  # in metres
  # transectGradient <- -0.20; transectIntercept <- 4700 # Planoles
  # transectGradient <- -0.24; transectIntercept <- 6000 # Planoles
  
  # transectGradient <- 3.125; transectIntercept <- 3412000 # Cadi
  # adjustXY <- FALSE; pointsCol <- FALSE; pieRadius <- 60; plotDem <- c(2,2); addPie <- T; verbose=T; poolSeqPlot=T
  # 
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
                       paste(paste("LG",linkageGroup,sep=" "),
                             paste("cM",cM,sep=" "),sep=", "),sep=", ")
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
  
  # New approach 2016! Using midpoints of square grid
  Xboundaries <- sort(seq(minX,maxX,demeSize))
  Yboundaries <- sort(seq(minY,maxY,demeSize))
  Xboundaries_midpoint <- Xboundaries+(demeSize/2)
  Yboundaries_midpoint <- Yboundaries+(demeSize/2)
 
  # Deme data stored in a list for each locus
  DemeData <- list()
  if (verbose==T) {
    cat(c(", Locus.."))
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
    goodHaps <- c("0,0","0,1","0,2","2,0","1,1","1,2","2,2")
    # note dropped 1,0 and 2,1 because they are extremely rare (perhaps an error?)
    # summaryfullData_thisLocus_backup <- table(backup)
    summaryfullData_thisLocus <- summaryfullData_thisLocus[names(summaryfullData_thisLocus)%in%goodHaps]
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
    DemeMatrix2D_r <- DemeMatrix2D
    # allele counts
    DemeMatrix2D_Ngenes <- DemeMatrix2D
    DemeMatrix2D_pCount <- DemeMatrix2D
    DemeMatrix2D_qCount <- DemeMatrix2D
    DemeMatrix2D_rCount <- DemeMatrix2D
    # haplotype counts
    DemeMatrix2D_FLA_FLA_Count <- DemeMatrix2D
    DemeMatrix2D_FLAr_FLA_Count <- DemeMatrix2D
    DemeMatrix2D_FLAr_FLAr_Count <- DemeMatrix2D
    DemeMatrix2D_FLA_fla_Count <- DemeMatrix2D
    DemeMatrix2D_FLAr_fla_Count <- DemeMatrix2D
    DemeMatrix2D_fla_fla_Count <- DemeMatrix2D
    # haplotype frequencies
    DemeMatrix2D_FLA_FLA <- DemeMatrix2D
    DemeMatrix2D_FLAr_FLA <- DemeMatrix2D
    DemeMatrix2D_FLAr_FLAr <- DemeMatrix2D
    DemeMatrix2D_FLA_fla <- DemeMatrix2D
    DemeMatrix2D_FLAr_fla <- DemeMatrix2D
    DemeMatrix2D_fla_fla <- DemeMatrix2D
    # genotype frequencies
    #DemeMatrix2D_pp <- DemeMatrix2D
    #DemeMatrix2D_qq <- DemeMatrix2D
    #DemeMatrix2D_pq <- DemeMatrix2D
    
    # Individuals 
    DemeMatrix2D_inds <- DemeMatrix2D
    # Haplotypes
    DemeMatrix2D_haplos <- DemeMatrix2D
    
    # Keep list of individuals
    DemeIndividuals <- matrix(0,ncol(DemeMatrix2D_p)*nrow(DemeMatrix2D_p),7)
    colnames(DemeIndividuals) <- c("Xmin","Xmax","Ymin","Ymax","numInds","Inds","Genos")
    demeCounter <- 0
    if (plotMap==T) {
      par(mfrow=c(1,1))
      # cores and flanks - 421000,426000,4685500,4687800
      plotHZbasic(HZdata,409000,437000,4684500,4692000,TRUE,FALSE,TitleStart=paste("Planoles ",thisLocusName,sep=""),plotAge=F)
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
          
          if (all(jLoc_kDeme%in%names(summaryfullData_thisLocus)==F)) {
            next
          }
          else {
            
            # cut back to good haplotypes
            TotalCount <- length(jLoc_kDeme)
            jLoc_kDeme <- jLoc_kDeme[jLoc_kDeme%in%goodHaps]
              # (sum(jLoc_kDeme==0) + sum(jLoc_kDeme==1) + sum(jLoc_kDeme==2) + sum(jLoc_kDeme=="-9"))
            Ncount <- length(jLoc_kDeme)
            Ngenes <- 2*length(jLoc_kDeme)
            # make r the recombinant haplotype
            pCount <- (2*sum(jLoc_kDeme=="0,0")) + sum(jLoc_kDeme=="0,1") + sum(jLoc_kDeme=="1,1")
            qCount <- (2*sum(jLoc_kDeme=="2,2")) + sum(jLoc_kDeme=="1,2") + sum(jLoc_kDeme=="1,1") 
            rCount <- (2*sum(jLoc_kDeme=="0,2")) + sum(jLoc_kDeme=="0,1") + sum(jLoc_kDeme=="1,2")
            missCount <- TotalCount-Ncount
            TotalGenotyped <- TotalCount
            # Haplotypes
            FLA_FLA_Count <- sum(jLoc_kDeme=="0,0")
            FLAr_FLA_Count <- sum(jLoc_kDeme=="0,1")
            FLAr_FLAr_Count <- sum(jLoc_kDeme=="0,2")
            FLA_fla_Count <- sum(jLoc_kDeme=="1,1")
            FLAr_fla_Count <- sum(jLoc_kDeme=="1,2")
            fla_fla_Count <- sum(jLoc_kDeme=="2,2")
            # frequencies
            p_freq <- pCount/Ngenes
            q_freq <- qCount/Ngenes
            r_freq <- rCount/Ngenes
            FLA_FLA_freq <- FLA_FLA_Count/Ncount
            FLAr_FLA_freq <- FLAr_FLA_Count/Ncount
            FLAr_FLAr_freq <- FLAr_FLAr_Count/Ncount
            FLA_fla_freq <- FLA_fla_Count/Ncount
            FLAr_fla_freq <- FLAr_fla_Count/Ncount
            fla_fla_freq <- fla_fla_Count/Ncount
            
            # send to 2D matrices
            DemeMatrix2D_p[thisY,thisX] <- p_freq
            DemeMatrix2D_q[thisY,thisX] <- q_freq
            DemeMatrix2D_r[thisY,thisX] <- r_freq
            # allele counts
            DemeMatrix2D_Ngenes[thisY,thisX] <- Ngenes
            DemeMatrix2D_pCount[thisY,thisX] <- pCount
            DemeMatrix2D_qCount[thisY,thisX] <- qCount
            DemeMatrix2D_rCount[thisY,thisX] <- rCount
            # haplotype counts
            DemeMatrix2D_FLA_FLA_Count[thisY,thisX] <- FLA_FLA_Count
            DemeMatrix2D_FLAr_FLA_Count[thisY,thisX] <- FLAr_FLA_Count
            DemeMatrix2D_FLAr_FLAr_Count[thisY,thisX] <- FLAr_FLAr_Count
            DemeMatrix2D_FLA_fla_Count[thisY,thisX] <- FLA_fla_Count
            DemeMatrix2D_FLAr_fla_Count[thisY,thisX] <- FLAr_fla_Count
            DemeMatrix2D_fla_fla_Count[thisY,thisX] <- fla_fla_Count
            # haplotype frequencies
            DemeMatrix2D_FLA_FLA[thisY,thisX] <- FLA_FLA_freq
            DemeMatrix2D_FLAr_FLA[thisY,thisX] <- FLAr_FLA_freq
            DemeMatrix2D_FLAr_FLAr[thisY,thisX] <- FLAr_FLAr_freq
            DemeMatrix2D_FLA_fla[thisY,thisX] <- FLA_fla_freq
            DemeMatrix2D_FLAr_fla[thisY,thisX] <- FLAr_fla_freq
            DemeMatrix2D_fla_fla[thisY,thisX] <- fla_fla_freq
            
            # record individuals
            DemeMatrix2D_inds[thisY,thisX] <- paste(as.vector(inds_kDeme),collapse="/")
            # Haplotypes
            DemeMatrix2D_haplos[thisY,thisX] <- paste(jLoc_kDeme,collapse=",")
            # DemeIndividuals[demeCounter,1:2] <- theseX
            # DemeIndividuals[demeCounter,3:4] <- theseY
            # DemeIndividuals[demeCounter,"Inds"] <- paste(as.vector(HZdata[theseRows2D,"PlantID_final"]),collapse="/")
            # DemeIndividuals[demeCounter,"numInds"] <- length(as.vector(HZdata[theseRows2D,"PlantID_final"]))
            # DemeIndividuals[demeCounter,"Genos"] <- paste(jLoc_kDeme,collapse=",")
            if (plotMap==T) {
              myCols <- rainbow(20) # 
              points(cbind(HZdata[theseRows2D,"Easting"],HZdata[theseRows2D,"Northing"]), pch = 21, cex=0.8,
                     col="black",bg=sample(myCols,1))
            }
          }
        }
      }
    }
    # write.csv(DemeIndividuals,"DemeIndividuals.csv")
    # individuals
    
    DemeIndividuals_melted <- melt(t(DemeIndividuals))
    
    # allele freq tables
    DemeMatrix2D_p_melted <- melt(t(DemeMatrix2D_p))
    colnames(DemeMatrix2D_p_melted) <- c("X","Y","p")
    DemeMatrix2D_q_melted <- melt(t(DemeMatrix2D_q))
    colnames(DemeMatrix2D_q_melted) <- c("X","Y","q")
    DemeMatrix2D_r_melted <- melt(t(DemeMatrix2D_r))
    colnames(DemeMatrix2D_r_melted) <- c("X","Y","r")
    # allele count tables
    DemeMatrix2D_pCount_melted <- melt(t(DemeMatrix2D_pCount))
    colnames(DemeMatrix2D_pCount_melted) <- c("X","Y","pCount")
    DemeMatrix2D_qCount_melted <- melt(t(DemeMatrix2D_qCount))
    colnames(DemeMatrix2D_qCount_melted) <- c("X","Y","qCount")
    DemeMatrix2D_rCount_melted <- melt(t(DemeMatrix2D_rCount))
    colnames(DemeMatrix2D_rCount_melted) <- c("X","Y","rCount")
    # haplotype count tables
    DemeMatrix2D_FLA_FLA_Count_melted <- melt(t(DemeMatrix2D_FLA_FLA_Count))
    colnames(DemeMatrix2D_FLA_FLA_Count_melted) <- c("X","Y","FLA_FLA_Count")
    DemeMatrix2D_FLAr_FLA_Count_melted <- melt(t(DemeMatrix2D_FLAr_FLA_Count))
    colnames(DemeMatrix2D_FLAr_FLA_Count_melted) <- c("X","Y","FLAr_FLA_Count")
    DemeMatrix2D_FLAr_FLAr_Count_melted <- melt(t(DemeMatrix2D_FLAr_FLAr_Count))
    colnames(DemeMatrix2D_FLAr_FLAr_Count_melted) <- c("X","Y","FLAr_FLAr_Count")
    DemeMatrix2D_FLA_fla_Count_melted <- melt(t(DemeMatrix2D_FLA_fla_Count))
    colnames(DemeMatrix2D_FLA_fla_Count_melted) <- c("X","Y","FLA_fla_Count")
    DemeMatrix2D_FLAr_fla_Count_melted <- melt(t(DemeMatrix2D_FLAr_fla_Count))
    colnames(DemeMatrix2D_FLAr_fla_Count_melted) <- c("X","Y","FLAr_fla_Count")
    DemeMatrix2D_fla_fla_Count_melted <- melt(t(DemeMatrix2D_fla_fla_Count))
    colnames(DemeMatrix2D_fla_fla_Count_melted) <- c("X","Y","fla_fla_Count")
    # haplotype frequency tables
    DemeMatrix2D_FLA_FLA_melted <- melt(t(DemeMatrix2D_FLA_FLA))
    colnames(DemeMatrix2D_FLA_FLA_melted) <- c("X","Y","FLA_FLA_freq")
    DemeMatrix2D_FLAr_FLA_melted <- melt(t(DemeMatrix2D_FLAr_FLA))
    colnames(DemeMatrix2D_FLAr_FLA_melted) <- c("X","Y","FLAr_FLA_freq")
    DemeMatrix2D_FLAr_FLAr_melted <- melt(t(DemeMatrix2D_FLAr_FLAr))
    colnames(DemeMatrix2D_FLAr_FLAr_melted) <- c("X","Y","FLAr_FLAr_freq")
    DemeMatrix2D_FLA_fla_melted <- melt(t(DemeMatrix2D_FLA_fla))
    colnames(DemeMatrix2D_FLA_fla_melted) <- c("X","Y","FLA_fla_freq")
    DemeMatrix2D_FLAr_fla_melted <- melt(t(DemeMatrix2D_FLAr_fla))
    colnames(DemeMatrix2D_FLAr_fla_melted) <- c("X","Y","FLAr_fla_freq")
    DemeMatrix2D_fla_fla_melted <- melt(t(DemeMatrix2D_fla_fla))
    colnames(DemeMatrix2D_fla_fla_melted) <- c("X","Y","fla_fla_freq")
   
    # Inds
    DemeMatrix2D_inds_melted <- melt(t(DemeMatrix2D_inds))
    colnames(DemeMatrix2D_inds_melted) <- c("X","Y","inds")
    # Genotypes
    DemeMatrix2D_haplos_melted <- melt(t(DemeMatrix2D_haplos))
    colnames(DemeMatrix2D_haplos_melted) <- c("X","Y","haplos")
    # nrow(DemeMatrix2D_pq_melted); nrow(DemeMatrix2D_inds_melted); nrow(DemeMatrix2D_genos_melted)
    # convert Easting to start from 0 (+1000). This starting at 1000 is to avoid tricky intercept problems for far West demes
    X_adj <- (DemeMatrix2D_pCount_melted[,"X"])-min(Xboundaries_midpoint)
    Y_adj <- (DemeMatrix2D_pCount_melted[,"Y"])-min(Yboundaries_midpoint)
    
    # combine
    DemeMatrix2D_combine <- cbind(DemeMatrix2D_pCount_melted[,1:2],X_adj,Y_adj,DemeMatrix2D_pCount_melted$pCount,
                                  DemeMatrix2D_qCount_melted$qCount,DemeMatrix2D_rCount_melted$rCount,
                                  DemeMatrix2D_p_melted$p,DemeMatrix2D_q_melted$q,DemeMatrix2D_r_melted$r,
                                  DemeMatrix2D_FLA_FLA_Count_melted$FLA_FLA_Count,
                                  DemeMatrix2D_FLAr_FLA_Count_melted$FLAr_FLA_Count,
                                  DemeMatrix2D_FLAr_FLAr_Count_melted$FLAr_FLAr_Count,
                                  DemeMatrix2D_FLA_fla_Count_melted$FLA_fla_Count,
                                  DemeMatrix2D_FLAr_fla_Count_melted$FLAr_fla_Count,
                                  DemeMatrix2D_fla_fla_Count_melted$fla_fla_Count,
                                  DemeMatrix2D_FLA_FLA_melted$FLA_FLA_freq,
                                  DemeMatrix2D_FLAr_FLA_melted$FLAr_FLA_freq,
                                  DemeMatrix2D_FLAr_FLAr_melted$FLAr_FLAr_freq,
                                  DemeMatrix2D_FLA_fla_melted$FLA_fla_freq,
                                  DemeMatrix2D_FLAr_fla_melted$FLAr_fla_freq,
                                  DemeMatrix2D_fla_fla_melted$fla_fla_freq,
                                  as.vector(DemeMatrix2D_inds_melted$inds),
                                  as.vector(DemeMatrix2D_haplos_melted$haplos))
    colnames(DemeMatrix2D_combine) <- c("X","Y","X_adj","Y_adj","pCount","qCount","rCount","p","q","r",
                                        "FLA_FLA_Count","FLAr_FLA_Count","FLAr_FLAr_Count","FLA_fla_Count","FLAr_fla_Count",
                                        "fla_fla_Count",
                                        "FLA_FLA_freq","FLAr_FLA_freq","FLAr_FLAr_freq","FLA_fla_freq","FLAr_fla_freq",
                                        "fla_fla_freq",
                                        "inds","haplos")
    DemeMatrix2D_combine_noZero <- DemeMatrix2D_combine[DemeMatrix2D_combine[,"pCount"]!=-9,]
    #head(DemeMatrix2D_pqCount)
    # DemeMatrix2D_combine[DemeMatrix2D_combine[,"pp"]!=-9,c(1:10)]
    # just in case the demes with no data vary among loci - we find the minX now.
    xMin <- min(DemeMatrix2D_combine[,"X_adj"])
    xMax <- max(DemeMatrix2D_combine[,"X_adj"])
    yMin <- min(DemeMatrix2D_combine[,"Y_adj"])
    yMax <- max(DemeMatrix2D_combine[,"Y_adj"])
    ##########################
    # Pie plot with transect #
    ##########################
    if (noPlots==F) {
      if (length(lociNames)>=4) {
        # now set up multipanel
        par(mfrow=c(2,1))
        par(mar=c(3,3,3,1))
        par(oma=c(3,3,3,1))
      }
      if (length(lociNames)<4) {
        # now set up multipanel
        par(mfrow=c(2,1))
        par(mar=c(3,3,3,1))
        par(oma=c(3,3,3,1))
      } 
      plot(DemeMatrix2D_combine_noZero[,"X_adj"],DemeMatrix2D_combine_noZero[,"Y_adj"],type='n',
           main=paste(thisLocusName,paste(demeSize,"m demes",sep=""),sep=", "),cex.main=1,bty="n",yaxt="n",xaxt="n",
           xlab = "",ylab = "",col = "black", bg = "black",xaxs = "i", yaxs = "i",xlim=c(xMin-2000,24300),
           ylim=c(yMin-4000,10050))
      axis(1,at=seq(xMin-2000,xMax,1000),labels=NA,col.axis="black",las=2,cex.axis=1.1,lwd=1.6,pos=yMin-4000, tck = -.02)
      axis(1,at=seq(-2000,xMax,2000),labels=seq(xMin-2000,xMax,2000)/1000,col.axis="black",tck=-0.3,las=1,cex.axis=1.1,lwd=0,pos=yMin-4000.5)
      axis(2,at=seq(yMin-4000,yMax+3800,1000),labels=NA,col.axis="black",las=1,cex.axis=1.1,lwd=1.6,pos=xMin-2000,tck = -.02)
      axis(2,at=seq(yMin-4000,yMax+3800,2000),labels=seq(yMin-4000,yMax+3800,2000)/1000,col.axis="black",tck=-.03,las=1,cex.axis=1.1,lwd=0,pos=xMin-2000.5)
      lines(x=c(24000,24000),y=c(yMin-4000,yMax+3000),lwd=1.6)
      lines(x=c(xMin-2000,24000),y=c(yMax+3000,yMax+3000),lwd=1.6)
      
      curve(((transectGradient*x)+(transectIntercept)), from=c(min(DemeMatrix2D_combine_noZero[,"X_adj"])-1000,
                                                               to=max(DemeMatrix2D_combine_noZero[,"X_adj"])+1000), col="darkgray",add = TRUE,type = "l",lwd=2,lty=1)
      mtext("East (km)",side=1,cex=1,line=2.5)
      mtext("North (km)",side=2,cex=1,line=2.5)
      if (addPie==T) {
        # plotting pie charts
        # minSample <- 5
        totalN <- (DemeMatrix2D_combine_noZero[,"pCount"]+DemeMatrix2D_combine_noZero[,"qCount"]+DemeMatrix2D_combine_noZero[,"rCount"])/2
        keep <- which(totalN>=minSample)
        DemeMatrix2D_combine_noZero <- DemeMatrix2D_combine_noZero[keep,]
        allthisP <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[,"pCount"]))
        allthisQ <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[,"qCount"]))
        allthisR <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[,"rCount"]))
        allSum <- allthisP+allthisQ+allthisR
        allFreq <- allthisP/allSum
        if (mean(allFreq[1:4])>0.6) {
          allthisP_temp <- allthisQ
          allthisQ_temp <- allthisP
          allthisQ <- allthisP_temp
          allthisP <- allthisQ_temp
        }
        for (thisDeme2D in 1:nrow(DemeMatrix2D_combine_noZero)) {
          # thisDeme2D <- 1
          thisXraw <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[thisDeme2D,"X_adj"]))
          thisYraw <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[thisDeme2D,"Y_adj"]))
          thisP <- as.numeric(allthisP[thisDeme2D])
          thisQ <- as.numeric(allthisQ[thisDeme2D])
          thisR <- as.numeric(allthisR[thisDeme2D])
          
          # pie plot size has to be tweaked depending on deme size
          # small pie plots, original was radius 69 for deme size 200
          #     if (sum(thisP,thisQ)>minNum) {
          #       add.pie(z=c(thisP,thisQ), x=thisXmid, y=thisYmid, radius=69,col=c((rgb(255,0,0,255,maxColorValue=255)),(rgb(255,255,0,255,maxColorValue=255))),labels="")
          #     }
          # larger pie plots - 120 works okay for deme size 400
          totalN <- thisP + thisQ + thisR
          #if ((totalN < minNum)==F) {
          add.pie(z=c(thisP,thisR,thisQ), x=thisXraw, y=thisYraw, radius=85,col=c((rgb(255,0,0,255,maxColorValue=255)),"blue",(rgb(255,255,0,255,maxColorValue=255))),labels="")
          #}
        }
      }
      if (addPie==F) {
        points(DemeMatrix2D_combine_noZero[,"X_adj"],DemeMatrix2D_combine_noZero[,"Y_adj"])
      }
      
      #######################################
      # Calculating distance along transect #
      #######################################
      # now find the intersection to transect and how far along the transect each deme sits
      distTransect <- matrix(0,nrow(DemeMatrix2D_combine_noZero),2)
      colnames(distTransect) <- c("DistanceRefaxis","DistanceRefaxis_adj")
      # calculate distance along transect
      xCross_start_set <- NULL
      for (thisRow in 1:nrow(distTransect)) {
        # thisRow <- 2
        # 1. find perpendicular linear equation
        # y - y1 = m(x - x1)
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
        if (plotPerp==T) {
          lines(x=c(thisX,xCross),y=c(thisY,yCross),lwd=0.8,col='blue')
          #points(x=xCross,y=yCross,pch=21,col='red',cex=0.6,lwd=0.5)
        }
      }
      # original in column 1
      distTransect[,1] <- distTransect[,2] 
      # rescale to min val (may need to keep an eye on this and make it further negative , <5000. Depends on where the perp hits)
      #distTransect[,2] <- distTransect[,2] - min(distTransect[,2])
      
      #######################
      # Pool Seq comparison #
      #######################
      # only do this once
      thisLocus <- 1
      if (poolSeqPlot==T & thisLocus==1) {
        # setup matrix for poolSeq comparison
        # 1. YP4 ??? 52 random plants from yellow parapatric
        # 2. yelPool20 ??? 20 yellow plants between yellow flank and parapatric region (individuals J1141-J1160)
        # 3. YP2 ??? 50 random plants from yellow lower flank
        # 4. YP1 ??? 50 random plants from yellow diagonal flank
        # 5. MP2 - 50 random plants from magenta lower flank
        # 6. MP4 - 50 random plants from magenta lower flank
        # 7. mgtPool20 ??? 20 magenta plants between magenta flank and parapatric region (individuals J1121-J1140)
        # 8. MP11 ??? 50 random plants from magenta parapatric
        
        poolLoc <- matrix(0,7,8)
        poolLoc <- as.data.frame(poolLoc)
        colnames(poolLoc) <- c("pool","Easting","Northing","XPos_deme","YPos_deme","X_adj","Y_adj","DistAlongTransect")
        poolLoc[1:7,1] <- c("DemeKASP_Start","YP4","YP2","YP1","MP2","MP4","MP11")
        poolLoc[1,4:5] <- DemeMatrix2D_combine_noZero[1,1:2]
        poolLoc[1,6:7] <- DemeMatrix2D_combine_noZero[1,3:4]
        poolLoc[2,2:3] <- c(411604.6046997,4690386.0014007)
        poolLoc[3,2:3] <- c(421903.0320033,4686484.3909220)
        poolLoc[4,2:3] <- c(422119.7797325,4686390.2520109)
        poolLoc[5,2:3] <- c(424462.3465657,4686112.3149515)
        poolLoc[6,2:3] <- c(425146.6755758,4686021.5773514)
        poolLoc[7,2:3] <- c(431621.6848923,4686847.1448457)
        poolLoc[,2] <- as.numeric(as.vector(poolLoc[,2]))
        poolLoc[,3] <- as.numeric(as.vector(poolLoc[,3]))
        poolLoc[,4] <- as.numeric(as.vector(poolLoc[,4]))
        poolLoc[,5] <- as.numeric(as.vector(poolLoc[,5]))
        poolLoc[,6] <- as.numeric(as.vector(poolLoc[,6]))
        poolLoc[,7] <- as.numeric(as.vector(poolLoc[,7]))
        Deme2D_pools <- DemeMatrix2D
        Deme2D_pools[,] <- -9
        
        # find pool locations and distance along transect
        # 1. YP4 ??? 52 random plants from yellow parapatric
        # 2. yelPool20 ??? 20 yellow plants between yellow flank and parapatric region (individuals J1141-J1160)
        # 3. YP2 ??? 50 random plants from yellow lower flank
        # 4. YP1 ??? 50 random plants from yellow diagonal flank
        # 5. MP2 - 50 random plants from magenta lower flank
        # 6. MP4 - 50 random plants from magenta lower flank
        # 7. mgtPool20 ??? 20 magenta plants between magenta flank and parapatric region (individuals J1121-J1140)
        # 8. MP11 ??? 50 random plants from magenta parapatric
        # find locations for poolSeq data (no genotype calculations)
        demeCounter <- 0
        theseRows2D <- NULL
        for (thisXcol in 1:ncol(Deme2D_pools)) {
          for (thisYrow in 1:nrow(Deme2D_pools)) {
            demeCounter <- demeCounter + 1
            # thisXcol <- 7
            # thisYrow <- 5
            thisDeme_mid_X <- as.numeric(colnames(Deme2D_pools)[thisXcol])
            thisDeme_mid_Y <- as.numeric(rownames(Deme2D_pools)[thisYrow])
            
            theseX <- c(as.numeric(colnames(Deme2D_pools)[thisXcol])-(demeSize/2),as.numeric(colnames(Deme2D_pools)[thisXcol])+(demeSize/2))
            theseY <- c(as.numeric(rownames(Deme2D_pools)[thisYrow])-(demeSize/2),as.numeric(rownames(Deme2D_pools)[thisYrow])+(demeSize/2))
            theseRows2D <- as.numeric(poolLoc[,"Easting"]) >= theseX[1] & as.numeric(poolLoc[,"Easting"]) <= theseX[2] &
              as.numeric(poolLoc[,"Northing"]) <= theseY[2] & as.numeric(poolLoc[,"Northing"]) >= theseY[1]
            #jLoc_kDeme[2] <- "NGY"
            #(all(jLoc_kDeme!="NGY"))
            #HZdataSetXYdata[,colnames(HZdataSetXYdata)==thisLocusName]
            xCross_start_set <- poolLoc[1,"X_adj"]
            yCross_start_set <- poolLoc[1,"Y_adj"]
            
            if (any(theseRows2D[2:7])==T) {
              poolLoc[which(theseRows2D==T),"XPos_deme"] <- thisDeme_mid_X
              poolLoc[which(theseRows2D==T),"YPos_deme"] <- thisDeme_mid_Y
              # Need to use the most upper left deme from KASP data as starting reference
              # the Easting & Northing was adjusted by min values as follows:
              # poolLoc[1,"XPos"]-min(Xboundaries_midpoint)
              # poolLoc[1,"YPos"]-min(Yboundaries_midpoint)
              poolLoc[which(theseRows2D==T),6] <- poolLoc[which(theseRows2D==T),4] - min(Xboundaries_midpoint)
              poolLoc[which(theseRows2D==T),7] <- poolLoc[which(theseRows2D==T),5] - min(Yboundaries_midpoint)
              
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
              interceptPerp <- ((gradientPerp*(-poolLoc[which(theseRows2D==T),"X_adj"]))+poolLoc[which(theseRows2D==T),"Y_adj"])
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
              poolLoc[which(theseRows2D==T),"DistAlongTransect"] <- as.vector(dist(cbind(xvals,yvals),method="euclidean"))
              #curve(((transectGradient*x)+(transectIntercept)), from=min(DemeMatrix2D_combine_noZero[,"X_adj"]),
              #      to=max(DemeMatrix2D_combine_noZero[,"X_adj"]), col="darkgray",add = TRUE,type = "l",lwd=3,lty=1)
              # An alternative is to use the solve function, gives the same result
              # linearEqn <- matrix(c(transectIntercept,transectGradient,interceptPerp,gradientPerp),2,2,byrow=T)
              # XYcross <- c(-solve(cbind(linearEqn[,2],-1)) %*% linearEqn[,1])
              if (plotPerp==T) {
                lines(x=c(poolLoc[which(theseRows2D==T),"X_adj"],xCross),y=c(poolLoc[which(theseRows2D==T),"Y_adj"],yCross),lwd=1.4,col='black')
                points(x=poolLoc[which(theseRows2D==T),"X_adj"],y=poolLoc[which(theseRows2D==T),"Y_adj"],pch=21,col='black',bg='black',cex=0.8)
                #points(x=xCross,y=yCross,pch=21,col='black',bg='black',cex=0.6,lwd=0.5)
              }
            }
          }
        }
      }
      # poolSeq within bounds of KASP genotypes? Yes
      # min(poolLoc[,"XPos_deme"])>minX; max(poolLoc[,"XPos_deme"])<maxX 
      # min(poolLoc[,"YPos_deme"])>minY; max(poolLoc[,"YPos_deme"])<maxY 
      
    }
   
    #######################################
    # Calculating distance along transect #
    #######################################
    
    #  find the intersection to transect and how far along the transect each deme sits
    distTransect <- matrix(0,nrow(DemeMatrix2D_combine_noZero),2)
    colnames(distTransect) <- c("DistanceRefaxis","DistanceRefaxis_adj")
    # calculate distance along transect
    xCross_start_set <- NULL
    for (thisRow in 1:nrow(distTransect)) {
      # thisRow <- 1
      # 1. find perpendicular linear equation
      # y - y1 = m(x - x1)
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
     
    SampleSize <- (DemeMatrix2D_combine_noZero[,"pCount"]+ DemeMatrix2D_combine_noZero[,"qCount"])/2
    SampleSizeMin <- SampleSize > minSample
    DemeMatrix2D_combine_noZero <- cbind(DemeMatrix2D_combine_noZero[,1:4],round(distTransect[,2],1),DemeMatrix2D_combine_noZero[,5:ncol(DemeMatrix2D_combine_noZero)])
    colnames(DemeMatrix2D_combine_noZero)[5] <- "distAlongTransect"
    
    ################################################
    # plot allele freqs by distance along transect #
    ################################################
    if (noPlots==F) {
      chartLabel <- markerLabelPieChart(thisLocusName,SNPlociListUpdated)
      # re-orientate p (magenta) & q (yellow)
      if (mean(DemeMatrix2D_combine_noZero[1:4,"p"])
          > mean(DemeMatrix2D_combine_noZero[(nrow(DemeMatrix2D_combine_noZero)-4):nrow(DemeMatrix2D_combine_noZero),"p"])) {
        DemeMatrix_temp <- DemeMatrix2D_combine_noZero
        DemeMatrix_temp[,"q"] <- DemeMatrix2D_combine_noZero[,"p"]
        DemeMatrix_temp[,"p"] <- DemeMatrix2D_combine_noZero[,"q"]
        DemeMatrix_temp[,"qCount"] <- DemeMatrix2D_combine_noZero[,"pCount"]
        DemeMatrix_temp[,"pCount"] <- DemeMatrix2D_combine_noZero[,"qCount"]
        #DemeMatrix_temp[,"ppCount"] <- DemeMatrix2D_combine_noZero[,"qqCount"]
        #DemeMatrix_temp[,"qqCount"] <- DemeMatrix2D_combine_noZero[,"ppCount"]
        #DemeMatrix_temp[,"pp"] <- DemeMatrix2D_combine_noZero[,"qq"]
        #DemeMatrix_temp[,"qq"] <- DemeMatrix2D_combine_noZero[,"pp"]
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
      
      plot(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],
           xlab="",ylab="",main=paste(chartLabel,textSec,sep="\n"),ylim=c(0,1),type='n',cex.main=0.9,bty="n",yaxt="n",xaxt="n",
           xlim=c(min(DemeMatrix2D_combine_noZero[,"X_adj"])-1000,max(DemeMatrix2D_combine_noZero[,"X_adj"])+2100))
      axis(1,at=seq(xMin-1000,xMax+1000,1000),labels=NA,col.axis="black",las=2,cex.axis=1.1,lwd=1.6,pos=-0.01, tck = -.02)
      axis(1,at=seq(0,xMax+1000,2000),labels=seq(0,xMax+1000,2000)/1000,col.axis="black",tck=-0.3,las=1,cex.axis=1.1,lwd=0,pos=-0.02)
      axis(2,at=seq(-0.01,1.01,0.1),labels=NA,col.axis="black",las=1,cex.axis=1.1,lwd=1.6,pos=xMin-1000,tck = -.02)
      axis(2,at=seq(0,1.0,0.2),labels=seq(0,1,0.2),col.axis="black",tck=-.03,las=1,cex.axis=1.1,lwd=0,pos=xMin-1005)
      lines(x=c(25000,25000),y=c(0,1),lwd=1.6)
      lines(x=c((xMin-1000),25000),y=c(1,1),lwd=1.6)
      
      if (plotAxisLab==T) {
        mtext("Distance along transect (m)",side=1,cex=1.0,line=2.5)
        mtext("frequency (p)",side=2,cex=1.0,line=2.5)
      }
      points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],DemeMatrix2D_combine_noZero[SampleSizeMin,"p"],
             pch=21,col='black',bg='red')
      points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],DemeMatrix2D_combine_noZero[SampleSizeMin,"q"],
             pch=21,col='black',bg='yellow')
      points(DemeMatrix2D_combine_noZero[SampleSizeMin,"distAlongTransect"],DemeMatrix2D_combine_noZero[SampleSizeMin,"r"],
             pch=21,col='black',bg='blue')
      
      
    }
    
    # Pie plots to file
    # all demes version
    if ("pdf" %in% plotType) {
      # plot all - pdf
      fileDemes <- paste(paste("DemesAll",thisLocusName,sep="_"),".pdf",sep="")
      pdf(fileDemes,width = 16, height = 7,pointsize=18)
      par(mfrow=c(1,1))
      par(mar=c(5,5,2,1))
      
      plot(DemeMatrix2D_combine_noZero[,"X_adj"],DemeMatrix2D_combine_noZero[,"Y_adj"],type='n',
           main=paste(thisLocusName,paste(demeSize,"m demes",sep=""),sep=", "),cex.main=1,bty="n",yaxt="n",xaxt="n",
           xlab = "",ylab = "",col = "black", bg = "black",xaxs = "i", yaxs = "i",xlim=c(xMin-2000,24300),
           ylim=c(yMin-4000,10050))
      axis(1,at=seq(xMin-2000,xMax,1000),labels=NA,col.axis="black",las=2,cex.axis=1.1,lwd=1.6,pos=yMin-4000, tck = -.02)
      axis(1,at=seq(-2000,xMax,2000),labels=seq(xMin-2000,xMax,2000)/1000,col.axis="black",tck=-0.3,las=1,cex.axis=1.1,lwd=0,pos=yMin-4000.5)
      axis(2,at=seq(yMin-4000,yMax+3800,1000),labels=NA,col.axis="black",las=1,cex.axis=1.1,lwd=1.6,pos=xMin-2000,tck = -.02)
      axis(2,at=seq(yMin-4000,yMax+3800,2000),labels=seq(yMin-4000,yMax+3800,2000)/1000,col.axis="black",tck=-.03,las=1,cex.axis=1.1,lwd=0,pos=xMin-2000.5)
      lines(x=c(24000,24000),y=c(yMin-4000,yMax+3000),lwd=1.6)
      lines(x=c(xMin-2000,24000),y=c(yMax+3000,yMax+3000),lwd=1.6)
      
      curve(((transectGradient*x)+(transectIntercept)), from=c(min(DemeMatrix2D_combine_noZero[,"X_adj"])-1000,
                                                               to=max(DemeMatrix2D_combine_noZero[,"X_adj"])+1000), col="darkgray",add = TRUE,type = "l",lwd=2,lty=1)
      mtext("East (km)",side=1,cex=1,line=2.5)
      mtext("North (km)",side=2,cex=1,line=2.5)
      # add lines to plot
      if (plotPerp==T) {
        for (thisRow in 1:nrow(distTransect)) {
          # thisRow <- 2
          # 1. find perpendicular linear equation
          # y - y1 = m(x - x1)
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
            #HZdataSetXYdata[,colnames(HZdataSetXYdata)==thisLocusName]
            #xCross_start_set <- poolLoc[1,"X_adj"]
            #yCross_start_set <- poolLoc[1,"Y_adj"]
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
          if (plotPerp==T) {
            lines(x=c(thisX,xCross),y=c(thisY,yCross),lwd=1.0,col='blue')
            #points(x=xCross,y=yCross,pch=21,col='red',cex=0.6,lwd=0.5)
          }
        }
      }
      if (addPie==T) {
        # plotting pie charts
        # minSample <- 5
        totalN <- (DemeMatrix2D_combine_noZero[,"pCount"]+DemeMatrix2D_combine_noZero[,"qCount"]+DemeMatrix2D_combine_noZero[,"rCount"])/2
        keep <- which(totalN>=minSample)
        DemeMatrix2D_combine_noZero <- DemeMatrix2D_combine_noZero[keep,]
        allthisP <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[,"pCount"]))
        allthisQ <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[,"qCount"]))
        allthisR <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[,"rCount"]))
        allSum <- allthisP+allthisQ+allthisR
        allFreq <- allthisP/allSum
        if (mean(allFreq[1:4])>0.6) {
          allthisP_temp <- allthisQ
          allthisQ_temp <- allthisP
          allthisQ <- allthisP_temp
          allthisP <- allthisQ_temp
        }
        for (thisDeme2D in 1:nrow(DemeMatrix2D_combine_noZero)) {
          # thisDeme2D <- 1
          thisXraw <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[thisDeme2D,"X_adj"]))
          thisYraw <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[thisDeme2D,"Y_adj"]))
          thisP <- as.numeric(allthisP[thisDeme2D])
          thisQ <- as.numeric(allthisQ[thisDeme2D])
          thisR <- as.numeric(allthisR[thisDeme2D])
          
          # pie plot size has to be tweaked depending on deme size
          # small pie plots, original was radius 69 for deme size 200
          #     if (sum(thisP,thisQ)>minNum) {
          #       add.pie(z=c(thisP,thisQ), x=thisXmid, y=thisYmid, radius=69,col=c((rgb(255,0,0,255,maxColorValue=255)),(rgb(255,255,0,255,maxColorValue=255))),labels="")
          #     }
          # larger pie plots - 120 works okay for deme size 400
          totalN <- thisP + thisQ + thisR
          #if ((totalN < minNum)==F) {
          add.pie(z=c(thisP,thisR,thisQ), x=thisXraw, y=thisYraw, radius=85,col=c((rgb(255,0,0,255,maxColorValue=255)),"blue",(rgb(255,255,0,255,maxColorValue=255))),labels="")
          #}
        }
      }
      if (poolSeqPlot==T & thisLocus==1) {
        # setup matrix for poolSeq comparison
        # 1. YP4 ??? 52 random plants from yellow parapatric
        # 2. yelPool20 ??? 20 yellow plants between yellow flank and parapatric region (individuals J1141-J1160)
        # 3. YP2 ??? 50 random plants from yellow lower flank
        # 4. YP1 ??? 50 random plants from yellow diagonal flank
        # 5. MP2 - 50 random plants from magenta lower flank
        # 6. MP4 - 50 random plants from magenta lower flank
        # 7. mgtPool20 ??? 20 magenta plants between magenta flank and parapatric region (individuals J1121-J1140)
        # 8. MP11 ??? 50 random plants from magenta parapatric
        
        poolLoc <- matrix(0,7,8)
        poolLoc <- as.data.frame(poolLoc)
        colnames(poolLoc) <- c("pool","Easting","Northing","XPos_deme","YPos_deme","X_adj","Y_adj","DistAlongTransect")
        poolLoc[1:7,1] <- c("DemeKASP_Start","YP4","YP2","YP1","MP2","MP4","MP11")
        poolLoc[1,4:5] <- DemeMatrix2D_combine_noZero[1,1:2]
        poolLoc[1,6:7] <- DemeMatrix2D_combine_noZero[1,3:4]
        poolLoc[2,2:3] <- c(411604.6046997,4690386.0014007)
        poolLoc[3,2:3] <- c(421903.0320033,4686484.3909220)
        poolLoc[4,2:3] <- c(422119.7797325,4686390.2520109)
        poolLoc[5,2:3] <- c(424462.3465657,4686112.3149515)
        poolLoc[6,2:3] <- c(425146.6755758,4686021.5773514)
        poolLoc[7,2:3] <- c(431621.6848923,4686847.1448457)
        poolLoc[,2] <- as.numeric(as.vector(poolLoc[,2]))
        poolLoc[,3] <- as.numeric(as.vector(poolLoc[,3]))
        poolLoc[,4] <- as.numeric(as.vector(poolLoc[,4]))
        poolLoc[,5] <- as.numeric(as.vector(poolLoc[,5]))
        poolLoc[,6] <- as.numeric(as.vector(poolLoc[,6]))
        poolLoc[,7] <- as.numeric(as.vector(poolLoc[,7]))
        Deme2D_pools <- DemeMatrix2D
        Deme2D_pools[,] <- -9
        
        # find pool locations and distance along transect
        # 1. YP4 ??? 52 random plants from yellow parapatric
        # 2. yelPool20 ??? 20 yellow plants between yellow flank and parapatric region (individuals J1141-J1160)
        # 3. YP2 ??? 50 random plants from yellow lower flank
        # 4. YP1 ??? 50 random plants from yellow diagonal flank
        # 5. MP2 - 50 random plants from magenta lower flank
        # 6. MP4 - 50 random plants from magenta lower flank
        # 7. mgtPool20 ??? 20 magenta plants between magenta flank and parapatric region (individuals J1121-J1140)
        # 8. MP11 ??? 50 random plants from magenta parapatric
        # find locations for poolSeq data (no genotype calculations)
        demeCounter <- 0
        theseRows2D <- NULL
        for (thisXcol in 1:ncol(Deme2D_pools)) {
          for (thisYrow in 1:nrow(Deme2D_pools)) {
            demeCounter <- demeCounter + 1
            # thisXcol <- 7
            # thisYrow <- 5
            thisDeme_mid_X <- as.numeric(colnames(Deme2D_pools)[thisXcol])
            thisDeme_mid_Y <- as.numeric(rownames(Deme2D_pools)[thisYrow])
            
            theseX <- c(as.numeric(colnames(Deme2D_pools)[thisXcol])-(demeSize/2),as.numeric(colnames(Deme2D_pools)[thisXcol])+(demeSize/2))
            theseY <- c(as.numeric(rownames(Deme2D_pools)[thisYrow])-(demeSize/2),as.numeric(rownames(Deme2D_pools)[thisYrow])+(demeSize/2))
            theseRows2D <- as.numeric(poolLoc[,"Easting"]) >= theseX[1] & as.numeric(poolLoc[,"Easting"]) <= theseX[2] &
              as.numeric(poolLoc[,"Northing"]) <= theseY[2] & as.numeric(poolLoc[,"Northing"]) >= theseY[1]
            #jLoc_kDeme[2] <- "NGY"
            #(all(jLoc_kDeme!="NGY"))
            #HZdataSetXYdata[,colnames(HZdataSetXYdata)==thisLocusName]
            if (any(theseRows2D[2:7])==T) {
              poolLoc[which(theseRows2D==T),"XPos_deme"] <- thisDeme_mid_X
              poolLoc[which(theseRows2D==T),"YPos_deme"] <- thisDeme_mid_Y
              # Need to use the most upper left deme from KASP data as starting reference
              # the Easting & Northing was adjusted by min values as follows:
              # poolLoc[1,"XPos"]-min(Xboundaries_midpoint)
              # poolLoc[1,"YPos"]-min(Yboundaries_midpoint)
              poolLoc[which(theseRows2D==T),6] <- poolLoc[which(theseRows2D==T),4] - min(Xboundaries_midpoint)
              poolLoc[which(theseRows2D==T),7] <- poolLoc[which(theseRows2D==T),5] - min(Yboundaries_midpoint)
              
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
              interceptPerp <- ((gradientPerp*(-poolLoc[which(theseRows2D==T),"X_adj"]))+poolLoc[which(theseRows2D==T),"Y_adj"])
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
              poolLoc[which(theseRows2D==T),"DistAlongTransect"] <- as.vector(dist(cbind(xvals,yvals),method="euclidean"))
              #curve(((transectGradient*x)+(transectIntercept)), from=min(DemeMatrix2D_combine_noZero[,"X_adj"]),
              #      to=max(DemeMatrix2D_combine_noZero[,"X_adj"]), col="darkgray",add = TRUE,type = "l",lwd=3,lty=1)
              # An alternative is to use the solve function, gives the same result
              # linearEqn <- matrix(c(transectIntercept,transectGradient,interceptPerp,gradientPerp),2,2,byrow=T)
              # XYcross <- c(-solve(cbind(linearEqn[,2],-1)) %*% linearEqn[,1])
              if (plotPerp==T) {
                lines(x=c(poolLoc[which(theseRows2D==T),"X_adj"],xCross),y=c(poolLoc[which(theseRows2D==T),"Y_adj"],yCross),lwd=1.4,col='black')
                points(x=poolLoc[which(theseRows2D==T),"X_adj"],y=poolLoc[which(theseRows2D==T),"Y_adj"],pch=21,col='black',bg='black',cex=0.7)
                #points(x=xCross,y=yCross,pch=21,col='black',bg='darkgrey',cex=0.7,lwd=0.9)
              }
            }
          }
        }
      }
      if (addPie==F) {
        points(DemeMatrix2D_combine_noZero[,"X_adj"],DemeMatrix2D_combine_noZero[,"Y_adj"])
      }
      dev.off()
    }
    # zoomed version
    if ("pdf" %in% plotType) {
      # plot zoomed in area
      fileDemes <- paste(paste("DemesZoom",thisLocusName,sep="_"),".pdf",sep="")
      pdf(fileDemes,width = 16, height = 8.2,pointsize=18)
      par(mfrow=c(1,1))
      par(mar=c(5,5,2,1))
      #xMin <- min(DemeMatrix2D_combine_noZero[,"X_adj"])
      #xMax <- max(DemeMatrix2D_combine_noZero[,"X_adj"])
      #yMin <- min(DemeMatrix2D_combine_noZero[,"Y_adj"])
      #yMax <- max(DemeMatrix2D_combine_noZero[,"Y_adj"])
      plot(DemeMatrix2D_combine_noZero[,"X_adj"],DemeMatrix2D_combine_noZero[,"Y_adj"],type='n',
           main=paste(thisLocusName,paste(demeSize,"m demes",sep=""),sep=", "),cex.main=1,bty="n",yaxt="n",xaxt="n",
           xlab = "",ylab = "",col = "black", bg = "black",xaxs = "i", yaxs = "i",xlim=c(7500,20005),
           ylim=c(-2000,5010))
      axis(1,at=seq(8000,20000,1000),labels=NA,col.axis="black",las=2,cex.axis=1.1,lwd=1.6,pos=-2000, tck = -.02)
      axis(1,at=seq(8000,20000,2000),labels=seq(8000,20000,2000)/1000,col.axis="black",tck=-0.3,las=1,cex.axis=1.1,lwd=0,pos=-2000.5)
      axis(2,at=seq(-2000,5000,1000),labels=NA,col.axis="black",las=1,cex.axis=1.1,lwd=1.6,pos=8000,tck = -.02)
      axis(2,at=seq(-2000,5000,2000),labels=seq(-2000,5000,2000)/1000,col.axis="black",tck=-.03,las=1,cex.axis=1.1,lwd=0,pos=7999)
      lines(x=c(20000,20000),y=c(-2000,5000),lwd=1.6)
      lines(x=c(8000,20000),y=c(5000,5000),lwd=1.6)
      curve(((transectGradient*x)+(transectIntercept)), from=c(8000),to=20000, col="darkgray",add = TRUE,type = "l",lwd=2,lty=1)
      mtext("East (km)",side=1,cex=1,line=2.5)
      mtext("North (km)",side=2,cex=1,line=2.5)
      # add lines to plot
      if (plotPerp==T) {
        for (thisRow in 1:nrow(distTransect)) {
          # thisRow <- 2
          # 1. find perpendicular linear equation
          # y - y1 = m(x - x1)
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
          if (plotPerp==T) {
            lines(x=c(thisX,xCross),y=c(thisY,yCross),lwd=0.8,col='blue')
            #points(x=xCross,y=yCross,pch=21,col='red',cex=0.6,lwd=0.5)
          }
        }
      }
      if (addPie==T) {
        # plotting pie charts
        # minSample <- 5
        totalN <- (DemeMatrix2D_combine_noZero[,"pCount"]+DemeMatrix2D_combine_noZero[,"qCount"]+DemeMatrix2D_combine_noZero[,"rCount"])/2
        keep <- which(totalN>=minSample)
        DemeMatrix2D_combine_noZero <- DemeMatrix2D_combine_noZero[keep,]
        allthisP <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[,"pCount"]))
        allthisQ <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[,"qCount"]))
        allthisR <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[,"rCount"]))
        allSum <- allthisP+allthisQ+allthisR
        allFreq <- allthisP/allSum
        if (mean(allFreq[1:4])>0.6) {
          allthisP_temp <- allthisQ
          allthisQ_temp <- allthisP
          allthisQ <- allthisP_temp
          allthisP <- allthisQ_temp
        }
        for (thisDeme2D in 1:nrow(DemeMatrix2D_combine_noZero)) {
          # thisDeme2D <- 1
          thisXraw <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[thisDeme2D,"X_adj"]))
          thisYraw <- as.numeric(as.vector(DemeMatrix2D_combine_noZero[thisDeme2D,"Y_adj"]))
          thisP <- as.numeric(allthisP[thisDeme2D])
          thisQ <- as.numeric(allthisQ[thisDeme2D])
          thisR <- as.numeric(allthisR[thisDeme2D])
          
          # pie plot size has to be tweaked depending on deme size
          # small pie plots, original was radius 69 for deme size 200
          #     if (sum(thisP,thisQ)>minNum) {
          #       add.pie(z=c(thisP,thisQ), x=thisXmid, y=thisYmid, radius=69,col=c((rgb(255,0,0,255,maxColorValue=255)),(rgb(255,255,0,255,maxColorValue=255))),labels="")
          #     }
          # larger pie plots - 120 works okay for deme size 400
          totalN <- thisP + thisQ + thisR
          #if ((totalN < minNum)==F) {
          add.pie(z=c(thisP,thisR,thisQ), x=thisXraw, y=thisYraw, radius=85,col=c((rgb(255,255,0,255,maxColorValue=255)),"blue",(rgb(255,0,0,255,maxColorValue=255))),labels="")
          #}
        }
      }
      if (poolSeqPlot==T & thisLocus==1) {
        # setup matrix for poolSeq comparison
        # 1. YP4 ??? 52 random plants from yellow parapatric
        # 2. yelPool20 ??? 20 yellow plants between yellow flank and parapatric region (individuals J1141-J1160)
        # 3. YP2 ??? 50 random plants from yellow lower flank
        # 4. YP1 ??? 50 random plants from yellow diagonal flank
        # 5. MP2 - 50 random plants from magenta lower flank
        # 6. MP4 - 50 random plants from magenta lower flank
        # 7. mgtPool20 ??? 20 magenta plants between magenta flank and parapatric region (individuals J1121-J1140)
        # 8. MP11 ??? 50 random plants from magenta parapatric
        
        poolLoc <- matrix(0,7,8)
        poolLoc <- as.data.frame(poolLoc)
        colnames(poolLoc) <- c("pool","Easting","Northing","XPos_deme","YPos_deme","X_adj","Y_adj","DistAlongTransect")
        poolLoc[1:7,1] <- c("DemeKASP_Start","YP4","YP2","YP1","MP2","MP4","MP11")
        poolLoc[1,4:5] <- DemeMatrix2D_combine_noZero[1,1:2]
        poolLoc[1,6:7] <- DemeMatrix2D_combine_noZero[1,3:4]
        poolLoc[2,2:3] <- c(411604.6046997,4690386.0014007)
        poolLoc[3,2:3] <- c(421903.0320033,4686484.3909220)
        poolLoc[4,2:3] <- c(422119.7797325,4686390.2520109)
        poolLoc[5,2:3] <- c(424462.3465657,4686112.3149515)
        poolLoc[6,2:3] <- c(425146.6755758,4686021.5773514)
        poolLoc[7,2:3] <- c(431621.6848923,4686847.1448457)
        poolLoc[,2] <- as.numeric(as.vector(poolLoc[,2]))
        poolLoc[,3] <- as.numeric(as.vector(poolLoc[,3]))
        poolLoc[,4] <- as.numeric(as.vector(poolLoc[,4]))
        poolLoc[,5] <- as.numeric(as.vector(poolLoc[,5]))
        poolLoc[,6] <- as.numeric(as.vector(poolLoc[,6]))
        poolLoc[,7] <- as.numeric(as.vector(poolLoc[,7]))
        Deme2D_pools <- DemeMatrix2D
        Deme2D_pools[,] <- -9
        
        # find pool locations and distance along transect
        # 1. YP4 ??? 52 random plants from yellow parapatric
        # 2. yelPool20 ??? 20 yellow plants between yellow flank and parapatric region (individuals J1141-J1160)
        # 3. YP2 ??? 50 random plants from yellow lower flank
        # 4. YP1 ??? 50 random plants from yellow diagonal flank
        # 5. MP2 - 50 random plants from magenta lower flank
        # 6. MP4 - 50 random plants from magenta lower flank
        # 7. mgtPool20 ??? 20 magenta plants between magenta flank and parapatric region (individuals J1121-J1140)
        # 8. MP11 ??? 50 random plants from magenta parapatric
        # find locations for poolSeq data (no genotype calculations)
        demeCounter <- 0
        theseRows2D <- NULL
        for (thisXcol in 1:ncol(Deme2D_pools)) {
          for (thisYrow in 1:nrow(Deme2D_pools)) {
            demeCounter <- demeCounter + 1
            # thisXcol <- 7
            # thisYrow <- 5
            thisDeme_mid_X <- as.numeric(colnames(Deme2D_pools)[thisXcol])
            thisDeme_mid_Y <- as.numeric(rownames(Deme2D_pools)[thisYrow])
            
            theseX <- c(as.numeric(colnames(Deme2D_pools)[thisXcol])-(demeSize/2),as.numeric(colnames(Deme2D_pools)[thisXcol])+(demeSize/2))
            theseY <- c(as.numeric(rownames(Deme2D_pools)[thisYrow])-(demeSize/2),as.numeric(rownames(Deme2D_pools)[thisYrow])+(demeSize/2))
            theseRows2D <- as.numeric(poolLoc[,"Easting"]) >= theseX[1] & as.numeric(poolLoc[,"Easting"]) <= theseX[2] &
              as.numeric(poolLoc[,"Northing"]) <= theseY[2] & as.numeric(poolLoc[,"Northing"]) >= theseY[1]
            #jLoc_kDeme[2] <- "NGY"
            #(all(jLoc_kDeme!="NGY"))
            #HZdataSetXYdata[,colnames(HZdataSetXYdata)==thisLocusName]
            if (any(theseRows2D[2:7])==T) {
              poolLoc[which(theseRows2D==T),"XPos_deme"] <- thisDeme_mid_X
              poolLoc[which(theseRows2D==T),"YPos_deme"] <- thisDeme_mid_Y
              # Need to use the most upper left deme from KASP data as starting reference
              # the Easting & Northing was adjusted by min values as follows:
              # poolLoc[1,"XPos"]-min(Xboundaries_midpoint)
              # poolLoc[1,"YPos"]-min(Yboundaries_midpoint)
              poolLoc[which(theseRows2D==T),6] <- poolLoc[which(theseRows2D==T),4] - min(Xboundaries_midpoint)
              poolLoc[which(theseRows2D==T),7] <- poolLoc[which(theseRows2D==T),5] - min(Yboundaries_midpoint)
              
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
              interceptPerp <- ((gradientPerp*(-poolLoc[which(theseRows2D==T),"X_adj"]))+poolLoc[which(theseRows2D==T),"Y_adj"])
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
              poolLoc[which(theseRows2D==T),"DistAlongTransect"] <- as.vector(dist(cbind(xvals,yvals),method="euclidean"))
              #curve(((transectGradient*x)+(transectIntercept)), from=min(DemeMatrix2D_combine_noZero[,"X_adj"]),
              #      to=max(DemeMatrix2D_combine_noZero[,"X_adj"]), col="darkgray",add = TRUE,type = "l",lwd=3,lty=1)
              # An alternative is to use the solve function, gives the same result
              # linearEqn <- matrix(c(transectIntercept,transectGradient,interceptPerp,gradientPerp),2,2,byrow=T)
              # XYcross <- c(-solve(cbind(linearEqn[,2],-1)) %*% linearEqn[,1])
              if (plotPerp==T) {
                lines(x=c(poolLoc[which(theseRows2D==T),"X_adj"],xCross),y=c(poolLoc[which(theseRows2D==T),"Y_adj"],yCross),lwd=1.4,col='black')
                points(x=poolLoc[which(theseRows2D==T),"X_adj"],y=poolLoc[which(theseRows2D==T),"Y_adj"],pch=21,col='black',bg='black',cex=0.7)
                #points(x=xCross,y=yCross,pch=21,col='black',bg='darkgrey',cex=0.7,lwd=0.9)
              }
            }
          }
        }
      }
      if (addPie==F) {
        points(DemeMatrix2D_combine_noZero[,"X_adj"],DemeMatrix2D_combine_noZero[,"Y_adj"])
      }
      dev.off()
    }
    
    # Send data to list
    #totalSampleSize2D <- (sum(DemeMatrix2D_combine_noZero[,"pCount"])+sum(DemeMatrix2D_combine_noZero[,"qCount"]))/2
    # send all 2D data to list
    DemeData[[thisLocusName]][["DemeMatrix2D_pCount"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_qCount"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_rCount"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_Ngenes"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_combine_noZero"]] <- list()
    
    DemeData[[thisLocusName]][["DemeMatrix2D_FLA_FLA_Count"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_FLAr_FLA_Count"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_FLAr_FLAr_Count"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_FLA_fla_Count"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_FLAr_fla_Count"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_FLAr_fla_Count"]] <- list()
    
    DemeData[[thisLocusName]][["DemeMatrix2D_FLA_FLA_freq"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_FLAr_FLA_freq"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_FLAr_FLAr_freq"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_FLA_fla_freq"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_FLAr_fla_freq"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_FLAr_fla_freq"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_Individuals"]] <- list()
    DemeData[[thisLocusName]][["DemeMatrix2D_DistanceRefLine"]] <- list()
    
    DemeData[[thisLocusName]][["DemeMatrix2D_pCount"]] <- DemeMatrix2D_pCount
    DemeData[[thisLocusName]][["DemeMatrix2D_qCount"]] <- DemeMatrix2D_qCount
    DemeData[[thisLocusName]][["DemeMatrix2D_rCount"]] <- DemeMatrix2D_rCount
    DemeData[[thisLocusName]][["DemeMatrix2D_Ngenes"]] <- DemeMatrix2D_Ngenes
    DemeData[[thisLocusName]][["DemeMatrix2D_combine_noZero"]] <- DemeMatrix2D_combine_noZero
    
    DemeData[[thisLocusName]][["DemeMatrix2D_FLA_FLA_Count"]] <- DemeMatrix2D_FLA_FLA_Count
    DemeData[[thisLocusName]][["DemeMatrix2D_FLAr_FLA_Count"]] <- DemeMatrix2D_FLAr_FLA_Count
    DemeData[[thisLocusName]][["DemeMatrix2D_FLAr_FLAr_Count"]] <- DemeMatrix2D_FLAr_FLAr_Count
    DemeData[[thisLocusName]][["DemeMatrix2D_FLA_fla_Count"]] <- DemeMatrix2D_FLA_fla_Count
    DemeData[[thisLocusName]][["DemeMatrix2D_FLAr_fla_Count"]] <- DemeMatrix2D_FLAr_fla_Count
    DemeData[[thisLocusName]][["DemeMatrix2D_FLAr_fla_Count"]] <- DemeMatrix2D_FLAr_fla_Count
    
    DemeData[[thisLocusName]][["DemeMatrix2D_FLA_FLA_freq"]] <- DemeMatrix2D_FLA_FLA
    DemeData[[thisLocusName]][["DemeMatrix2D_FLAr_FLA_freq"]] <- DemeMatrix2D_FLAr_FLA
    DemeData[[thisLocusName]][["DemeMatrix2D_FLAr_FLAr_freq"]] <- DemeMatrix2D_FLAr_FLAr
    DemeData[[thisLocusName]][["DemeMatrix2D_FLA_fla_freq"]] <- DemeMatrix2D_FLA_fla
    DemeData[[thisLocusName]][["DemeMatrix2D_FLAr_fla_freq"]] <- DemeMatrix2D_FLAr_fla
    DemeData[[thisLocusName]][["DemeMatrix2D_FLAr_fla_freq"]] <- DemeMatrix2D_FLAr_fla
    DemeData[[thisLocusName]][["DemeMatrix2D_Individuals"]] <- DemeIndividuals
    DemeData[[thisLocusName]][["DemeMatrix2D_DistanceRefLine"]] <- distTransect
    
  }  # end locus loop
  return(DemeData)
}
