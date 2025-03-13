multiClinePlot_gplots <- function(thisData,numColFacet,hapData) {
  # thisData <- clineFits_test
  # numColFacet <- 5
  #rownames(thisData) <- as.vector(thisData[,"numberPlot"])
  # Sigmoid curve function
  thisData <- thisData[thisData[,"N_total"]!=0,]
  myCurveModelFinalFigs <- function(XCordinate) {
    ((this_pL + ((this_pR-this_pL)/(1 + exp((-4*(XCordinate-thisCentre))/thisGradient)))))
  }
  # functions for expected line
  explineVals <- function(theseParams) {
    XtoEval <- seq(0,25000,100)
    expVals <- matrix(0,length(XtoEval),2)
    expVals[,1] <- XtoEval
    expVals[,2] <- myCurveModelFinalFigs(XtoEval)
    return(expVals)
  }
  explineValsDemes <- function(theseParams,thisDemeDataEval) {
    # thisDemeDataEval <- locusClineMatrix_p
    XtoEval <- as.numeric(as.vector(thisDemeDataEval[,"Xbin"]))
    expVals <- matrix(0,length(XtoEval),2)
    expVals[,1] <- XtoEval
    expVals[,2] <- myCurveModelFinalFigs(XtoEval)
    return(expVals)
  }
  explineValsDemesLinearReg <- function(theseParams,thisDemeDataEval) {
    # thisDemeDataEval <- locusClineMatrix_p
    XtoEval <- as.numeric(as.vector(thisDemeDataEval[,"Xbin"]))
    expVals <- matrix(0,length(XtoEval),2)
    expVals[,1] <- XtoEval
    thisGradient <- theseParams[1]
    thisIntercept <- theseParams[2]
    myCurveModellm <- function(XCordinate) {
      ((thisGradient*XCordinate) + thisIntercept)
    }
    expVals[,2] <- myCurveModellm(XtoEval)
    return(expVals)
  }
  getMarkerLabelNew <- function(thisLocusName,summaryClines,SampleSize) {
    # SampleSize <- 1514
    # summaryClines <- summaryModelsAlllociFigs_bothRoads
    markerType <- as.vector(summaryClines[summaryClines[,"locusName"]==thisLocusName,"Type"])
    thisDataRow <- as.vector(summaryClines[summaryClines[,"locusName"]==thisLocusName,]) 
    markerLabelMake <- NULL
    if (markerType=="Parentage_nDNA" | markerType=="Parentage_cpDNA") {
      markerLabelMake <- "P"
    }
    if (markerType=="Diagnostic") {
      markerLabelMake <- "D"
    }
    if (markerType=="Diagnostic_cpDNA") {
      markerLabelMake <- "D"
    }  
    if (markerType=="Diagnostic_ROSel") {
      markerLabelMake <- "D_ROSel"
    }
    if (markerType=="SULF") {
      markerLabel <- "SULF"
    }
    if (markerType=="DEF") {
      markerLabel <- "DEF"
    }
    if (markerType=="DICH") {
      markerLabelMake <- "DICH"
    }
    if (markerType=="VENOSA") {
      markerLabelMake <- "VENOSA"
    }
    if (markerType=="MIXTA") {
      markerLabelMake <- "MIXTA"
    }
    if (markerType=="FLO") {
      markerLabelMake <- "FLO"
    }
    if (markerType=="EL") {
      markerLabelMake <- "EL"
    }
    if (markerType=="barrier") {
      markerLabelMake <- "barrier"
    }
    if (markerType=="ROS1") {
      markerLabelMake <- "ROS1"
    }
    if (markerType=="ROS2") {
      markerLabelMake <- "ROS2"
    }
    if (markerType=="ROS3") {
      markerLabelMake <- "ROS3"
    }
    if (markerType=="Sulf_linked") {
      markerLabelMake <- "SULF_LK"
    }
    SampleLabel <- NULL
    SampleLabel <- paste("n=",SampleSize,sep="")
    markerLabelMake <- paste(NumCline,paste(markerLabelMake,paste("LG",(thisDataRow["LG"]),sep=""),paste(thisDataRow["cM"],"cM",sep=""),SampleLabel,sep=";"),sep=", ")
    return(markerLabelMake)
  } 
  # first pass to collect all the possible Xbins
  XbinCollect <- NULL
  for (thisLocus in 1:nrow(thisData)) {
    # thisLocus <- 1
    XbinCollect <- c(XbinCollect,as.numeric(as.vector(unlist(strsplit(as.vector(thisData[thisLocus,"Xvals"]),";")))))
    # locusClineMatrix[,"Xbin"]
  }
  Xcordinates_full <- unique(XbinCollect)
  # insert some extra Xbins to ensure smooth curves
  Xcordinates_full <- c(Xcordinates_full,seq(200,800,100),seq(1000,1800,100),seq(1900,2600,100),seq(2800,13000,100))
  # setup matrices for multifacet plots
  locusClineMatrix <- matrix(NA,nrow=length(Xcordinates_full),ncol=(nrow(thisData)+1))
  colnames(locusClineMatrix) <- c("Xbin",as.vector(thisData[,"locusName"]))
  locusClineMatrix[,"Xbin"] <- Xcordinates_full
  locusClineMatrix <- locusClineMatrix[order(locusClineMatrix[,"Xbin"]),]
  # Now makes separate matrices for different data
  locusClineMatrix_p <- as.data.frame(locusClineMatrix) # allele freq p
  locusClineMatrix_N <- locusClineMatrix_p # N samples
  locusClineMatrix_exp <- locusClineMatrix_p # exp based on best cline fit
  locusClineMatrix_col <- locusClineMatrix_p # colour of curve
  locusClineMatrix_p_logit <- locusClineMatrix_p # allele freq p logit scale
  locusClineMatrix_exp_logit <- locusClineMatrix_p # exp freq based on best cline fit (logit scale)
  # send data to new matrices
  clineDetailHeaderAll <- NULL
  locusHeaderAll <- NULL
  for (thisLocus in 1:nrow(thisData)) {
    # thisLocus <- 1
    locusName <- as.vector(thisData[thisLocus,"locusName"])
    #position_thisScaff <- as.vector(thisData[thisLocus,"position_thisScaff"])
    #scaffold <- as.vector(thisData[thisLocus,"scaffold"])
    LG <- as.vector(thisData[thisLocus,"LG"])
    Xbins <- as.numeric(as.vector(unlist(strsplit(as.vector(thisData[thisLocus,"Xvals"]),";"))))
    p_adj <- as.numeric(as.vector(unlist(strsplit(as.vector(thisData[thisLocus,"p"]),";"))))
    Nind <- as.numeric(as.vector(unlist(strsplit(as.vector(thisData[thisLocus,"Ngenes"]),";"))))/2
    p_logit <- log(p_adj/(1-p_adj))
    locusClineMatrix_p[match(Xbins,locusClineMatrix[,"Xbin"]),colnames(locusClineMatrix)==locusName] <- p_adj
    locusClineMatrix_N[match(Xbins,locusClineMatrix[,"Xbin"]),colnames(locusClineMatrix)==locusName] <- Nind
    locusClineMatrix_p_logit[match(Xbins,locusClineMatrix[,"Xbin"]),colnames(locusClineMatrix)==locusName] <- p_logit
    locusClineMatrix_col[match(Xbins,locusClineMatrix[,"Xbin"]),colnames(locusClineMatrix)==locusName] <- "gray"
    #plot(locusClineMatrix[,"Xbin"],locusClineMatrix_p[,"s203_228663"])
    if (thisData[thisLocus,"centre"]==0) {
      thisCentre <- -9
      thisGradient <- 13000
      this_pL <- -9
      this_pR <- +9
      thisLmax <- 0
      clineDetailHeader <- "no fit"
      if (LG==4) {
        locusHeader <- paste(thisData[thisLocus,"numberPlot"],locusName,sep=": ")
      }
      if (scaffold=="ros_assembly") {
        locusHeader <- paste(thisData[thisLocus,"numberPlot"],position_thisScaff,sep=": ")
      }
    }
    if (thisData[thisLocus,"centre"]!=0) {
      thisCentre <- as.numeric(as.vector(thisData[thisLocus,"centre"]))
      thisGradient <- as.numeric(as.vector(thisData[thisLocus,"width"]))
      this_pL <- as.numeric(as.vector(thisData[thisLocus,"pL"]))
      this_pR <- as.numeric(as.vector(thisData[thisLocus,"pR"]))
      locusClineMatrix_exp[,colnames(locusClineMatrix)==locusName] <- explineValsDemes(c(thisCentre,thisGradient,this_pL,this_pR),locusClineMatrix_p)[,2]
      locusClineMatrix_p_logit[,colnames(locusClineMatrix)==locusName] <- log(locusClineMatrix_exp[,colnames(locusClineMatrix)==locusName]/(1-locusClineMatrix_exp[,colnames(locusClineMatrix)==locusName]))
      locusClineMatrix_exp_logit[,colnames(locusClineMatrix)==locusName] <- log(locusClineMatrix_exp[,colnames(locusClineMatrix)==locusName]/(1-locusClineMatrix_exp[,colnames(locusClineMatrix)==locusName]))
      clineDetailHeader <-paste(c(paste("c=", round(thisCentre,0),sep=""),
                        paste("w=",round(thisGradient,0),sep=""),
                        paste("pL=",round(this_pL,1),sep=""),
                        paste("pR=",round(this_pR,1))),collapse=",")
      # if (LG==4) {
      #   locusHeader <- paste(thisData[thisLocus,"numberPlot"],locusName,sep=": ")
      # }
      # if (scaffold=="ros_assembly") {
      #   locusHeader <- paste(thisData[thisLocus,"numberPlot"],position_thisScaff,sep=": ")
      # }
    }
    clineDetailHeaderAll <- c(clineDetailHeaderAll,clineDetailHeader)
    locusHeaderAll <- c(locusHeaderAll,locusName)
  }

  # adjust column names for figures
  colnames(locusClineMatrix_p) <- c("Xbin",locusHeaderAll)
  colnames(locusClineMatrix_N) <- c("Xbin",locusHeaderAll)
  colnames(locusClineMatrix_exp) <- c("Xbin",locusHeaderAll)
  colnames(locusClineMatrix_col) <- c("Xbin",locusHeaderAll)
  colnames(locusClineMatrix_p_logit) <- c("Xbin",locusHeaderAll)
  colnames(locusClineMatrix_exp_logit) <- c("Xbin",locusHeaderAll)
  
  # reshape matrices 
  locusClineMatrix_p_melt <- melt(locusClineMatrix_p,id.vars="Xbin")
  colnames(locusClineMatrix_p_melt)[2:3] <- c("locus","p")
  locusClineMatrix_p_exp_melt <- melt(locusClineMatrix_exp,id.vars="Xbin")
  colnames(locusClineMatrix_p_exp_melt)[2:3] <- c("locus","p_exp")
  locusClineMatrix_p_logit_melt <- melt(locusClineMatrix_p_logit,id.vars="Xbin")
  colnames(locusClineMatrix_p_logit_melt)[2:3] <- c("locus","p_logit")
  locusClineMatrix_exp_logit_melt <- melt(locusClineMatrix_exp_logit,id.vars="Xbin")
  colnames(locusClineMatrix_exp_logit_melt)[2:3] <- c("locus","p_logit_exp")
  locusClineMatrix_col_melt <- melt(locusClineMatrix_col,id.vars="Xbin")
  colnames(locusClineMatrix_col_melt)[2:3] <- c("locus","col")
  
  # combine matrices
  # check locus names in same spot
  #length(locusClineMatrix_p_Melt[,1])==length(locusClineMatrix_exp_Melt[,1])
  locusClineMatrix_All_Melt <- cbind(locusClineMatrix_p_melt,locusClineMatrix_p_exp_melt[,"p_exp",drop=F],
                                     locusClineMatrix_p_logit_melt[,"p_logit",drop=F],
                                     locusClineMatrix_exp_logit_melt[,"p_logit_exp",drop=F],
                                     locusClineMatrix_col_melt[,"col",drop=F])
  locusClineMatrix_All_Melt[,"Xbin"] <- locusClineMatrix_All_Melt[,"Xbin"]/1000
  plot.clines <- ggplot(locusClineMatrix_All_Melt, aes(Xbin, p)) + 
                geom_line(aes(x=Xbin, y=p_exp),colour="blue",size = 1.8,na.rm = TRUE) + # curve
                geom_point(colour="black",fill="gray",size=3,shape=21,alpha = 0.7,na.rm = TRUE) + # transparent points
                facet_wrap(~locus,ncol=numColFacet) + # multi-faceted with 5 per column
                #background_grid(major = 'xy', minor = "none",size.major = 0.6) + # add thin horizontal lines 
                #theme_light() +
                theme_bw() +
                #theme(panel.border = element_blank())+
                panel_border() + # and a border around each panel
                theme(axis.title.x = element_text(size=20,vjust=-0.5)) + 
                theme(axis.title.y = element_text(size=20,vjust=1.2)) + 
                theme(strip.text.x = element_text(size=15)) +
                theme(axis.text=element_text(size=18)) +
                ylab("allele frequency (p)")  + xlab("distance (km)") + theme(legend.position="none") 
  return(plot.clines)
}