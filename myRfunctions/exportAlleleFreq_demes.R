# export allele frequency data from demes for set of loci
exportAlleleFreqDemes <- function(demeData) {
  # demeData <- HZ_200mDemes
  # minDemeSize <- 10
  # Find all X coordinates across all loci 
  # x coordinates available at each locus may differ slightly from the full set
  # i.e. some loci missing more data and thus deme number may differ
  XcordinateCollect <- NULL
  for (thisLocus in 1:length(names(demeData))) {
    # thisLocus <- 1
    thisLocusName <- names(demeData)[thisLocus]
    thisData <- demeData[[thisLocus]]
    # all X coorindates
    XCordinates <- thisData[["DemeMatrix2D_DistanceRefLine"]][,2]
    cat("\n ",thisLocusName)
    cat(".. num demes =", length(XCordinates))
    # x coordinates available at each locus may differ slightly from the full set
    # i.e. some loci missing more data and thus deme number may differ
    XcordinateCollect <- c(XcordinateCollect,round(XCordinates,1))
  }
  allXordinates <- sort(unique(XcordinateCollect))
  colnames_alleleFreq <- c("distAlongTransect",paste(names(demeData), rep(c("p"),length(names(demeData))),sep="_"),
                           paste(names(demeData), rep(c("N"),length(names(demeData))),sep="_"))
  alleleFreq <- matrix(0,nrow=length(allXordinates),length(colnames_alleleFreq))
  alleleFreq[,2:3] <- NA
  colnames(alleleFreq) <- colnames_alleleFreq
  alleleFreq[,"distAlongTransect"] <- allXordinates
  
  for (thisLocus in 1:length(names(demeData))) {
    # thisLocus <- 1
    thisLocusName <-names(demeData)[thisLocus]
    thisData <- demeData[[thisLocus]]
    #names(thisData)
    thisData_details <- thisData[["DemeMatrix2D_combine_noZero"]]
    for (thisDeme in 1:nrow(thisData_details)) {
      # thisDeme <- 1
      thisDemeDat <- thisData_details[thisDeme,]
      inds <- as.vector(apply(thisDemeDat[,6:7],1,sum)/2)
      putRow <- which(alleleFreq[,"distAlongTransect"]==thisDemeDat[,"distAlongTransect"])
      alleleFreq[putRow,which(colnames(alleleFreq)==paste(thisLocusName,"p",sep="_"))] <- thisDemeDat[,"p"]
      alleleFreq[putRow,which(colnames(alleleFreq)==paste(thisLocusName,"N",sep="_"))] <- inds
    }
  }
  return(alleleFreq)
}
