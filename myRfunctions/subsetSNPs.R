# subsets the SNP data. provided a main list of loci to select, this function can also randomly sample down (randomNum) any number of 
# SNP loci. IF the randomNum is equal or more that the list, it will simply use all of the loci available.

subsetSNPs <- function(thisData,locusList,randomNum,verbose) {
  # thisData=GenoEcol_Planoles; locusList = loci_95; randomNum = 95; locusList=loci_105
  pedLoci <- which(colnames(thisData)%in%locusList)
  # check if any of loci in 2009-2014 set are not in the pedigreeLoci_strict
  #SNPlociListUpdated[(SNPlociListUpdated[,"genos_2009_2014_only"]==1 |  SNPlociListUpdated[,"genos_allYears"]==1),"LocusName"]
  
  # subset markers to pedigree set
  thisData_cut <- subset(thisData, select=c(PlantID_final:GenoScore,pedLoci,DataSet:propMissing))
  #GenoEcol_Planoles_ped <- GenoEcol_Planoles[,c(1:13,pedLoci,236:ncol(GenoEcol_Planoles))]
  
  # remove cpDNA markers
  dropCP <- c("scpDNA_seq_21466","scpDNA_seq_58467")
  thisData_cut <- thisData_cut[,(colnames(thisData_cut)%in%dropCP)==F]
  
  # re-establish which columns have locus data
  lociColumns <- (which(colnames(thisData_cut)=="GenoScore")+1):(which(colnames(thisData_cut)=="DataSet")-1)
  lociNames <- colnames(thisData_cut)[lociColumns]
  #lociColumns_No_s91 <- lociColumns[-which(lociNames=="s91_78256")] # exclude s91 from missing data calcs
  
  # another set excluding clinal markers in LD with each other
  dropClines <- c("s1187_290152","ros_assembly_543443", "ros_assembly_567004", "ros_assembly_575837","ros_assembly_576271", "ros_assembly_620992", "ros_assembly_653015", "ros_assembly_670530", "ros_assembly_715001","s829_8371", "s829_204463","s316_93292","s316_257789","s91_78256")
  
  lociNames_noClines <- lociNames[lociNames%in%dropClines==F]
  lociColumns_noClines <- which(colnames(thisData_cut)%in%lociNames_noClines)
  
  # randomly select loci
  if (randomNum<length(locusList)) {
    theseLoci <- sample(lociNames,randomNum,replace=F)
    pedLoci_cols <- which(colnames(thisData_cut)%in%theseLoci)
    # subset markers to pedigree set
    GenoEcol_Planoles_ped <- cbind(thisData_cut[,which(colnames(thisData_cut)=="PlantID_final"):which(colnames(thisData_cut)=="GenoScore")],thisData_cut[,pedLoci_cols],thisData_cut[,which(colnames(thisData_cut)=="DataSet"):which(colnames(thisData_cut)=="propMissing")])
  }
  # recalculate totalLoci, missing etc
  if (randomNum>=length(locusList)) {
    theseLoci <- lociNames
    pedLoci_cols <- which(colnames(thisData_cut)%in%theseLoci)
    # subset markers to pedigree set
    GenoEcol_Planoles_ped <- cbind(thisData_cut[,which(colnames(thisData_cut)=="PlantID_final"):which(colnames(thisData_cut)=="GenoScore")],thisData_cut[,pedLoci_cols],thisData_cut[,which(colnames(thisData_cut)=="DataSet"):which(colnames(thisData_cut)=="propMissing")])
  }
  
  # totals_rawSNPs <- as.vector(apply(thisData_cut[,pedLoci_cols],1, FUN = getTotalGenotyped))
  # missing <- as.vector(apply(thisData_cut[,pedLoci_cols],1, FUN = getMissing))
  # NGY <- as.vector(apply(thisData_cut[,pedLoci_cols],1, FUN = getNGY))
  # totals <- as.vector(apply(thisData_cut[,pedLoci_cols],1, FUN = getTotalGenotyped))
  # totals_Available <- as.vector(apply(thisData_cut[,pedLoci_cols],1, FUN = getTotalGenotypesAvailable))
  # 
  # propMissing <- missing/totals_Available
  # 
  # # check problem with 2015 & 2016 plants -10 issue
  # #colnames(GenoEcol_Planoles_ped)
  # #GenoEcol_Planoles_ped[GenoEcol_Planoles_ped[,"FirstYear"]==2015,"s899_432971"]
  # 
  # # updating
  # GenoEcol_Planoles_ped[,"missing"] <- missing
  # GenoEcol_Planoles_ped[,"propMissing"] <- propMissing
  # GenoEcol_Planoles_ped[,"totals"] <- totals
  # GenoEcol_Planoles_ped[,"NGY"] <- NGY
  # GenoEcol_Planoles_ped[,"totals_available"] <- totals_Available
  if (verbose==T) {
    cat("\n Loci remaining: ",length(pedLoci_cols))
  }
  return(GenoEcol_Planoles_ped)
}