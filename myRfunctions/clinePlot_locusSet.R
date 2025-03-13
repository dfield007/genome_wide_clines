# Cline plots locus set (with and without cline fitting)
# For repeating figure generations over a vector of loci names
# David Field 16/6/2017
clinePlot_locusSet <- function(clineModel,locusSet,clineFitsSummary,thisTrail,thisPop,plotTypes,refLocus=NA,refPoints=F) {
  # HZ_data <- HZ_Planoles
  # clineFitsSummary <- ros_assembly_rep3
  # locusSet=as.vector(clineFitsSummary[,"locusName"])
  # thisTrail <- "test"; thisPop <- "Planoles"; plotTypes <- c("pdf")
  # refLocus="ros_assembly_541834"; refPoints=T
  cat("\n Generating cline plots.... ") 
  cat("\n refLocus= ",refLocus) 
  cat("\n Cline fits present for loci: ") 
  cat("... ",as.vector(clineFitsSummary[which(clineFitsSummary[,"centre"]!=0),"locusName"])) 
  cat("\n  ") 
  cat("\n Only allele frequencies (no cline fits) for loci: ") 
  cat("... ",as.vector(clineFitsSummary[(clineFitsSummary[,"centre"]==0 & as.numeric(as.vector(clineFitsSummary[,"N_total"]))>0),"locusName"])) 
  for (thisLocus in 1:length(locusSet)) {
    # thisLocus <- 1
    locusName <- locusSet[thisLocus]
    thisDrop <- which(clineFitsSummary[,"locusName"]==locusName)
    #clineFitsSummary[thisDrop,]
    if (clineFitsSummary[thisDrop,"N_total"]==0 | clineFitsSummary[thisDrop,"N_total"]=="noData") {
      next
    }
    if (clineFitsSummary[thisDrop,"centre"]==0 & as.numeric(as.vector(clineFitsSummary[thisDrop,"N_total"]))!=0) {
      # plot cline
      
   #   clinePlot(clineModel="no fit",locusName,clineFitsSummary,curveFit=F,cexMain=1.2,cexAxisPlot= 1.85,
  #              axisLineWdth=2.5,cexPoints=1.8,lwdCurve=4,refLocus,refPoints)      
      if ("pdf"%in%plotTypes) {
        # as pdf
        fileClineFitName <- paste(paste(paste(locusName,"ClineFit",sep="_"),thisTrail,sep="_"),".pdf",sep="")
        pdf(fileClineFitName,width = 15, height = 15,pointsize=12)
        clinePlot(clineModel="no fit",locusName,clineFitsSummary,curveFit=F,cexMain=1.2,cexAxisPlot= 2.7,
                  axisLineWdth=3.8,cexPoints=3,lwdCurve=4.3,refLocus,refPoints)
        dev.off()
      }
      if ("png"%in%plotTypes) {
        # as png
        fileClineFitName <- paste(paste(paste(locusName,"ClineFit",sep="_"),thisTrail,sep="_"),".png",sep="")
        png(fileClineFitName,width = 1000, height = 1000, pointsize = 14,antialias="default")
        clinePlot(clineModel="no fit",locusName,clineFitsSummary,curveFit=F,cexMain=1.4,cexAxisPlot= 2,
                  axisLineWdth=2,cexPoints=1.8,lwdCurve=4,refLocus,refPoints)
        dev.off()
      }
    }
    if (as.numeric(as.vector(clineFitsSummary[thisDrop,"centre"]!=0))) {
      # plot cline
     # clinePlot(clineModel="symm_5p",locusName,clineFitsSummary,curveFit=T,cexMain=1.2,cexAxisPlot= 1.85,
    #            axisLineWdth=2.5,cexPoints=1.8,lwdCurve=4,refLocus,refPoints)      
      
      if ("pdf"%in%plotTypes) {
        # as pdf
        fileClineFitName <- paste(paste(paste(locusName,"ClineFit",sep="_"),thisTrail,sep="_"),".pdf",sep="")
        pdf(fileClineFitName,width = 15, height = 15,pointsize=12)
        clinePlot(clineModel,locusName,clineFitsSummary,curveFit=T,cexMain=1.2,cexAxisPlot= 2.7,
                  axisLineWdth=3.8,cexPoints=3,lwdCurve=4.3,refLocus,refPoints)
        dev.off()
      }
      if ("png"%in%plotTypes) {
        # as png
        fileClineFitName <- paste(paste(paste(locusName,"ClineFit",sep="_"),thisTrail,sep="_"),".png",sep="")
        png(fileClineFitName,width = 1000, height = 1000, pointsize = 14,antialias="default")
        clinePlot(clineModel,locusName,clineFitsSummary,curveFit=T,cexMain=1.4,cexAxisPlot= 2,
                  axisLineWdth=2,cexPoints=1.8,lwdCurve=4,refLocus,refPoints)
        dev.off()
      }
    }
  }
}
