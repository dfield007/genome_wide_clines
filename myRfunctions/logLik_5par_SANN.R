# likelihood function for ML cline fitting with SANN
# David Field (16/6/2017)
logLik_5par_SANN <- function(theseParams,Xvals,obsVals) {
  #   theseParams <- c(13000,2700,0.08,0.9,15)
    pCount <- obsVals[,1]
    qCount <- obsVals[,2]
    Ngenes <- obsVals[,3]
  #   Xvals <- XCordinates
  #   Xvals <- XCordinates_good 
  #   pCount <- thisLoc_DemeMatrix2D_pqCountNoZero_good$pCount
  #   qCount <- thisLoc_DemeMatrix2D_pqCountNoZero_good$qCount
  #   Ngenes <- pCount+qCount
  centre <- theseParams[1]
  width <-  theseParams[2]
  pL <- theseParams[3]
  pR <- theseParams[4]
  alphaTerm <- theseParams[5]
  #   if (width<=800) {
  #     width <- 800
  #   }
  #   if (pL<0) {
  #     pL <- 0
  #   }
  #   if (pR>1) {
  #     pR <- 1
  #   }
  
  # pochMpfr(100, 4)
  # choose(100,4)
  # chooseZ(100, 4)
  # factorialZ(100,4)
  # log(pochMpfr((alphaTerm*cline_model_sig_symm(centre, width, pL, pR,Xvals))[1:5],pCount[1:5]))
  # log first?
  # 
  # chooseZ((alphaTerm*cline_model_sig_symm(centre, width, pL, pR,Xvals)),pCount)
  # test2 <- mpfr((alphaTerm*cline_model_sig_symm(centre, width, pL, pR,Xvals)),pCount)
  # 
  
  #part1 <- log(choose(Ngenes,pCount))
  #part2 <- log(pochMpfr((alphaTerm*cline_model_sig_symm(centre, width, pL, pR,Xvals)),pCount))
  
  #log(choose(Ngenes,pCount)) == choose(log(Ngenes),log(pCount))
  log_likelihood <- as.numeric(sum(log(choose(Ngenes,pCount)) +
                          log(pochMpfr((alphaTerm*cline_model_sig_symm(centre, width, pL, pR,Xvals)),pCount)) +
                          log(pochMpfr((alphaTerm*(1-cline_model_sig_symm(centre, width, pL, pR,Xvals))),qCount)) -
                          log(pochMpfr(alphaTerm,Ngenes))))
  return(log_likelihood)
}
