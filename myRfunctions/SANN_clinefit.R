# Simulated Annealing (random walk) algorithm for obtaining estimates of Likelihood surface
# David Field (16/6/2017)
SANN_clinefit <- function (fn, positions, obsVals,initialParam, deltap = NULL, scale = 1, rptfreq = -1,
                           acceptscale = 1.01, rejectscale = 0.99, nmax = 10000, retvals = FALSE, 
                           retfreq = 10, verbose = T, ...) {
  # fn <- logLik_5par_SANN; startParams <- c(13000,2700,0.08,0.9,15)
  # initialParam <- startParams
  # deltap = NULL ; scale = 0.4 ; rptfreq = -1 ; acceptscale = 1.01 ; rejectscale = 0.99 ; nmax = 1000 ; 
  # retvals = TRUE;  retfreq = 100 ; verbose = TRUE
  # thisSummary <- summaryModelsFinal 
  # thisData <- dat2;  
  # positions <- XCordinates
  # logLik5_mcmc_4par_U_cwpLpR_test(c(13000,2700, 0.06,0.9),dat2,summaryModelsFinal)
  positions <- XCordinates
  ndim <- length(initialParam)
  paramSet <- initialParam
  minp <- paramSet
  # val <- fn(c(13000,2700, 0.06,0.9),thisData,thisSummary)
  val <- fn(initialParam,Xvals=positions,obsVals)
  minval <- val
  if (retvals) {
    info <- matrix(nrow = round(nmax/retfreq), ncol = 3 * ndim + 3)
    dimnames(info) <- list(NULL, c(paste("p", 1:ndim, sep = ""), 
                        paste("minp", 1:ndim, sep = ""), paste("deltap", 
                        1:ndim, sep = ""), "val", "minval", "accept"))
  }
  it <- 1
  if (is.null(deltap)) {
    deltap <- initialParam * 0.05
  }
  deltap[4] <- deltap[3]
  while (it <= nmax) {
    oldp <- paramSet
    oldval <- val
    pPotential <- paramSet + runif(ndim, -1, 1) * deltap
    while (pPotential[1]<500 | pPotential[1]>20000 | pPotential[2]<500 | pPotential[2]>20000 | pPotential[3]<0 | pPotential[3]>0.8 | pPotential[4]<0.3 | pPotential[4]>1.0 | pPotential[5]<1 | pPotential[5]>50) {
      pPotential <- paramSet + runif(ndim, -1, 1) * deltap
    }
    paramSet <- pPotential
    val <- fn(paramSet,Xvals=positions,obsVals)
    dval <- val - oldval
    saveinfo <- (retvals && it%%retfreq == 0)
    savecount <- it%/%retfreq
    if (saveinfo) {
      info[savecount, -ncol(info)] <- c(paramSet, minp, deltap, val, minval)
    }
    if (verbose) 
      if (it%%100 == 0) 
        cat(".")
    if (it%%500 == 0) 
      cat(it)
    # modified dval to >0, otherwise its finding the valleys NOT the peaks! 
    # Also changed the acceptance ratio to (L'/L'')*scale
    # new values always accepted if lnL'' > lnL'.
    # if lnL'' < lnL', accept with probability (lnL'/lnL'')*scale
    # scale seems to be required to reduce the probability, otherwise accept rate >> 0.5 
    if (dval > 0 || (((oldval/val)*scale) > runif(1))) {
      if (saveinfo) 
        info[savecount, ncol(info)] <- 1
      if (val > minval) {
        minval <- val
        minp <- paramSet
      }
      deltap <- deltap * acceptscale
    }
    else {
      if (saveinfo) 
        info[savecount, ncol(info)] <- 0
      val <- oldval
      paramSet <- oldp
      deltap <- deltap * rejectscale
    }
    it <- it + 1
    if ((rptfreq > 0) && (it%%rptfreq) == 1) 
      cat("it=", it, " value=", val, " min. val.=", minval, 
          "\n")
  }
  if (retvals) 
    return(list(minimum = minval, estimate = minp, funcalls = it, 
                retvals = info))
    else return(list(minimum = minval, estimate = minp, funcalls = it))
}
