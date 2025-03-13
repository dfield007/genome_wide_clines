# cline model - symetric sigmoid 4 parameters. 
# centre, width, p0 (pL), p1 (pR)
cline_model_sig_symm <- function(centre, width, pL, pR,Xvals) {
  ((pL + ((pR-pL)/(1 + exp((-4*(Xvals-centre))/width)))))
}

# cline model - 3 part cline with xx parameters
cline_model_sig_symm <- function(centre, width, pL, pR,Xvals) {
  
  LeftSide <- 
  
  ((pL + ((pR-pL)/(1 + exp((-4*(Xvals-centre))/width)))))
  
  
}