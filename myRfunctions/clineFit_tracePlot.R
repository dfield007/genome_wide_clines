# trace plot of SANN results. Panels include:
# 1. -lnL profile for accepted values across iterations
# 2. step size for centre
# 3. centre parameter across iterations
# 4. width parameter across iterations
# 5. p0 parameter across iterations
# 6. p1 parameter across iterations
# 7. Fst parameter across iterations
# 8. Overall parameter values accepts for width and centre
# 9. Overall parameter values accepts for width and p0

tracePlot <- function(clineFit,locusName,acceptRate) {
  par(mfrow=c(3,3),mar=c(3.9, 3.9, 2.7, 1),oma=c(1,1,1,1))  
  plot((1:nrow(clineFit$retvals)),clineFit$retvals[,"val"], main=paste(locusName,acceptRate,sep=", Accept rate="),type="l",xlab="iteration",ylab="lnL")
  mtext(paste("max -lnL = ",round(clineFit$minimum,3)),3,cex=0.8)
  valsSubHead <- c(paste("c=",round(clineFit$estimate[1],0),sep=""),paste("w=",round(clineFit$estimate[2],0),sep=""),
                   paste("p0=",round(clineFit$estimate[3],2),sep=""), paste("p1=",round(clineFit$estimate[4],2),sep=""),
                   paste("Fst=",round((1/(clineFit$estimate[5]+1)),2),sep=""))
  plot(((1:nrow(clineFit$retvals))),clineFit$retvals[,"deltap1"], main="Step size - centre",type="l",xlab="iteration",ylab="delta_centre")
  mtext(paste(valsSubHead,collapse=", "),cex=0.8)
  # plot(((1:nrow(clineFit$retvals))),clineFit$retvals[,"deltap2"], main="Step size - width",type="l",xlab="iteration",ylab="delta_width")
  # mtext(paste(valsSubHead,collapse=", "),cex=0.8)
  plot((1:nrow(clineFit$retvals)),clineFit$retvals[,"p1"], main="Values centre",type="l",xlab="iteration",ylab="centre",col="black",lwd=0.2)
  plot(((1:nrow(clineFit$retvals))),clineFit$retvals[,"p2"], main="Values width",type="l",xlab="iteration",ylab="width",col="black",lwd=0.2)
  plot(((1:nrow(clineFit$retvals))),clineFit$retvals[,"p3"], main="Values p0",type="l",xlab="iteration",ylab="p0",col="black",lwd=0.2)
  plot(((1:nrow(clineFit$retvals))),clineFit$retvals[,"p4"], main="Values p1",type="l",xlab="iteration",ylab="p1",col="black",lwd=0.2)
  plot(((1:nrow(clineFit$retvals))*10),(1/(clineFit$retvals[,"p5"]+1)), main="Values Fst",type="l",xlab="iteration",ylab="Fst",col="black",lwd=0.2)
  #plot(clineFit$retvals[,"p1"],clineFit$retvals[,"p2"], main="centre / width",xlab="centre",ylab="width",col="black")
  #plot(clineFit$retvals[,"p3"],clineFit$retvals[,"p2"], main="p0 / width",xlab="p0",ylab="width",col="black")
  smoothScatter(clineFit$retvals[,"p1"],clineFit$retvals[,"p2"], main="centre / width",xlab="centre",ylab="width")
  smoothScatter(clineFit$retvals[,"p3"],clineFit$retvals[,"p2"],main="p0 / width",xlab="p0",ylab="width")
  
}      
