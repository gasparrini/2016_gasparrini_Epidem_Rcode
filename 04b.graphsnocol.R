################################################################################
# Updated version of the code for the analysis in:
#
#   "Modelling lagged associations in environmental time series data:
#     a simulation study"
#   Antonio Gasparrini
#   Epidemiology 2016
#   http://www.ag-myresearch.com/2016_gasparrini_epidem.html
#
# Update: 15 January 2017
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   https://github.com/gasparrini/2016_gasparrini_Epidem_Rcode
################################################################################

################################################################################
# GRAPHS
################################################################################

################################################################################
# BASELINE TREND AND SEASONAL COMPONENT

pdf("figure1nocol.pdf",height=4,width=7)

par(mar=c(5,4,3,1),mgp=c(2.5,1,0))

plot(as.Date(data$date),exp(timeeff),type="l",xlab="Time",ylab="Daily deaths",
  ylim=c(100,160),lwd=1.5,frame.plot=T,col=1)

dev.off()

################################################################################
# GRAPHICAL REPRESENTATION OF THE SIMULATED EFFECT SURFACES

# ARGUMENTS FOR 3D PLOTS
arg3D <- list(y=0:20,ticktype="detailed",theta=230,ltheta=150,phi=40,lphi=30,
  ylab="Lag",zlab="RR",shade=0.75,r=sqrt(3),d=5,cex.axis=0.7,cex.lab=0.8,
  border=grey(0.3),col=grey(0.99))

pdf("figure2nocol.pdf",height=8,width=7)

layout(matrix(1:4,ncol=2,byrow=T))

for(j in seq(true)) {
  
  par(mar=c(2,1,2,1))
  x <- if(j%in%1:2) 0:50 else -20:35
  zlim <- if(j%in%1:2) c(0.99,1.06) else c(0.99,1.25)
  xlab <- if(j%in%1:2) "Ozone (PPM)" else "Temperature (C)"
  d3 <- do.call(persp,modifyList(arg3D,list(x=x,z=exp(true[[j]]),
    zlim=zlim,xlab=xlab)))
  title(names(combsim[j]),cex.main=1.2)
  if(j%in%1:2) {
    lines (trans3d(x=10,y=0:20,
      z=exp(true[[j]][as.character(10),paste0('lag',0:20)]),
      pmat=d3),col=1,lwd=2)
  } else {
    lines (trans3d(x=30,y=0:20,
      z=exp(true[[j]][as.character(30),paste0('lag',0:20)]),
      pmat=d3),col=1,lwd=2)
    lines (trans3d(x=-15,y=0:20,
      z=exp(true[[j]][as.character(-15),paste0('lag',0:20)]),
      pmat=d3),col=1,lwd=2)
  }
}

dev.off()

################################################################################
# GRAPH OF SIMULATION RESULTS

pdf("figure3nocol.pdf",height=3.5,width=9)

layout(matrix(1:4,2),heights=c(0.1,1))

for(j in 1:2) {
  
  par(mar=c(0,4,0,1))
  plot(1:10,type="n",frame.plot=F,axes=F,xlab="",ylab="")
  text(5,5,names(combsim)[j],cex=1.5,font=2)
  
  par(mar=c(4,4,2,1),mgp=c(2.5,1,0))
  
  plot(0:20,exp(true[[j]][as.character(10),]),type="n",xlab="Lag",
    ylab="RR",ylim=c(0.995,1.01),bty="l")
  for(m in sample[[j]][[1]]) lines(0:20,exp(m[as.character(10),]),
    lwd=1.5,col=grey(0.85))
  abline(h=1)
  lines(0:20,exp(true[[j]][as.character(10),]),lwd=1.5)
  lines(0:20,exp(pred[[j]][[1]][as.character(10),]),lty=5,lwd=1.5,col=1)
  mtext("Lag-response at 10 PPM",cex=0.8)
  legend("top",c("True","DLM (average)","DLM (sample)"),lty=c(1,5,1),
    lwd=1,col=c(1,1,grey(0.85)),cex=0.7,bty="n",inset=0.05)
}

dev.off()


pdf("figure4nocol.pdf",height=5,width=8)

layout(matrix(1:6,3),heights=c(0.1,1,1))

for(j in 3:4) {
  
  par(mar=c(0,4,0,1))
  plot(1:10,type="n",frame.plot=F,axes=F,xlab="",ylab="")
  text(5,5,names(combsim)[j],cex=1.5,font=2)
  
  par(mar=c(4,4,2,1),mgp=c(2.5,1,0))
  
  plot(0:20,exp(true[[j]][as.character(30),]),type="n",xlab="Lag",
    ylab="RR",ylim=c(0.98,1.10),bty="l")
  for(m in sample[[j]][[1]]) lines(0:20,exp(m[as.character(30),]),
    lwd=1.5,col=grey(0.85))
  abline(h=1)
  lines(0:20,exp(true[[j]][as.character(30),]),lwd=1.5)
  lines(0:20,exp(pred[[j]][[1]][as.character(30),]),lty=5,lwd=1.5,col=1)
  mtext("Lag-response at 30C",cex=0.8)
  legend("top",c("True","DLNM (average)","DLNM (sample)"),lty=c(1,5,1),
    lwd=1,col=c(1,1,grey(0.85)),cex=0.7,bty="n",inset=0.05)
  
  plot(0:20,exp(true[[j]][as.character(-15),]),type="n",xlab="Lag",
    ylab="RR",ylim=c(0.98,1.05),bty="l")
  for(m in sample[[j]][[1]]) lines(0:20,exp(m[as.character(-15),]),
    lwd=1.5,col=grey(0.85))
  abline(h=1)
  lines(0:20,exp(true[[j]][as.character(-15),]),lwd=1.5)
  lines(0:20,exp(pred[[j]][[1]][as.character(-15),]),lty=5,lwd=1.5,col=1)
  mtext("Lag-response at -15C",cex=0.8)
  legend("top",c("True","DLNM (average)","DLNM (sample)"),lty=c(1,5,1),
    lwd=1,col=c(1,1,grey(0.85)),cex=0.7,bty="n",inset=0.05)
}

dev.off()


pdf("figure5nocol.pdf",height=6.5,width=7)

layout(matrix(1:4,2,byrow=T))

for(j in seq(true)) {
  
  par(mar=c(4,4,4,1),mgp=c(2.5,1,0))
  
  x <- if(j%in%1:2) 0:50 else -20:35
  ylim <- if(j%in%1:2) c(0.95,1.25) else c(0.9,1.50)
  xlab <- if(j%in%1:2) "Ozone (PPM)" else 
    expression(paste("Temperature (",degree,"C)"))
  
  plot(x,exp(truecum[[j]]),type="n",xlab=xlab,ylab="RR",ylim=ylim,bty="l")
  abline(h=1)
  lines(x,exp(truecum[[j]]),lwd=1.5)
  lines(x,exp(cumpred[[j]][[1]]),lty=5,lwd=1.5,col=1)
  lines(x,exp(cumpred[[j]][[2]]),lty=2,lwd=1.5,col=1)
  lines(x,exp(cumpred[[j]][[4]]),lty=4,lwd=1.5,col=1)
  title(names(combsim)[j],cex=0.8)
  legend("top",c("True",mod[c(1,2,4)]),lty=c(1,5,2,4),lwd=1,col=1,cex=0.7,
    bty="n",inset=0.01,ncol=2)

}

dev.off()

#
