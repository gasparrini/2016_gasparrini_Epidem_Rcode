################################################################################
# SUPPLEMENTAL MATERIAL of the article:
#   "Modelling lagged associations in environmental time series data:
#     a simulation study"
#   Antonio Gasparrini
#   Epidemiology 2016
#
# This code reproduces the simulation results illustrated in the article
#
# 12 July 2016
# * an updated version of this code, (hopefully) compatible with future
#   versions of the software, is available at the personal website of the
#   first author (www.ag-myresearch.com)
################################################################################

################################################################################
# SUPPLEMENTARY MATERIAL
################################################################################

################################################################################
# MOVING AVERAGE MODEL WITH LAG 0-7

bmp("efigure1.bmp",width=1950,height=1800,pointsize=48)

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
  lines(x,exp(cumpred[[j]][[3]]),lty=5,lwd=1.5,col=2)
  title(names(combsim)[j],cex=0.8)
  legend("top",c("True",mod[3]),lty=c(1,5),lwd=1,col=1:2,cex=0.7,
    bty="n",inset=0.01,ncol=2)
  
}

dev.off()

################################################################################
# COMPARISON OF MOVING AVERAGE DEFINITIONS

bmp("efigure2.bmp",width=1950,height=1800,pointsize=48)

layout(matrix(1:4,2))

for(j in 3:4) {
  
  par(mar=c(4,4,4,1),mgp=c(2.5,1,0))
  
  xlab <- expression(paste("Temperature (",degree,"C)"))
  
  plot(-20:35,exp(truecum[[j]]),type="n",xlab=xlab,ylab="RR",
    ylim=c(0.9,1.50),bty="l")
  abline(h=1)
  lines(-20:35,exp(cumpred[[j]][[3]]),lwd=1.5)
  lines(xpred07b,exp(cumpredsens[[2*j+2]]),lty=5,lwd=1.5,col=2)
  title(names(combsim)[j],cex=0.8)
  legend("top",paste(mod[3],c("Eq. 3.1","Eq. 3.2")),lty=c(5,1),
    lwd=1,col=2:1,cex=0.7,bty="n",inset=0.05)
  
  plot(-20:35,exp(truecum[[j]]),type="n",xlab=xlab,ylab="RR",
    ylim=c(0.9,1.50),bty="l")
  abline(h=1)
  lines(-20:35,exp(cumpred[[j]][[4]]),lwd=1.5)
  lines(xpred020b,exp(cumpredsens[[2*j+3]]),lty=5,lwd=1.5,col=2)
  title(names(combsim)[j],cex=0.8)
  legend("top",paste(mod[4],c("Eq. 3.1","Eq. 3.2")),lty=c(5,1),
    lwd=1,col=2:1,cex=0.7,bty="n",inset=0.05)
}

dev.off()

################################################################################
# ASSUMPTIONS ON THE BI-DIMENSIONAL SHAPE

bmp("efigure3.bmp",width=1650,height=2750,pointsize=48)


layout(matrix(1:8,4),heights=c(0.1,1,1))

for(j in c(2,4)) {
  
  par(mar=c(0,4,0,1))
  plot(1:10,type="n",frame.plot=F,axes=F,xlab="",ylab="")
  text(5,5,names(combsim)[j],cex=1.3,font=2)
  
  par(mar=c(2,1,2,1))
  
  x <- if(j==2) 0:50 else -20:35
  zlim <- if(j==2) c(0.99,1.06) else c(0.99,1.25)
  xlab <- if(j==2) "Ozone (PPM)" else "Temperature (C)"
  
  do.call(persp,modifyList(arg3D,list(x=x,z=exp(true[[j]]),zlim=zlim,
    xlab=xlab,main="True")))
  do.call(persp,modifyList(arg3D,list(x=x,z=exp(pred[[j]][[1]]),zlim=zlim,
    xlab=xlab,main=ifelse(j==2,"DLM","DLNM"))))
  do.call(persp,modifyList(arg3D,list(x=x,z=exp(pred[[j]][[4]]),zlim=zlim,
    xlab=xlab,main="MA 0-20")))
  
}

dev.off()

################################################################################
# SENSITIVITY ANALYSES

bmp("efigure4.bmp",width=1950,height=1800,pointsize=48)


layout(matrix(1:4,2))

par(mar=c(4,4,4,1),mgp=c(2.5,1,0))


plot(0:50,exp(truecum[[1]]),type="n",xlab="Ozone (PPM)",ylab="RR",
  ylim=c(0.95,1.25),bty="l")
abline(h=1)
lines(0:50,exp(truecum[[1]]),lwd=1.5)
lines(0:50,exp(cumpredsens[[1]]),lty=5,lwd=1.5,col=2)
lines(0:50,exp(cumpredsens[[2]]),lty=4,lwd=1.5,col=4)
title(names(combsimsens)[1],cex=0.8)
legend("top",c("True","Unconstrained DLM, lag 0-3","DLM with 2 knots, lag 0-7"),
  lty=c(1,5,4),lwd=1,col=c(1,2,4),cex=0.7,bty="n",inset=0.01,ncol=1)

xlab <- expression(paste("Temperature (",degree,"C)"))

plot(-20:35,exp(truecum[[3]]),type="n",xlab=xlab,ylab="RR",
  ylim=c(0.9,1.50),bty="l")
abline(h=1)
lines(-20:35,exp(truecum[[3]]),lwd=1.5)
lines(-20:35,exp(cumpredsens[[3]]),lty=5,lwd=1.5,col=2)
lines(-20:35,exp(cumpredsens[[4]]),lty=4,lwd=1.5,col=4)
title(names(combsimsens)[3],cex=0.8)
legend("top",c("True","Unconstrained DLNM, lag 0-3","DLNM with 2 knots, lag 0-7"),
  lty=c(1,5,4),lwd=1,col=c(1,2,4),cex=0.7,bty="n",inset=0.01,ncol=1)

plot(-20:35,exp(truecum[[4]]),type="n",xlab=xlab,ylab="RR",
  ylim=c(0.9,1.50),bty="l")
abline(h=1)
lines(-20:35,exp(truecum[[4]]),lwd=1.5)
lines(-20:35,exp(cumpredsens[[5]]),lty=5,lwd=1.5,col=2)
lines(-20:35,exp(cumpredsens[[6]]),lty=4,lwd=1.5,col=4)
title(names(combsimsens)[3],cex=0.8)
legend("top",c("True","DLNM, w(l) with 2 knots","DLNM, w(l) with 4 knots"),
  lty=c(1,5,4),lwd=1,col=c(1,2,4),cex=0.7,bty="n",inset=0.01,ncol=1)

plot(-20:35,exp(truecum[[4]]),type="n",xlab=xlab,ylab="RR",
  ylim=c(0.9,1.50),bty="l")
abline(h=1)
lines(-20:35,exp(truecum[[4]]),lwd=1.5)
lines(-20:35,exp(cumpredsens[[7]]),lty=5,lwd=1.5,col=2)
title(names(combsimsens)[3],cex=0.8)
legend("top",c("True with overdispersion","DLNM with quasi-Poisson family"),
  lty=c(1,5),lwd=1,col=c(1,2),cex=0.7,bty="n",inset=0.01,ncol=1)

dev.off()

#
