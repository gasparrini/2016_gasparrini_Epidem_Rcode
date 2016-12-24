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
# EXAMPLE OF SIMULATION OF THE DATA AND FIT OF THE MODELS
################################################################################

################################################################################
# DEFINE THE EXPOSURE AND THE CONFOUNDER MODEL

# LOAD PACKAGES
library(dlnm) ; library(splines) ; library(tsModel) ; library(MASS)

# REAL DATA
data <- chicagoNMMAPS

# REAL TEMPERATURE SERIES FOR NON-LINEAR SCENARIO, STANDARDIZED IN -20 TO 35C
tmean <- ((data$temp-min(data$temp))/diff(range(data$temp)))*55-20

# REAL MORTALITY SERIES
y <- data$death

# SPLINES OF TIME
timebasis <- ns(data$time,df=10*14)

# SIMULATE BASELINE MORTALITY
m <- glm(y~timebasis,family=poisson,na.action=na.exclude)
timeeff <- predict(m,newdata=list(timebasis=timebasis))

# MATRIX OF EXPOSURE HISTORIES
Qtmean <- Lag(tmean,0:20)

################################################################################
# DEFINE THE BI-DIMENSIONAL ASSOCIATION AND SIMULATE THE OUTCOME
# EXAMPLE RELATED TO THE NON-LINEAR LONG-LAG SCENARIO
# (SEE THE SCRIPTS REPRODUCING THE FULL ANALYSIS FOR OTHER FUNCTIONS/SCENARIOS)

# BASIC FUNCTIONS TO SIMULATE UNIDIMENSIONAL SHAPES
fcurv <- function(x) ifelse(x>20,((x-20)/2.4)^2.5/450,((20-x)/200)^1.5)
wdecay <- function(lag,scale) exp(-lag/scale)
wpeak2 <- function(lag) 4*dchisq(lag,6)

# FUNCTION TO SIMULATE THE BI-DIMENSIONAL EXPOSURE-LAG-RESPONSE
fcurvlong <- function(x,lag) fcurv(x) * ifelse(x>20,wdecay(lag,0.7),wpeak2(lag))

# TRUE BI-DIMENSIONAL EXPOSURE-LAG-RESPONSE ASSOCIATION 
true <- outer(-20:35,0:20,fcurvlong)
dimnames(true) <- list(-20:35,paste("lag",0:20,sep=""))

# TRUE OVERALL CUMULATIVE EXPOSURE-RESPONSE ASSOCIATION
truecum <- rowSums(true)

# DEFINE THE TRUE EFFECT
fcumeff <- function(hist,lag,fun) sum(do.call(fun,list(hist,lag)))
cumeff <- apply(Qtmean,1,fcumeff,0:20,fcurvlong)

# SET THE SEED AND SIMULATE THE OVERDISPERSED RESPONSE (WARNING DUE TO MISSING)
# DEFINE THE OVERDISPERSION PHI
phi <- 1.3
# DEFINE THETA (PHI = 1 + MU/THETA, SO THETA = MU/(PHI-1) )
theta <- exp(timeeff+cumeff)/(phi-1)
# SIMULATE FROM A NEGATIVE BINOMIAL DISTRIBUTION
set.seed(1)
y <- rnegbin(length(cumeff),exp(timeeff+cumeff),theta)

# TO SIMULATE FROM A POISSON DISTRIBUTION WITH NO OVERDISPERSION
#y <- rpois(length(cumeff),exp(timeeff+cumeff))

################################################################################
# RUN THE MODELS

# DLNM
vk <- quantile(tmean,c(10,75,90)/100)
cb1 <- crossbasis(tmean,lag=20,argvar=list("bs",degree=2,knots=vk),
  arglag=list(knots=logknots(20,nk=3)))
m1 <- glm(y~cb1+timebasis,family=quasipoisson)
cp1 <- crosspred(cb1,m1,at=-20:35,cen=20)

# MOVING AVERAGE (MA) MODEL WITH LAG 0-20
cb2 <- crossbasis(tmean,lag=20,argvar=list("bs",degree=2,knots=vk),
  arglag=list("strata",df=1))
m2 <- glm(y~cb2+timebasis,family=quasipoisson)
cp2 <- crosspred(cb2,m2,at=-20:35,cen=20)

# MOVING AVERAGE (MA) MODEL WITH LAG 0-20 (STANDARD DEFINITION)
tmean020 <- rowMeans(Lag(tmean,0:20))
cb3 <- onebasis(tmean020,"bs",degree=2,knots=vk)
m3 <- glm(y~cb3+timebasis,family=quasipoisson)
cp3 <- crosspred(cb3,m3,cen=20)

# COMPARISON OF TRUE AND ESTIMATED EXPOSURE-LAG-RESPONSE SURFACES
persp(x=-20:35,y=0:20,z=exp(true),ylab="Lag",zlab="RR",zlim=c(0.98,1.25),
  ticktype="detailed",theta=230,ltheta=150,phi=40,shade=0.75,d=5,
  col="lightskyblue",main="True")
plot(cp1,zlim=c(0.98,1.25),theta=230,ltheta=150,phi=40,main="DLNM")
plot(cp2,zlim=c(0.98,1.25),theta=230,ltheta=150,phi=40,main="MA 0-20")

# COMPARISON OF TRUE AND ESTIMATED OVERALL CUMULATIVE EXPOSURE-RESPONSE CURVE
plot(-20:35,exp(truecum),type="l",xlab="Temperature (C)",ylab="RR",
  ylim=c(0.9,1.50),bty="l",main="Comparison")
abline(h=1)
lines(cp1,"overall",col=2,lty=5)
lines(cp2,"overall",col=4,lty=2)
lines(cp3,"overall",col=3,lty=4)
legend("top",c("True","DLNM","MA0-20","MA0-20 (std)"),col=c(1,2,4,3),
  lty=c(1,5,2,4),bty="n")

#
