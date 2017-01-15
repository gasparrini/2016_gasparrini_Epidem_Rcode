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
# SET THE PARAMETERS AND DEFINE OBJECTS FOR THE SIMULATIONS
################################################################################

# LOAD PACKAGES
library(dlnm) ; library(splines) ; library(tsModel) ; library(MASS)

################################################################################
# DEFINE THE EXPOSURE AND THE CONFOUNDER MODEL

# REAL DATA
data <- chicagoNMMAPS

# REAL OZONE SERIES FOR NON-LINEAR SCENARIO, STANDARDIZED IN 0 - 50 PPB
ozone <- ((data$o3-min(data$o3))/diff(range(data$o3)))*50

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
Qozone <- Lag(ozone,0:20)

################################################################################
# DEFINE THE BI-DIMENSIONAL ASSOCIATIONS USED FOR SIMULATING DATA

# BASIC FUNCTIONS TO SIMULATE UNIDIMENSIONAL SHAPES
flin <- function(x) 0.001*x
fcurv <- function(x) ifelse(x>20,((x-20)/2.4)^2.5/450,((20-x)/200)^1.5)
wdecay <- function(lag,scale) exp(-lag/scale)
wpeak1 <- function(lag) 5*dnorm(lag,2,5)
wpeak2 <- function(lag) 4*dchisq(lag,6)

# FUNCTIONS TO SIMULATE THE BI-DIMENSIONAL EXPOSURE-LAG-RESPONSE
flinshort <- function(x,lag) flin(x) * wdecay(lag,0.7)
flinlong <- function(x,lag) flin(x) * wpeak1(lag)
fcurvshort <- function(x,lag) fcurv(x) * wdecay(lag,0.7)
fcurvlong <- function(x,lag) fcurv(x) * ifelse(x>20,wdecay(lag,0.7),wpeak2(lag))

# COMBINATIONS OF FUNCTIONS USED TO SIMULATE DATA
combsim <- c("flinshort","flinlong","fcurvshort","fcurvlong")
names(combsim) <- c("Linear short-lag scenario","Linear long-lag scenario",
  "Non-linear short-lag scenario","Non-linear long-lag scenario")

# LIST WITH TRUE BI-DIMENSIONAL EXPOSURE-LAG-RESPONSE FOR EACH COMBINATION
true <- lapply(seq(combsim), function(j) {
  xpred <- if(j%in%1:2) 0:50 else -20:35
  temp <- outer(xpred,0:20,combsim[j])
  dimnames(temp) <- list(xpred,paste("lag",0:20,sep=""))
  return(temp)
})
names(true) <- combsim

# LIST OF TRUE OVERALL CUMULATIVE EXPOSURE-RESPONSE FOR EACH COMBINATION
truecum <- lapply(true,rowSums)
names(truecum) <- combsim

# FUNCTION TO COMPUTE THE CUMULATIVE EFFECT GIVEN AN EXPOSURE HISTORY
fcumeff <- function(hist,lag,fun) sum(do.call(fun,list(hist,lag)))

################################################################################
# SETTINGS FOR THE SIMULATION

# NUMBER OF ITERATIONS
# NB: CHANGE TO 5000 TO REPRODUCE THE RESULTS (BUT LONG COMPUTING TIME)
nsim <- 3

# NUMBER OF SAMPLES USED AS EXAMPLES OF INDIVIDUAL ESTIMATES
nsample <- min(nsim,20)

# NOMINAL VALUE
qn <- qnorm(0.975)

#
