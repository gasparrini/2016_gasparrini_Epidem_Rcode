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
# RUN THE SENSITIVITY ANALYSIS
################################################################################

################################################################################
# CREATE SCENARIOS AND MODELS

sens <- paste0("sens",seq(11))

# SIMULATION SCENARIOS
combsimsens <- combsim[c(1,1,3,3,4,4,4,3,3,4,4)]

# LAG
laglist <- c(3,7,3,7,20,20,20,3,3,20,20)

# CROSS-BASIS ARGUMENTS (ONLY FOR DL MODELS)
varknots <- quantile(tmean,c(10,75,90)/100,na.rm=T)
argvarlist <- rep(list(list("lin"),list("bs",degree=2,knots=varknots)),c(2,5))
arglaglist <- list(
  list("integer"),
  list(knots=logknots(7,nk=2)),
  list("integer"),
  list(knots=logknots(7,nk=2)),
  list(knots=logknots(20,nk=2)),
  list(knots=logknots(20,nk=4)),
  list(knots=logknots(20,nk=3))
)
  
# STANDARD MOVING AVERAGE WITH RESTRICTED RANGE
ma07b <- rowMeans(Lag(tmean,0:7))
b07b <- onebasis(ma07b,"bs",degree=2,knots=varknots)
xpred07b <- seq(ceiling(attr(b07b,"range")[1]),floor(attr(b07b,"range")[2]))

ma020b <- rowMeans(Lag(tmean,0:20))
b020b <- onebasis(ma020b,"bs",degree=2,knots=varknots)
xpred020b <- seq(ceiling(attr(b020b,"range")[1]),floor(attr(b020b,"range")[2]))

################################################################################
# CREATE THE OBJECT TO STORE RESULTS

# LISTS FOR STORING CUMULATIVE PREDICTION, BIAS, COVERAGE, RMSE AND SAMPLES
# ALL THE MODELS
cumpredsens <- lapply(seq(combsimsens), function(i) 0)
names(cumpredsens) <- sens
cumbiassens <- cumcovsens <- cumrmsesens <- cumpredsens

################################################################################
# START THE ITERATIONS
# 2 LEVELS:
#   - TYPE OF SURFACE (COMBINATIONS OF FUNCTIONS USED FOR SIMULATE EFFECTS)
#   - RANDOM GENERATED RESPONSE

# LOOP ACROSS SIMULATED EXPOSURE-LAG-RESPONSES
for(j in seq(sens)) {
  
  # PRINT
  cat("\n\n",sens[j],"\n")
  
  # DEFINE THE TRUE EFFECT (DEPENDENT ON SCENARIO)
  cumeff <- apply(if(j%in%1:2) Qozone else Qtmean,1,fcumeff,0:20,combsimsens[j])
  
################################################################################
# LOOP ACROSS RANDOMLY SIMULATED DATA
  
  for(i in seq(nsim)) {
    
    # PRINT
    cat(i,"")
    
    # SET THE SEED
    seed <- 10000 + i
    set.seed(seed)
    
    # SIMULATE THE DATA (SUPPRESS THE WARNING DUE TO MISSING VALUES)
    #   - IN SENS 7, NEGATIVE BINOMIAL
    #   - OTHERWISE POISSON
    if(j!=7) {
      suppressWarnings(y <- rpois(length(cumeff),exp(timeeff+cumeff)))
    } else {
      # DEFINE THE OVERDISPERSION PHI
      phi <- 1.3
      # DEFINE THETA (PHI = 1 + MU/THETA, SO THETA = MU/(PHI-1) )
      theta <- exp(timeeff+cumeff)/(phi-1)
      # SIMULATE FROM A NEGATIVE BINOMIAL DISTRIBUTION
      suppressWarnings(y <- rnegbin(length(cumeff),exp(timeeff+cumeff),theta))
    }
    
################################################################################
# MODELS
    
    if(j %in% seq(argvarlist)) {
      
      # ARGUMENTS OF THE CROSSBASIS (DEPENDING ON THE SCENARIO AND MODEL)
      x <- if(j%in%1:2) ozone else tmean
      argvar <- argvarlist[[j]]
      lag <- laglist[j]
      arglag <- arglaglist[[j]]
      family <- if(j!=7) "poisson" else "quasipoisson"
      
      # OBTAIN THE CROSS-BASIS
      cb <- crossbasis(x,lag=lag,argvar=argvar,arglag=arglag)
      
      # RUN THE MODEL AND PREDICT
      model <- glm(y~cb+timebasis,family=family)
      at <- if(j%in%1:2) 0:50 else -20:35
      cen <- ifelse(j%in%1:2,0,20)
      cp <- crosspred(cb,model,at=at,cen=cen)
      
      # SELECT THE TRUE VALUE
      truecumsens <- truecum[[combsimsens[j]]]
      
    } else if(j%in%c(8,10)){
      
      # RUN THE MODEL AND PREDICT
      model <- glm(y~b07b+timebasis,family=poisson)
      cp <- crosspred(b07b,model,at=xpred07b,cen=20)
      
      # SELECT THE TRUE VALUE
      truecumsens <- truecum[[combsimsens[j]]][as.character(xpred07b)]
      
    } else if(j%in%c(9,11)){
    
      # RUN THE MODEL AND PREDICT
      model <- glm(y~b020b+timebasis,family=poisson)
      cp <- crosspred(b020b,model,at=xpred020b,cen=20)
      
      # SELECT THE TRUE VALUE
      truecumsens <- truecum[[combsimsens[j]]][as.character(xpred020b)]
    }

    # STORE THE RESULTS FOR OVERALL CUMULATIVE SUMMARY
    cumpredsens[[j]] <- cumpredsens[[j]] + cp$allfit
    cumbiassens[[j]] <- cumbiassens[[j]] + (cp$allfit - truecumsens)
    cumcovsens[[j]] <- cumcovsens[[j]] + (truecumsens >= 
        cp$allfit-qn*cp$allse & truecumsens <= cp$allfit+qn*cp$allse)
    cumrmsesens[[j]] <- cumrmsesens[[j]] + (cp$allfit - truecumsens)^2

  }
  
################################################################################
# COMPUTE THE AVERAGE

  cumpredsens[[j]] <- cumpredsens[[j]]/nsim
  cumbiassens[[j]] <- cumbiassens[[j]]/nsim
  cumcovsens[[j]] <- cumcovsens[[j]]/nsim
  cumrmsesens[[j]] <- sqrt(cumrmsesens[[j]]/nsim)
  
}

################################################################################
# SAVE

save.image("sens.RData")

#
