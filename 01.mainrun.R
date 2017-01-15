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
# RUN THE MAIN ANALYSIS
################################################################################

################################################################################
# CREATE THE OBJECT TO STORE RESULTS

# MODELS
mod <- c('DL','MA0-3','MA0-7','MA0-20')

# LISTS FOR STORING PREDICTIONS AND SAMPLES FOR THE WHOLE SURFACE
# DLNM ONLY
pred <- lapply(seq(combsim),function(j) lapply(mod, function(m) 0))
sample <- lapply(seq(combsim),function(j) vector('list',length(mod)))
names(pred) <- names(sample) <- combsim
for(i in seq(pred)) names(pred[[i]]) <- names(sample[[i]]) <- mod

# LISTS FOR STORING CUMULATIVE PREDICTION, BIAS, COVERAGE, RMSE AND SAMPLES
# ALL THE MODELS
cumpred <- pred
cumsample <- sample
cumbias <- cumcov <- cumrmse <- cumpred

################################################################################
# START THE ITERATIONS
# 3 LEVELS:
#   - TYPE OF SURFACE (COMBINATIONS OF FUNCTIONS USED TO SIMULATE EFFECTS)
#   - RANDOM GENERATED RESPONSE
#   - MODEL

# LOOP ACROSS SIMULATED EXPOSURE-LAG-RESPONSES
for(j in seq(combsim)) {
  
  # PRINT
  cat("\n\n",names(combsim)[j],"\n")
 
  # DEFINE THE TRUE EFFECT (DEPENDENT ON SCENARIO)
  cumeff <- apply(if(j%in%1:2) Qozone else Qtmean,1,fcumeff,0:20,combsim[j])

################################################################################
# LOOP ACROSS RANDOMLY SIMULATED DATA
  
  for(i in seq(nsim)) {
    
    # PRINT
    cat(i,"")
    
    # SET THE SEED
    seed <- 10000 + i
    set.seed(seed)
    
    # SIMULATE THE DATA (SUPPRESS THE WARNING DUE TO MISSING VALUES)
    suppressWarnings(y <- rpois(length(cumeff),exp(timeeff+cumeff)))
    
################################################################################
# LOOP ACROSS MODELS
    
    for(m in seq(mod)) {
      
      # ARGUMENTS OF THE CROSSBASIS (DEPENDING ON THE SCENARIO AND MODEL)
      x <- if(j%in%1:2) ozone else tmean
      argvar <- if(j%in%1:2) list("lin") else 
        list("bs",degree=2,knots=quantile(x,c(10,75,90)/100,na.rm=T))
      lag <- c(20,3,7,20)[m]
      arglag <- if(m==1) list(knots=logknots(lag,nk=3)) else list("strata",df=1)
      
      # OBTAIN THE CROSS-BASIS
      cb <- crossbasis(x,lag=lag,argvar=argvar,arglag=arglag)
      
      # RUN THE MODEL AND PREDICT
      model <- glm(y~cb+timebasis,family=poisson)
      at <- if(j%in%1:2) 0:50 else -20:35
      cen <- ifelse(j%in%1:2,0,20)
      cp <- crosspred(cb,model,at=at,cen=cen)
      
      # MEASURE SIGNAL-TO-NOISE
      #cor(na.omit(log(y+1)),log(predict(model,type="response")+1))
      
      # STORE THE RESULTS FOR WHOLE SURFACE
      pred[[j]][[m]] <- pred[[j]][[m]] + cp$matfit
      if(i<=nsample) sample[[j]][[m]] <- c(sample[[j]][[m]],list(cp$matfit))
      
      # STORE THE RESULTS FOR OVERALL CUMULATIVE SUMMARY
      cumpred[[j]][[m]] <- cumpred[[j]][[m]] + cp$allfit
      cumbias[[j]][[m]] <- cumbias[[j]][[m]] + (cp$allfit - truecum[[j]])
      cumcov[[j]][[m]] <- cumcov[[j]][[m]] + (truecum[[j]] >= 
          cp$allfit-qn*cp$allse & truecum[[j]] <= cp$allfit+qn*cp$allse)
      cumrmse[[j]][[m]] <- cumrmse[[j]][[m]] + (cp$allfit - truecum[[j]])^2
      if(i<=nsample) cumsample[[j]][[m]] <- c(cumsample[[j]][[m]],
        list(cp$allfit))

    }
  }
    
################################################################################
# COMPUTE THE AVERAGE
  
  for(m in seq(mod)) {
    pred[[j]][[m]] <- pred[[j]][[m]]/nsim
    cumpred[[j]][[m]] <- cumpred[[j]][[m]]/nsim
    cumbias[[j]][[m]] <- cumbias[[j]][[m]]/nsim
    cumcov[[j]][[m]] <- cumcov[[j]][[m]]/nsim
    cumrmse[[j]][[m]] <- sqrt(cumrmse[[j]][[m]]/nsim)
  }
}

################################################################################
# SAVE

save.image("main.RData")

#
