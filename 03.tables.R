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
# TABLES
################################################################################

################################################################################
# TABLE 1

tab1 <- formatC(
  do.call(cbind,lapply(seq(combsim), function(j) {
    covsub <- !(if(j%in%1:2) 0:50 %in% 0 else -20:35 %in% 20)
    t(sapply(seq(mod)[1:4], function(i) c(
      mean(abs(cumbias[[j]][[i]]))*100,
      mean(cumcov[[j]][[i]][covsub]),
      mean(cumrmse[[j]][[i]])*100
    )))
  })),format="f",digits=2)

rownames(tab1) <- mod
colnames(tab1) <- c(t(outer(combsim,c('bias','cov','rmse'),paste,sep='-')))

write.csv(tab1,file="table1.csv")

################################################################################
# TABLE 2

tab2 <- formatC(
  do.call(rbind,lapply(seq(combsimsens[1:7]), function(j) {
    covsub <- !(names(cumcovsens[[j]]) %in% if(j%in%1:2) 0 else  20)
    c(
      mean(abs(cumbiassens[[j]]))*100,
      mean(cumcovsens[[j]]),
      mean(cumrmsesens[[j]])*100
    )
  })),format="f",digits=2)
tab2 <- cbind(combsimsens[1:7],"",tab2)

rownames(tab2) <- sens[1:7]
colnames(tab2) <- c('Scenario','Model','bias','cov','rmse')

write.csv(tab2,file="table2.csv")

#
