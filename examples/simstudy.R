# This script runs the simulation study in the FDR regression paper:
# 	James G. Scott, Ryan C. Kelly, Matthew A. Smith, Pengcheng Zhou,
# 	and Robert E. Kass (2015).  False discovery rate regression:
# 	application to neural synchrony detection in primary visual
# 	cortex.  Journal of the American Statistical Association, DOI:
# 	10.1080/01621459.2014.990973. arXiv:1307.3495 [stat.ME].

# Please first install the package from GitHub.
# See the Readme file at https://github.com/jgscott/FDRreg/

# Then run the script in batch mode using `R CMD BATCH simstudy.R &' (without the quotes)
# This script will produce the output file 'sim2results-single.Rdata'

library(FDRreg)
library(locfdr)
library(foreach)
library(doMC)
library(parallel)
library(mvtnorm)
library(faraway)	# for logit transforms

registerDoMC(cores=parallel::detectCores())

#options(error=utils::recover)
options(error=NULL)

# Global simulation settings
nobserve = 10000
nsims = 100
nmc = 2500
burn = 1000
fdr_level = 0.10

# Some functions of (x1,x2) that will generate prior logs odds
flist = list()
flist[[1]] = function(x1, x2) return(-3.0 + 1.5*x1 + 1.5*x2)
flist[[2]] = function(x1, x2) return(-3.25 + 3.5*x1^2 - 3.5*x2^2)
flist[[3]] = function(x1, x2) return(0 - 1.5*(x1-0.5)^2 - 5*abs(x2) )
flist[[4]] = function(x1, x2) return(-4.25 + 2*x1^2 + 2*x2^2 - 2*x1*x2)
flist[[5]] = function(x1, x2) return(-3)

# Some parameter settings for the mixture model that will generate signals
parlist = list()
parlist[[1]] = list(weights=c(0.48,0.04,0.48), mu=c(-2,0,2), tau2 = c(1,16,1))
parlist[[2]] = list(weights=c(0.4,0.2,0.4), mu=c(-1.25,0,1.25), tau2 = c(2,4,2))
parlist[[3]] = list(weights=c(0.3,0.4,0.3), mu=-c(0,0,0), tau2 = c(0.1,1,9))
parlist[[4]] = list(weights=c(0.2,0.3,0.3,0.2), mu=c(-3,-1.5,1.5,3), tau2 = c(0.01,.01,.01,0.01))


# # Plot logodds versus covariates
par(mar=c(3,3,3,1))
mylevs = c(10,10,8,10)
layout(matrix(1:4,2,2))
 for(j in 1:4) {
 myxfunction = flist[[j]]
 x1 = seq(-1,1,length=100)
 x2 = seq(-1,1,length=100)
 D = expand.grid(x1,x2)
 logodds = matrix(myxfunction(D[,1],D[,2]), nrow=100, ncol=100)
 contour(x=x1,y=x2,z=logodds, nlevels=mylevs[j], las=1, col='black',
    method='edge', labcex=0.7, main=paste0('Function ', LETTERS[j]))

 x1 = runif(nobserve, -1,1)
 x2 = runif(nobserve, -1,1)
 X = cbind(x1,x2)
 log_odds = myxfunction(x1,x2)
 success_prob = faraway::ilogit(log_odds)
 is_signal = rbinom(nobserve, 1, success_prob)
 cat(sum(is_signal)/nobserve, "\n")
 #points(x1,x2, col=rgb(.1,.1,.1,.1))
 # points(x1[is_signal==1], x2[is_signal==1], pch=19, col=rgb(0,0,0,.05), cex=0.5)
}




# Plot signal distributions
par(mar=c(4,3,3,1))
layout(matrix(1:4,2,2))
for(j in 1:4) {
   mypars = parlist[[j]]
   curve(dnormix(x, mypars$weights, mypars$mu, mypars$tau2), from=-5,to=5, n=1001L,
      bty='n', xlab=expression(theta), ylab='Density of signals', las=1,
      main=paste0('Case ', j))
}



results_list = list()
counter = 1

for(i in 3:length(flist)) {
   
   cat("Function", i, "\n")
   myxfunction = flist[[i]]
   
   for(j in seq_along(parlist)) {
      cat("\t Signal Density", j, "\n")
      mypars = parlist[[j]]

      # Parallel loop over simulations
      #thissim = matrix(0, nrow=nsims, ncol=10)
      #thissim = matrix(0, nrow=nsims, ncol=8)
      thissim = foreach(mc = 1:nsims, .combine='rbind') %dopar% {
         
         # Generate data
         x1 = runif(nobserve, -1,1)
         x2 = runif(nobserve, -1,1)
         X = cbind(x1,x2)
         log_odds = myxfunction(x1,x2)
         success_prob = faraway::ilogit(log_odds)
         is_signal = rbinom(nobserve, 1, success_prob)
         mu = FDRreg::rnormix(nobserve, mypars$weights, mypars$mu, mypars$tau2)
         mu[is_signal==0] = 0
         y = mu + rnorm(nobserve,0,1)
         
         # Spline design matrix
         Xs=model.matrix( ~ ns(X[,1], df=3, intercept=FALSE) + ns(X[,2], df=3, intercept=FALSE) - 1)
                  
         # Benjamini-Hochberg
         is_findingBH = BenjaminiHochberg(y, fdr_level)
         rates_BH = GetErrorRates(is_signal, is_findingBH)
         
         # Efron's method
         out1 = locfdr(y, nulltype=0, plot=0)
         efronFDR = getFDR(1-out1$fdr)
         is_findingEfron = 0 + {efronFDR$FDR <= fdr_level}
         rates_Efron = GetErrorRates(is_signal, is_findingEfron)

         # EB FDR regression with marginalization
         out2 = FDRreg(y, Xs, nulltype='theoretical', method='efron')
         findings = which(out2$FDR <= fdr_level)
         is_finding = rep(0, nobserve)
         is_finding[findings] = 1
         rates_FDRreg = GetErrorRates(is_signal, is_finding)
         
         # EB FDR regression
         out3 = FDRreg(y, Xs, nulltype='theoretical', method='pr')
         findings = which(out3$FDR <= fdr_level)
         is_finding = rep(0, nobserve)
         is_finding[findings] = 1
         rates_EBayesFDRreg = GetErrorRates(is_signal, is_finding)
         
         # Full Bayes FDR regression
         out4 = BayesFDRreg(y, Xs, nmc=nmc, nburn=burn, ncomps=3, nulltype='theoretical',
            priorpars = list(PriorMean = rep(0, ncol(Xs) + 1), PriorPrec = diag(1/10, ncol(Xs) + 1)))
         findings = which(out4$FDR <= fdr_level)
         is_finding = rep(0, nobserve)
         is_finding[findings] = 1
         rates_FullBayesFDRreg = GetErrorRates(is_signal, is_finding)
   
         c(   unlist(rates_BH[1:2])
            , unlist(rates_Efron[1:2])
            , unlist(rates_FDRreg[1:2])
            , unlist(rates_EBayesFDRreg[1:2])
            , unlist(rates_FullBayesFDRreg[1:2])
         )
      }
      results_list[[counter]] = list(func=i, signaldens=j, result=thissim)
      counter = counter+1
   }
}



save(results_list, file='sim2results-single.Rdata')

evens = seq(2, ncol(results_list[[1]]$result), by=2)
odds = seq(1, ncol(results_list[[1]]$result)-1, by=2)
for(j in seq_along(results_list)) {
   cat(round(colMeans(results_list[[j]]$result[,c(evens, odds)]),2), sep=" ", eol='\n')
}

