library(FDRreg)

# Read in data (relative file location on GitHub repo)
synchrony_smithkohn2008 = read.csv('../data/synchrony_smithkohn2008.csv', header=TRUE)
z = synchrony_smithkohn2008$z

# Histogram of the z scores
hist(z, 250)

# Run Benjamini-Hochberg at the 0.1 level
BH1 = BenjaminiHochberg(z, 0.1)
sum(BH1 == 1 & z > 0)

# Estimate an empirical null using Efron's central matching procedure
e1 = efron(z, nulltype='empirical', df=10)

# Using the empirical null, fit the two-groups model via predictive recursion
mu0 = e1$mu0
sigma0 = e1$sig0
fhat = prfdr(z, mu0=mu0, sig0=sigma0, control=list(npasses=20,gridsize = 500))

# Posterior probabilities and FDR from two-groups model
postprob = fhat$postprob
BFDR = getFDR(fhat$postprob)
plot(z, postprob)
sum(BFDR$FDR <= 0.1 & z > 0)

####
# Fit FDR regression models
####

X = cbind(synchrony_smithkohn2008$Dist, synchrony_smithkohn2008$TuningCor)

# Set up spline basis functions (df=3)
df = 3
b1 = bs(synchrony_smithkohn2008$Dist, df=df)
b2 = bs(synchrony_smithkohn2008$TuningCor, df=df)
Xs = model.matrix( ~  b1 + b2  - 1)

# Empirical Bayes FDR regression
fdr1 = FDRreg(z, Xs, nulltype='empirical', control=list(lambda = 1))
FDR1 = getFDR(fdr1$postprob)
plotFDR(fdr1, breaks=150, mar=c(2,4,1,1))
sum(FDR1$FDR <= 0.1 & z > 0)

# Covariate effects
# distance
ind_dist = 2:(df+1)
beta_dist = (fdr1$model)$coef[ind_dist]
fhat_dist = b1 %*% beta_dist
fhat_dist_se = sqrt(diag(b1 %*% solve(fdr1$model$hessian)[ind_dist, ind_dist] %*% t(b1)))

# tuning-curve correlation
ind_tcc = (df+2):(2*df+1)
beta_tcc = (fdr1$model)$coef[ind_tcc]
fhat_tcc = b2 %*% beta_tcc
fhat_tcc_se = sqrt(diag(b2 %*% solve(fdr1$model$hessian)[ind_tcc, ind_tcc] %*% t(b2)))


# Full Bayes FDR regression
nmc=10000
burn=2000
N = length(z)
fdr3 = BayesFDRreg(z, Xs, mu0=mu0, sig0 = rep(sigma0, N), nmc=nmc, nburn=burn, ncomps=3,
            priorpars = list(PriorMean = rep(0, ncol(Xs) + 1), PriorPrec = diag(1, ncol(Xs) + 1)),
            control=list(verbose=500))

plotFDR(fdr3)
FDR3 = getFDR(fdr3$postprob)
sum(FDR3$FDR <= 0.1)


#####
##### Publication plots
##### Compare ordinary versus regression-based methods for FDR control
#####

mymat = matrix(0, nrow=6, ncol=2)
mymat[1:2,1] = 1:2
mymat[3:4,1] = 3
mymat[5:6,1] = 4
mymat[,2] = c(5,5,6,6,7,7)

layout(mymat)

# Initial evidence for a distance effect
synchrony_smithkohn2008$NaiveClass = (z > 2)
synchrony_smithkohn2008$NaiveClass = factor(synchrony_smithkohn2008$NaiveClass, labels=c('Z < 2', 'Z >= 2'))
par(mar=c(2,4,1,1))
hist(synchrony_smithkohn2008$Dist[synchrony_smithkohn2008$NaiveClass == 'Z < 2'], breaks=25,
	xlim=c(0,4600), col='lightgrey', border='darkgrey', axes=FALSE, 
	xlab='', ylab='Frequency', main='')
legend('top', 'Z < 2', bty='n')
axis(2, tick=FALSE, las=1)
par(mar=c(2,4,0,1))
hist(synchrony_smithkohn2008$Dist[synchrony_smithkohn2008$NaiveClass == 'Z >= 2'], breaks=20,
	xlim=c(0,4600), col='lightgrey', border='darkgrey', axes=FALSE, 
	xlab='Distance (micrometers)', ylab='Frequency', main='')
axis(2, tick=FALSE, las=1)
axis(1, tick=FALSE, line=-1)
legend('top', 'Z >= 2', bty='n')
mtext('Pairwise distance (micrometers)', 1, line=1.25, cex=0.7)

par(mar=c(4,4,5,4))
## Fitted densities
hist(z, 200, prob=TRUE, col='lightgrey', border='grey', xlim=c(-3,9),
	main='Z scores', axes=FALSE, xlab='')
lines(fhat$x_grid, (1-fhat$pi0)*fhat$f1_grid, col='blue', lwd=2, lty='dotted')
curve(fhat$pi0*dnorm(x, mu0, sigma0), col='red', add=TRUE, lwd=2, lty='solid')
legend('topright', bty='n', cex=1,
	legend=c( expression( c %.% f[1](z) ), expression( (1-c) %.% f[0](z) ) ),
	col=c('blue', 'red'), lty=c('dotted', 'solid'), lwd=2)
axis(1, tick=FALSE, line=-1)
axis(2, las=1, tick=FALSE)


# Compare discoveries under the methods
mybreaks = seq(min(fhat$x_grid) - 0.01, max(fhat$x_grid) + 0.01, length=80)
hist(subset(z, FDR1$FDR <= 0.1), mybreaks, col=rgb(0.9,0,0,0.7), border=NA,
	 main='', xlab='', ylab='Frequency', axes=FALSE, xlim=c(-3,9), ylim=c(0,90))
axis(1, tick=FALSE, line=-1)
axis(2, las=1, tick=FALSE)
par(new=TRUE)
hist(subset(z, BFDR$FDR <= 0.1), mybreaks, col=rgb(0,0,1,0.2), border=NA,
	 main='Discoveries at FDR = 0.1', xlab='', ylab='', axes=FALSE,
	 ylim=c(0,90), xlim=c(-3,9))
par(new=TRUE)
hist(subset(z, BFDR$FDR <= 0.1 & FDR1$FDR <= 0.1), mybreaks, col=rgb(.5,0,.5,0.8), border=NA,
	 main='Discoveries at FDR = 0.1', xlab='', ylab='', axes=FALSE,
	 ylim=c(0,90), xlim=c(-3,9))
par(new=FALSE)
legend("topright", c("FDRR only", "Two-groups only", "Both methods"),
	col=c(rgb(0.9,0,0,0.5), rgb(0,0,1,0.2), rgb(0.5,0,0.5,0.8)),
	bty='n', cex=0.85, pch=15)

# Extra discoveries under FDRR
disagree = which(BFDR$FDR >= 0.1 & FDR1$FDR <= 0.1)
plot(TuningCor ~ jitter(Dist,2), data=synchrony_smithkohn2008[which(FDR1$FDR <= 0.1),],
	pch=19, col=rgb(0,0,0,0.25), cex=0.3, xlim=c(0, 4500),
	main='Extra discoveries under FDRR',
	type='n',
	bty='n', ylab='Tuning Curve Correlation', xlab='Distance', axes=FALSE)
points(TuningCor ~ jitter(Dist,2), data=synchrony_smithkohn2008[disagree,],
	pch=19, col=rgb(1,0,0,0.6), cex=0.65)
axis(1, line=0, tick=FALSE)
axis(2, las=1, tick=FALSE)


# Covariate plots
# distance
oj = order(X[,1])
jind = 2:4
fj = tcrossprod({fdr3$X}[,jind],  fdr3$coefficients[,jind])
fjmean = rowMeans(fj)
fjquantile = apply(fj, 1, quantile, prob=c(0.025, 0.975))
plot(0, xlim=range(X[,1]), ylim=range(fjquantile), type='n', 
	main='Partial regression function: distance',
	bty='n', axes=FALSE, xlab='', ylab=expression(s[1](distance)))
lowerborder = fjquantile[1,oj]
upperborder = fjquantile[2,oj]
polygon(c(X[oj,1], rev(X[oj,1])), c(lowerborder, rev(upperborder)), border=NA, col='lightgrey')
lines(X[oj,1], fjmean[oj], lwd=2, col='black')
axis(1, tick=FALSE, line=-1)
axis(2, las=1, tick=FALSE)
mtext('Pairwise distance (micrometers)', 1, line=1.25, cex=0.7)

# tuning-curve correlation
oj = order(X[,2])
jind = 5:7
fj = tcrossprod({fdr3$X}[,jind],  fdr3$coefficients[,jind])
fjmean = rowMeans(fj)
fjquantile = apply(fj, 1, quantile, prob=c(0.025, 0.975))
plot(0, xlim=range(X[,2]), ylim=range(fjquantile), type='n', 
	bty='n', main='Partial regression function: tuning curve',
	axes=FALSE, xlab='', ylab=expression(s[2](correlation)))
lowerborder = fjquantile[1,oj]
upperborder = fjquantile[2,oj]
polygon(c(X[oj,2], rev(X[oj,2])), c(lowerborder, rev(upperborder)), border=NA, col='lightgrey')
lines(X[oj,2], fjmean[oj], lwd=2, col='black')
axis(1, tick=FALSE, line=-1)
axis(2, las=1, tick=FALSE)
mtext('Tuning Curve Correlation', 1, line=1.25, cex=0.7)	


#####
##### End Publication plots
#####
