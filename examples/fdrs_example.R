library(FDRreg)


# Now try logit
flogit = function(x) log(x/(1-x))
ilogit = function(x) 1/{1+exp(-x)}

n = 5000

# How frequent are signals in the "on" and "off" regions?
w_on = 0.9
w_off = 0.01

# Z scores
beta_true = rep(flogit(w_off),n) 
beta_true[2251:2750] = flogit(w_on)
w_true = ilogit(beta_true)
wbar = mean(w_true)
gamma_true = rbinom(n, 1, w_true)
z = rnorm(n)

# Z scores under alternative
z[gamma_true == 1] = rnorm(sum(gamma_true), 0, 3)


out2 = FDRsmooth1D(z)
plot(out2$lambda, out2$bic, log='xy')
jbest = which.min(out2$bic)

par(mfrow=c(1,2))
plot(ilogit(out2$theta[,jbest]), ylim=c(0,1), type='l')
lines(w_true, col='blue')
abline(h=wbar, lty='dotted')

plot(out2$prfit$x_grid, out2$prfit$f1_grid) 
