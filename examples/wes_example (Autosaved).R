library(FDRreg)

# Data
z_mat = read.csv('../data/0.csv', header=FALSE)
z = as.numeric(as.matrix(z_mat))

hist(z, 100, prob=TRUE)
curve(dnorm(x, 0, 3), add=TRUE, n=1000, col='red')

# Deconvolution under theoretical null
out1 = prfdr(z, control=list(npasses=25, gridsize=500))
plot(out1$x_grid, out1$pitheta_grid/(1-out1$pi0), type='l')

# Empirical null
f0 = efron(z, nulltype='empirical')
mu0 = f0$mu0
sig0 = f0$sig0

# Deconvolution under empirical null
out2 = prfdr(z, mu0, sig0, control=list(npasses=25, gridsize=250))
plot(out2$x_grid, out2$pitheta_grid/(1-out2$pi0), type='l')


# Compare convolution of estimate for pi(theta) with true f1
plot(out1$x_grid, out1$f1_grid)