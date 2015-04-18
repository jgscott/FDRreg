// Copyight James G. Scott 2014
// Published under GPL3 licence

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(Rcpp,RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

using namespace Rcpp;

// ***************************************
// Basic utitilies used by the main functions
// ***************************************

// [[Rcpp::export]]
int mysample(NumericVector probs) {
  NumericVector csum = cumsum(probs);
  double u = Rf_runif(0.0,csum[csum.size()-1]);
  int k=0;
  while(u > csum[k]) {
    k++;
  }
  return k;
}

// [[Rcpp::export]]
double trapezoid(NumericVector x, NumericVector y) {
  int n = x.size();
  double result = 0.0;
  for(int i=1; i<n; i++) {
    result += (x[i] - x[i-1])*(y[i] + y[i-1])/2.0;
  }
  return result;
}

// ***************************************
// The main workhorse functions are below
// ***************************************


// [[Rcpp::export]]
NumericVector dnormix(NumericVector y, NumericVector weights,  NumericVector mu, NumericVector tau2) {
  // y is an n-vectors of observations (y[i])
  // weights, mu, and tau2 are p-vectors of weights, means, and variances in a K-component mixture model
  // returns the density of y[i] under the mixture model
  int ncases = y.size();
  int ncomps = mu.size();
  NumericVector normalized_weights = weights/sum(weights);
  NumericVector density(ncases, 0.0);
  NumericVector thisdens(ncases);
  NumericVector sigma = sqrt(tau2);
  for(int j=0; j < ncomps; j++) {
    thisdens = normalized_weights[j]*dnorm(y, mu[j], sigma[j]);
    for(int i=0; i < ncases; i++) {
      density[i] += thisdens[i];
    }
  }
  return density;
}


// [[Rcpp::export]]
NumericVector marnormix(NumericVector y, NumericVector sigma2, NumericVector weights,  NumericVector mu, NumericVector tau2) {
  // y and sigma2 are n-vectors of observations (y[i]) and observation-level variances (sigma2[i])
  // weights, mu, and tau2 are p-vectors of weights, means, and variances in a K-component mixture model
  // returns the marginal/predictive density of y[i] under a mixture of normals prior and Gaussian observation model
  int ncases = y.size();
  int ncomps = mu.size();
  NumericVector normalized_weights = weights/sum(weights);
  NumericVector density(ncases, 0.0);
  NumericVector thisdens(ncases, 0.0);
  NumericVector mysd(ncomps, 0.0);
  NumericVector z(ncases);
  for(int j=0; j < ncomps; j++) {
    NumericVector mysd = sqrt(sigma2 + tau2[j]);
    z = (y-mu[j])/mysd;
    thisdens = (normalized_weights[j])*(dnorm(z)/mysd);
    for(int i=0; i < ncases; i++) {
      density[i] += thisdens[i];
    }
  }
  return density;
}


// [[Rcpp::export]]
NumericVector rnormix(int n, NumericVector weights,  NumericVector mu, NumericVector tau2) {
  // n is the desired number of draws
  // weights, mu, and tau2 are p-vectors of weights, means, and variances in a K-component mixture model
  // returns a sample of size n from the corresponding mixture model
  int j;
  NumericVector result(n);
  NumericVector tau = sqrt(tau2);
  for(int i=0; i < n; i++) {
    j = mysample(weights);
    result[i] = Rf_rnorm(mu[j], tau[j]);
  }
  return result;
}


// [[Rcpp::export]]
IntegerVector draw_mixture_component(NumericVector y, NumericVector sigma2, NumericVector weights,  NumericVector mu, NumericVector tau2) {
  // y and sigma2 are n-vectors of observations (y[i]) and observation-level variances (sigma2[i])
  // weights, mu, and tau2 are p-vectors of weights, means, and variances in a K-component mixture model
  int ncases = y.size();
  int ncomps = mu.size();
  IntegerVector choices = seq_len(ncomps);
  NumericVector normalized_weights = weights/sum(weights);
  NumericVector mysd(ncomps, 0.0);
  IntegerVector ModelIndicator(ncases, 0);
  NumericVector thisprob(ncomps,0.0);
  NumericVector z(ncomps);
  for(int i=0; i < ncases; i++) {
    mysd = sqrt(tau2 + sigma2[i]);
    z = (mu - y[i])/mysd;
    thisprob = (normalized_weights)*(dnorm(z)/mysd);
    ModelIndicator[i] = mysample(thisprob);
  }
  return ModelIndicator;
}


// [[Rcpp::export]]
List PredictiveRecursionFDR(NumericVector z, IntegerVector sweeporder,
    NumericVector grid_x, NumericVector theta_guess,
    double mu0 = 0.0, double sig0 = 1.0, double nullprob=0.95, double decay = -0.67) {
  // z: z statistics
  // sweeporder: an ordering of the points in z, usually 5-10 stacked permutations of 1 ... n
  // grid_x: a grid of points at which the alternative density will be approximated
  // theta_guess: an initial guess for the sub-density under the alternative hypothesis
  // nullprob: an initial guess for the fraction of null cases
  // decay: the stochastic-approximation decay parameter, should be in (-1, -2/3)

  // Set-up
  int n = sweeporder.size();
  int k, gridsize = grid_x.size();
  NumericVector theta_subdens(clone(theta_guess));
  double pi0 = nullprob;
  NumericVector joint1(gridsize);
  NumericVector ftheta1(gridsize);
  double m0, m1, mmix, cc;

  // Begin sweep through the data
  for(int i=0; i<n; i++) {
    k = sweeporder[i];
    if(i % 200 == 0) Rcpp::checkUserInterrupt();  
    cc = pow(3.0+(double)i, decay);
    joint1 = dnorm(grid_x, z[k]-mu0, sig0) * theta_subdens;
    m0 = pi0*R::dnorm(z[k] - mu0, 0.0, sig0, 0);
    m1 = trapezoid(grid_x, joint1);
    mmix = m0 + m1;
    pi0 = (1.0-cc)*pi0 + cc*(m0/mmix);
    ftheta1 = joint1/mmix;
    theta_subdens = (1.0-cc)*theta_subdens + cc*ftheta1;
  }

  // Now calculate marginal distribution along the grid points
  NumericVector y_mix(gridsize);
  NumericVector y_signal(gridsize);
  for(int j=0; j<gridsize; j++) {
    joint1 = dnorm(grid_x, grid_x[j] - mu0, sig0)*theta_subdens;
    m0 = pi0*R::dnorm(grid_x[j], mu0, sig0, 0);
    m1 = trapezoid(grid_x, joint1);
    y_mix[j] = m0 + m1;
    y_signal[j] = m1/(1.0-pi0);
  }

  return Rcpp::List::create(Rcpp::Named("grid_x")=grid_x,
          Rcpp::Named("theta_subdens")=theta_subdens,
          Rcpp::Named("pi0")=pi0,
          Rcpp::Named("y_mix")=y_mix,
          Rcpp::Named("y_signal")=y_signal
          );
}



// [[Rcpp::export]]
List PredictiveRecursion_DifferentSigma(NumericVector z, double mu0, NumericVector sig0, IntegerVector sweeporder,
    NumericVector grid_x, NumericVector theta_guess,
    double nullprob=0.95, double decay = -0.67) {
  // z: vector of observed data, vector of length N
  // sig0: vector of standard errors se(z_i) under the null, of same length as z
  // sweeporder: a vector of indices in [0...N-1] representing the order in which the z are processed
  //     length(sweeporder) = Npasses*length(z)
  // grid_x: a grid of points at which the alternative density will be approximated
  // theta_guess: an initial guess for the sub-density under the alternative hypothesis
  // nullprob: an initial guess for the fraction of null cases
  // decay: the stochastic-approximation decay parameter, should be in (-1, -2/3)

  // Set-up
  int n = sweeporder.size();
  int k, gridsize = grid_x.size();
  NumericVector theta_subdens(clone(theta_guess));
  double pi0 = nullprob;
  NumericVector joint1(gridsize);
  NumericVector ftheta1(gridsize);
  double m0, m1, mmix, cc;

  // Begin sweep through the data
  for(int i=0; i<n; i++) {
    if(i % 200 == 0) Rcpp::checkUserInterrupt();  
    k = sweeporder[i];
    cc = pow(3.0+(double)i, decay);
    joint1 = dnorm(grid_x, z[k] - mu0, sig0[k]) * theta_subdens;
    m0 = pi0*R::dnorm(z[k] - mu0, 0.0, sig0[k], 0);
    m1 = trapezoid(grid_x, joint1);
    mmix = m0 + m1;
    pi0 = (1.0-cc)*pi0 + cc*(m0/mmix);
    ftheta1 = joint1/mmix;
    theta_subdens = (1.0-cc)*theta_subdens + cc*ftheta1;
  }

  return Rcpp::List::create(Rcpp::Named("grid_x")=grid_x,
          Rcpp::Named("theta_subdens")=theta_subdens,
          Rcpp::Named("pi0")=pi0
          );
}

// [[Rcpp::export]]
NumericVector GaussianConvolution(NumericVector x, NumericVector fx, double sigma) {
  // x: vector of grid points
  // fx: a vector of density evaluations f(x) at the grid points
  // sigma: the standard deviation of the Gaussian with which f(x) is to be convolved
  // this is a very naive O(n^2) Gaussian convolution suitable for small cheap problems
  // Use something cleverer (e.g. FFT) for bigger problems!

  // Set-up
  int n = x.size();
  NumericVector result(n ,0.0);

  // Sweep through the grid points
  for(int i=0; i<n; i++) {
    result[i] += sum(fx*dnorm(x, x[i], sigma))/sum(dnorm(x,x[i],sigma));
  }

  return result;
}


// [[Rcpp::export]]
List eval_pr_dens(NumericVector z, double mu0, NumericVector sig0, NumericVector grid_x, NumericVector grid_theta) {
  // z: z statistics
  // mu0: mean of Gaussian
  // sig0: vector of standard errors, se(z[i])
  // grid_x: a grid of points at which the alternative density pi(theta) is evaluated
  // grid_theta: the (unnormalized) mixing density pi(theta) at each point in grid_x
  // this function will evaluate the predictive density of each z point

  // Set-up
  int n = z.size();
  int gridsize = grid_x.size();
  NumericVector fsignal_z(n);
  NumericVector joint1(gridsize);
  double norm_constant = trapezoid(grid_x, grid_theta);

  // Begin sweep through the data, each time integrating by trap rule
  for(int i=0; i<n; i++) {
    if(i % 200 == 0) Rcpp::checkUserInterrupt();  
    joint1 = dnorm(grid_x, z[i] - mu0, sig0[i]) * grid_theta;
    fsignal_z[i] = trapezoid(grid_x, joint1)/norm_constant;
  }
  return Rcpp::List::create(Rcpp::Named("fsignal_z")=fsignal_z
          );
}


