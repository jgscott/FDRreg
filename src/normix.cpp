#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(Rcpp,RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

using namespace Rcpp;

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
    z = (y[i] - mu)/mysd;
    thisprob = (normalized_weights)*(dnorm(z)/mysd);
    ModelIndicator[i] = mysample(thisprob);
  }
  return ModelIndicator;
}


