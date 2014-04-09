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
    z = (mu - y[i])/mysd;
    thisprob = (normalized_weights)*(dnorm(z)/mysd);
    ModelIndicator[i] = mysample(thisprob);
  }
  return ModelIndicator;
}


// Dirk E.'s utility functions for subsetting

const double flagval = __DBL_MIN__; // works

// simple double value 'flagging' function
inline double flag(double a, bool b) { return b ? a : flagval; }

// [[Rcpp::export]]
NumericVector subsetter(NumericVector x, LogicalVector b) {
  // We use the flag() function to mark values of 'a' 
  // for which 'b' is false with the 'flagval'
  NumericVector a(clone(x));
  std::transform(a.begin(), a.end(), b.begin(), a.begin(), flag);

  // We use sugar's sum to compute how many true values to expect
  NumericVector res = NumericVector(sum(b));

  // And then copy the ones different from flagval from a into
  // res using the remove_copy function from the STL
  std::remove_copy(a.begin(), a.end(), res.begin(), flagval);
  return res;    
}


// // [[Rcpp::export]]
// List marmixEM(NumericVector y, int ncomps, double mu0, double sig0, double reltol) {
//   // y is an n-vectors of observations (y[i])
//   // k is the number of mixture components to be estimated
//   int ncases = y.size();

//   NumericVector weights(ncomps+1);
//   NumericVector means = rnorm(ncomps+1);
//   NumericVector vars(ncomps+1, 2.0);
//   double p0 = 0.1;
//   weights[ncomps] = p0;
//   means[ncomps] = mu0;
//   vars[ncomps] = sig0*sig0;
//   for(int k=0; k<ncomps; k++) {
//     weights[k] = (1.0-p0)/ncomps;
//   }

//   // Initialize the means
//   if(ncomps==1) means[0] = 0.0;
//   else {
//     means[0] = mean(subsetter(y, (y < -2.0)));
//     means[1] = mean(subsetter(y, (y > 2.0)));
//   }
  
//   NumericMatrix marglike(ncases, ncomps+1);
//   marglike(_,ncomps) = dnorm(y, mu0, sig0);
//   NumericMatrix postprob(ncases, ncomps+1);

//   double loglike = sum(log(dnormix(y, weights, means, vars+1.0)));
//   double loglikenew = loglike;
//   double rel_change;
//   bool converged=FALSE;

//   while(!converged) {

//     // E step
//     for(int k=0; k<ncomps; k++) {
//       marglike(_,k) = dnorm(y, means[k], sqrt(vars[k]));
//       postprob(_,k) = weights[k] * marglike(_,k);
//     }
//     postprob(_,ncomps) = weights[ncomps] * marglike(_,ncomps);
//     for(int i=0; i<ncases;  i++) {
//       postprob(i,_) = postprob(i,_)/sum( postprob(i,_) );
//     }

//     // M step
//     for(int k=0; k <= ncomps; k++) {
//       weights[k] = sum(postprob(_,k))/ncases;
//     }

//     for(int k=0; k < ncomps; k++) {
//       means[k] = sum(y*postprob(_,k))/sum(postprob(_,k));
//       NumericVector eps = (y - means[k]);
//       vars[k] = sum( eps*eps*postprob(_,k) ) / sum(postprob(_,k));
//       vars[k] = vars[k] > 1.0 ? vars[k] : 1.0;
//     }

//     loglikenew = sum(log(dnormix(y, weights, means, vars+1.0)));
//     rel_change = fabs(loglikenew - loglike)/fabs(loglike + reltol);
//     converged = (rel_change < reltol);
//     loglike = loglikenew;
//   }
  
//   return Rcpp::List::create(Rcpp::Named("weights") = weights,
// 			    Rcpp::Named("means") = means,
// 			    Rcpp::Named("vars") = vars );
// }

