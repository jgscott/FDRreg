#include "RcppArmadillo.h"
#include <R_ext/Utils.h>
#include <FDRreg.h>
#include <iostream>
#include <exception>
#include "RNG.h"
#include "PolyaGamma.h"

// Rcpp::depends(RcppArmadillo)
// [[Rcpp::interfaces(r, cpp)]]

using namespace Rcpp;
using namespace arma;
using namespace FDRreg;

template <typename T>
T square( const T& x) {
  return x * x;
}

// simulate a single mean-zero multivariate normal vector
vec rmvnormArma(mat sigma) {
   int ncols = sigma.n_cols;
   vec z = arma::randn(ncols);
   return chol(sigma) * z;
}


// [[Rcpp::export]]
double rtgamma_once(double shape, double rate, double lb, double ub) {
  // Draws from a truncated Gamma(shape, rate) distribution on (lb, ub)
  // we pass in the shape and rate parameters of the gamma
  // these are converted internally to scale parameters used by Rf_{dpq}gamma
  double lub = Rf_pgamma(lb, shape, 1.0/rate, 1,0);
  double uub = Rf_pgamma(ub, shape, 1.0/rate, 1,0);
  double u = Rf_runif(lub, uub);
  double draw = Rf_qgamma(u, shape, 1.0/rate, 1,0);
  return draw;
}

// Draw from a Dirichlet distribution with parameter alpha
// [[Rcpp::export]]
NumericVector rdirichlet_once(NumericVector alpha) {
  int d = alpha.size();
  NumericVector draw(d);
  for(int i=0; i<d; i++) {
    draw[i] = Rf_rgamma(alpha[i], 1.0);
  }
  return draw/sum(draw);
}



// compute the inverse logit transform for an arma::vec
colvec InvLogit(arma::vec x) {
  int d = x.n_elem;
  colvec result(d);
  for(int i=0; i<d; i++) {
    result(i) = 1.0/(1.0+exp(-x(i)));
  }
  return result;
}


// priorprob = vector of prior probabilities
// mfull = vector of marginal likelihoods under global predictive (mixture of null + alternative)
// mnull = vector of likelihoods under the null predictive
colvec PriorToPostprobEBayes(const arma::vec & priorprob, NumericVector mfull, NumericVector mnull, double p0) {
  int n = priorprob.n_elem;
  colvec postprob(n);
  for(int i=0; i<n; i++) {
    if(p0*mnull[i] >= mfull[i]) {
      postprob(i) = 0.0;
    } else {
      postprob(i) =  std::max(0.0, 1.0-(1.0-priorprob(i))*(mnull[i]/mfull[i]));
    }
  }
  return postprob;
}

// priorprob = vector of prior probabilities
// m1 = vector of marginal likelihoods under alternative hypothesis
// m0 = vector of likelihoods under the null predictive
colvec PriorToPostprobFullBayes(const arma::vec & priorprob, NumericVector m1, NumericVector m0) {
  int n = priorprob.n_elem;
  colvec postprob(n);
  double f0, f1;
  for(int i=0; i<n; i++) {
    f1 = priorprob(i) * m1[i];
    f0 = (1.0-priorprob(i)) * m0[i];
    postprob(i) = f1/(f0+f1);
  }
  return postprob;
}


colvec rpg(colvec ntrials, colvec logodds) {
// Draws random PG variances from arma::vectors of ntrials and logodds

  // Set up the RNG
  RNG r;
  PolyaGamma pg;
#ifdef USE_R
  GetRNGstate();
#endif
  int d = ntrials.n_elem;
  colvec result(d);
  for(int i=0; i<d; i++) {
    result[i] = pg.draw(ntrials(i), logodds(i), r);
  }

#ifdef USE_R
  PutRNGstate();
#endif
  
  return result;
}


// [[Rcpp::export]]
SEXP FDRregCPP(NumericVector z, const arma::mat & X, NumericVector M0, NumericVector MTot, const arma::mat & PriorPrecision, const arma::vec & PriorMean, int nmc, int nburn, double p0, const arma::vec & betaguess) {
  // z is an n-vector of test statistics
  // X is a design matrix in the FDRR problem, including the intercept
  // M0 is a vector of marginal likelihoods, f_0(z), under the null hypothesis
  // MTot is a vector of marginal likelihoods, f(z), under the mixture of null and alternative
  // PriorPrecision and PriorMean are the parameters of the multivariate normal prior for the regression coefficients
  // nmc and nburn say how long to run the Markov chain

  RNGScope scope;

  // Precompute what we can
  int ncases = X.n_rows;
  int nfeatures = X.n_cols;
  colvec PriorPrecxMean = PriorPrecision * PriorMean;
  colvec pgscale(ncases, fill::ones);
  colvec onehalf(ncases);
  onehalf.fill(0.5);

  // Storage for draws
  mat betasave(nmc, nfeatures);
  colvec priorprobsave(ncases, fill::zeros);
  colvec postprobsave(ncases, fill::zeros);

  // Initialize MCMC
  int nsamples = nmc + nburn;
  mat PosteriorPrecision(nfeatures, nfeatures);
  mat PosteriorVariance(nfeatures, nfeatures);
  colvec beta = betaguess;
  colvec logodds(ncases, fill::zeros);
  colvec priorprob(ncases);
  colvec postprob(ncases, fill::randu);
  colvec omega(ncases);
  colvec PosteriorMean(nfeatures, fill::zeros);

  // Main MCMC loop
  // bool interrupt = false;
  for(int i=0; i<nsamples; i++) {

 //    // Check for R user interrupt
 //    if(i % 100 == 0) {
 //      if (check_interrupt()) {
	// interrupt = true;
 //      }
 //      // throw exception if interrupt occurred
 //      if (interrupt) {
	// throw interrupt_exception("The FDRR model fit was interrupted.");
 //      }
 //    }

    if(i % 100 == 0) Rcpp::checkUserInterrupt();

    // compute prior and posterior probability of alternative hypothesis
    logodds = X*beta;
    priorprob = InvLogit(logodds);
    postprob = PriorToPostprobEBayes(priorprob, MTot, M0, p0);

    // Draw Polya-Gamma latents, given current regression estimate of log odds
    omega = rpg(pgscale, logodds);

    // Update regression parameters
    PosteriorPrecision = X.t() * diagmat(omega) * X + PriorPrecision;
    PosteriorVariance = inv(PosteriorPrecision);
    PosteriorMean = PosteriorVariance*(X.t() * (postprob - onehalf) + PriorPrecxMean);
    beta = PosteriorMean + rmvnormArma(PosteriorVariance);
    if(i >= nburn) {
      betasave.row(i-nburn) = beta.t();
      priorprobsave = priorprobsave + priorprob;
      postprobsave = postprobsave + postprob;
    }
  }
  postprobsave = (1.0/nmc)*postprobsave;
  priorprobsave = (1.0/nmc)*priorprobsave;

  return Rcpp::List::create(Rcpp::Named("priorprob")=priorprobsave,
			    Rcpp::Named("postprob")=postprobsave,
			    Rcpp::Named("betasave")=betasave
			    );
}



// [[Rcpp::export]]
SEXP EmpiricalBayesFDRregCPP(NumericVector z, const arma::mat & X, NumericVector M0, NumericVector M1, const arma::mat & PriorPrecision, const arma::vec & PriorMean, int nmc, int nburn, const arma::vec & betaguess) {
  // z is an n-vector of test statistics
  // X is a design matrix in the FDRR problem, including the intercept
  // M0 is a vector of marginal likelihoods, f_0(z), under the null hypothesis
  // M1 is a vector of marginal likelihoods, f_1(z), under the alternative hypothesis
  // PriorPrecision and PriorMean are the parameters of the multivariate normal prior for the regression coefficients
  // nmc and nburn say how long to run the Markov chain

  GetRNGstate();
  RNGScope scope;

  // Precompute what we can
  int ncases = X.n_rows;
  int nfeatures = X.n_cols;
  colvec PriorPrecxMean = PriorPrecision * PriorMean;
  colvec pgscale(ncases, fill::ones);
  colvec onehalf(ncases);
  onehalf.fill(0.5);
  
  // Storage for draws
  mat betasave(nmc, nfeatures);
  colvec priorprobsave(ncases, fill::zeros);
  colvec postprobsave(ncases, fill::zeros);

  // Initialize MCMC
  mat PosteriorPrecision(nfeatures, nfeatures);
  mat PosteriorVariance(nfeatures, nfeatures);
  colvec beta(nfeatures, fill::zeros);
  colvec omega(ncases);
  colvec PosteriorMean(nfeatures, fill::zeros);

  // Main MCMC loop
  int nsamples = nmc + nburn;
  for(int i=0; i<nsamples; i++) {

    // // Check for R user interrupt
    // if(i % 100 == 0) {
    //   if (check_interrupt()) {
	   //     interrupt = true;
    //   }
    //   // throw exception if interrupt occurred
    //   if (interrupt) {
    //     throw interrupt_exception("The FDRR model fit was interrupted.");
    //   }
    // }

    if(i % 100 == 0) Rcpp::checkUserInterrupt();

    // Update logodds/prior/posterior probability from regression model, M0, and M1
    colvec logodds = X*beta;
    colvec priorprob = InvLogit(logodds);
    colvec postprob = PriorToPostprobFullBayes(priorprob, M1, M0);
    
    // Draw Polya-Gamma latents, given current regression estimate of log odds
    omega = rpg(pgscale, logodds);

    // Update regression parameters
    PosteriorPrecision = X.t() * diagmat(omega) * X + PriorPrecision;
    PosteriorVariance = inv(PosteriorPrecision);
    PosteriorMean = PosteriorVariance*(X.t() * (postprob - onehalf) + PriorPrecxMean);
    beta = PosteriorMean + rmvnormArma(PosteriorVariance);

    if(i >= nburn) {
      betasave.row(i-nburn) = beta.t();
      priorprobsave = priorprobsave + priorprob;
      postprobsave = postprobsave + postprob;
    }
  }
  postprobsave = (1.0/nmc)*postprobsave;
  priorprobsave = (1.0/nmc)*priorprobsave;

  PutRNGstate();

  return Rcpp::List::create(Rcpp::Named("priorprob")=priorprobsave,
			    Rcpp::Named("postprob")=postprobsave,
			    Rcpp::Named("betasave")=betasave
			    );
}






// [[Rcpp::export]]
SEXP FullyBayesFDRregCPP(NumericVector z, const arma::mat & X, NumericVector M0, NumericVector M1, const arma::mat & PriorPrecision, const arma::vec & PriorMean, int nmc, int nburn, const arma::vec & betaguess) {
  // z is an n-vector of test statistics
  // X is a design matrix in the FDRR problem, including the intercept
  // M0 is a vector of marginal likelihoods, f_0(z), under the null hypothesis
  // ncomps the number of components you want in the Gaussian mixture model for the alternative
  // PriorPrecision and PriorMean are the parameters of the multivariate normal prior for the regression coefficients
  // nmc and nburn say how long to run the Markov chain

  GetRNGstate();
  RNGScope scope;

  // Precompute what we can
  int ncases = X.n_rows;
  int nfeatures = X.n_cols;
  colvec PriorPrecxMean = PriorPrecision * PriorMean;
  colvec pgscale(ncases, fill::ones);
  colvec onehalf(ncases);
  onehalf.fill(0.5);

  // Storage for draws
  mat betasave(nmc, nfeatures);
  colvec priorprobsave(ncases, fill::zeros);
  colvec postprobsave(ncases, fill::zeros);

  // NumericMatrix weightssave(nmc, ncomps);
  // NumericMatrix musave(nmc, ncomps);
  // NumericMatrix varsave(nmc, ncomps);

  // Initialize MCMC
  mat PosteriorPrecision(nfeatures, nfeatures);
  mat PosteriorVariance(nfeatures, nfeatures);
  colvec beta(nfeatures, fill::zeros);
  colvec omega(ncases);
  colvec PosteriorMean(nfeatures, fill::zeros);

  // // Mixture model parameters for alternative hypothesis
  // LogicalVector Gamma(ncases);
  // NumericVector comp_alpha(ncomps,1.0); // Dirichlet prior for mixing weifhts
  // NumericVector comp_means = muguess;
  // NumericVector comp_tau2(ncomps, 1.0);
  // NumericVector comp_marginalvariance = comp_tau2 + 1.0;
  // NumericVector comp_weights(ncomps, 1.0);
  // comp_weights = comp_weights/sum(comp_weights);
  // NumericVector M1 = dnorm(z, 0.0, 5.0);
  
  // Main MCMC loop
  int nsamples = nmc + nburn;
  for(int i=0; i<nsamples; i++) {

    colvec logodds = X*beta;
    colvec priorprob = InvLogit(logodds);
    colvec postprob = PriorToPostprobFullBayes(priorprob, M1, M0);
    
    // Draw Polya-Gamma latents, given current regression estimate of log odds
    omega = rpg(pgscale, logodds);

    // Update regression parameters
    PosteriorPrecision = X.t() * diagmat(omega) * X + PriorPrecision;
    PosteriorVariance = inv(PosteriorPrecision);
    PosteriorMean = PosteriorVariance*(X.t() * (postprob - onehalf) + PriorPrecxMean);
    beta = PosteriorMean + rmvnormArma(PosteriorVariance);

    // // Component-level sufficient statistics: just redeclare each pass through
    // NumericVector sumsignals(ncomps, 0.0);
    // NumericVector nsignals(ncomps, 0.0);
    // NumericVector totsumsq(ncomps, 0.0);
    // int whichcomponent;
    // double muhat, muvar, resid;

    // // Draw signal versus noise indicator
    // NumericVector udraws = Rcpp::runif(ncases);
    // //Gamma = (udraws <= Rcpp::as<NumericVector>(wrap(postprob)));
    // for(int j=0; j<ncases; ++j) {
    //   Gamma[i] = (udraws[i] < postprob(i)) ? TRUE : FALSE;
    // }
    
    // NumericVector Signals = subsetter(z, Gamma);
    // int totsignals = Signals.size();
    // NumericVector sigma2(totsignals, 1.0);
    // IntegerVector mycomps = FDRreg::draw_mixture_component(Signals, sigma2, comp_weights, comp_means, comp_tau2);

    // // Update sufficient statistics for each component
    // for(int j=0; j<Signals.size(); j++) {
    //   nsignals[mycomps[j]] += 1.0;
    //   sumsignals[mycomps[j]] += Signals[j];
    //   resid = (Signals[j] - comp_means[mycomps[j]]);
    //   totsumsq[mycomps[j]] += std::pow(resid, 2.0);
    // }

    // // Update component-level parameters using the sufficient stats we've aggregated
    // NumericVector eps = rnorm(ncomps);
    // for(int k=0; k < ncomps; k++) {
    //   comp_marginalvariance[k] = 1.0/rtgamma_once( (nsignals[k] + 1.0)/2.0, (totsumsq[k] + 1.0)/2.0, 0.0, 1.0);
    //   comp_tau2[k] =  comp_marginalvariance[k] - 1.0;
    //   muvar = comp_tau2[k]/(nsignals[k] + comp_tau2[k]*0.01);
    //   muhat = sumsignals[k]/(nsignals[k] + 0.01*comp_tau2[k]);
    //   comp_means[k] = muhat + sqrt(muvar)*eps[k];
    // }
    // comp_weights = rdirichlet_once(comp_alpha + nsignals);

    // // Recalculate marginal likelihood under mixture-of-normals alternative
    // M1 = FDRreg::dnormix(z, comp_weights, comp_means, comp_marginalvariance);

    if(i >= nburn) {
      betasave.row(i-nburn) = beta.t();
      priorprobsave = priorprobsave + priorprob;
      postprobsave = postprobsave + postprob;
    }
  }
  postprobsave = (1.0/nmc)*postprobsave;
  priorprobsave = (1.0/nmc)*priorprobsave;

  PutRNGstate();

  return Rcpp::List::create(Rcpp::Named("priorprob")=priorprobsave,
			    Rcpp::Named("postprob")=postprobsave,
			    Rcpp::Named("betasave")=betasave
			    // Rcpp::Named("musave")=musave,
			    // Rcpp::Named("weightssave")=weightssave,
			    // Rcpp::Named("varsave")=varsave
			    );
}


