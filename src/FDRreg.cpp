#include "RcppArmadillo.h"
#include <R_ext/Utils.h>
#include <FDRreg.h>
#include <iostream>
#include <exception>
#include "RNG.h"
#include "PolyaGamma.h"

// Rcpp::depends(RcppArmadillo)

using namespace Rcpp;
using namespace arma;

template <typename T>
T square( const T& x) {
  return x * x;
}

class interrupt_exception : public std::exception {
public:
    /**
     * Constructor.
     * @param[in] message A description of event that
     *  caused this exception.
     */
    interrupt_exception(std::string message)
	: detailed_message(message)
    {};

    /**
     * Virtual destructor. Needed to avoid "looser throw specification" errors.
     */
    virtual ~interrupt_exception() throw() {};

    /**
     * Obtain a description of the exception.
     * @return Description.
     */
    virtual const char* what() const throw() {
	return detailed_message.c_str();
    }

    /**
     * String with details on the error.
     */
    std::string detailed_message;
};


static inline void check_interrupt_impl(void* /*dummy*/) {
    R_CheckUserInterrupt();
}

inline bool check_interrupt() {
    return (R_ToplevelExec(check_interrupt_impl, NULL) == FALSE);
}

// simulate a single mean-zero multivariate normal vector
mat rmvnormArma(mat sigma) {
   int ncols = sigma.n_cols;
   mat Y = randn(ncols, 1);
   return chol(sigma) * Y;
}


// Draw from a truncated gamma distribution

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

// Draw from a Dirichlet distribution
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
colvec InvLogit(const arma::vec & x) {
  int d = x.n_elem;
  colvec result(d);
  for(int i=0; i<d; i++) {
    result[i] = 1.0/(1.0+exp(-x(i)));
  }
  return result;
}

// priorprob = vector of prior probabilities
// mfull = vector of marginal likelihoods under global predictive (mixture of null + alternative)
// mnull = vector of likelihoods under the null predictive
colvec PriorToPostprobEBayes(const arma::vec & priorprob, NumericVector mfull, NumericVector mnull) {
  int n = priorprob.n_elem;
  colvec postprob(n);
  for(int i=0; i<n; i++) {
    if(mnull[i] >= mfull[i]) {
      postprob(i) = 0.0;
    } else {
      postprob(i) =  1.0-(1.0-priorprob(i))*(mnull[i]/mfull[i]);
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
    f1 = priorprob[i] * m1[i];
    f0 = (1.0-priorprob[i]) * m0[i];
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
SEXP FDRregCPP(NumericVector z, const arma::mat & X, NumericVector M0, NumericVector MTot, const arma::mat & PriorPrecision, const arma::vec & PriorMean, int nmc, int nburn) {
  // z is an n-vector of test statistics
  // X is a design matrix in the FDRR problem, including the intercept
  // M0 is a vector of marginal likelihoods, f_0(z), under the null hypothesis
  // MTot is a vector of marginal likelihoods, f(z), under the mixture of null and alternative
  // PriorPrecision and PriorMean are the parameters of the multivariate normal prior for the regression coefficients
  // nmc and nburn say how long to run the Markov chain

  // Precompute what we can
  int ncases = X.n_rows;
  int nfeatures = X.n_cols;
  colvec PriorPrecxMean = PriorPrecision * PriorMean;
  colvec pgscale(ncases, fill::ones);
  colvec onehalf(ncases);
  onehalf.fill(0.5);
  bool interrupt = false;

  // Storage for draws
  mat betasave(nmc, nfeatures);
  colvec priorprobsave(ncases, fill::zeros);
  colvec postprobsave(ncases, fill::zeros);

  // Initialize MCMC
  int nsamples = nmc + nburn;
  mat PosteriorPrecision(nfeatures, nfeatures);
  mat PosteriorVariance(nfeatures, nfeatures);
  colvec beta(nfeatures, fill::zeros);
  colvec logodds(ncases, fill::zeros);
  colvec priorprob(ncases);
  colvec postprob(ncases, fill::randu);
  colvec omega(ncases);
  colvec PosteriorMean(nfeatures, fill::zeros);
  
  // Initialize intercept at something reasonable
  beta[0] = 0.0;
  for(int i=0; i<ncases; i++) {
    if(fabs(z[i])  > 2.0) {
      beta[0] += 1.0;
    }
  }
  beta[0] = beta[0]/ncases;
  beta[0] = log(beta[0]/(1.0-beta[0]));

  // Main MCMC loop
  for(int i=0; i<nsamples; i++) {

    // Check for R user interrupt
    if(i % 100 == 0) {
      if (check_interrupt()) {
	interrupt = true;
      }
      // throw exception if interrupt occurred
      if (interrupt) {
	throw interrupt_exception("The FDRR model fit was interrupted.");
      }
    }

    // compute prior and posterior probability of alternative hypothesis
    logodds = X*beta;
    priorprob = InvLogit(logodds);
    postprob = PriorToPostprobEBayes(priorprob, MTot, M0);
    
    // Draw Polya-Gamma latents, given current regression estimate of log odds
    omega = rpg(pgscale, logodds);

    // Update regression parameters
    PosteriorPrecision = X.t() * diagmat(omega) * X + PriorPrecision;
    PosteriorVariance = inv(PosteriorPrecision);
    PosteriorMean = PosteriorVariance*(X.t() * (postprob - onehalf) + PriorPrecxMean);
    beta = PosteriorMean + rmvnormArma(PosteriorVariance);
    if(i >= nburn) {
      betasave.row(i-nburn) = PosteriorMean.t();
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
SEXP BayesFDRregCPP(NumericVector z, const arma::mat & X, NumericVector M0, int ncomps, const arma::mat & PriorPrecision, const arma::vec & PriorMean, int nmc, int nburn, NumericVector grid) {
  // z is an n-vector of test statistics
  // X is a design matrix in the FDRR problem, including the intercept
  // M0 is a vector of marginal likelihoods, f_0(z), under the null hypothesis
  // ncomps the number of components you want in the Gaussian mixture model for the alternative
  // PriorPrecision and PriorMean are the parameters of the multivariate normal prior for the regression coefficients
  // nmc and nburn say how long to run the Markov chain

  // Precompute what we can
  int ncases = X.n_rows;
  int nfeatures = X.n_cols;
  colvec PriorPrecxMean = PriorPrecision * PriorMean;
  colvec pgscale(ncases, fill::ones);
  colvec onehalf(ncases);
  NumericVector sigma2(ncases, 1.0);
  onehalf.fill(0.5);
  bool interrupt = false;

  // Storage for draws
  mat betasave(nmc, nfeatures);
  colvec priorprobsave(ncases, fill::zeros);
  colvec postprobsave(ncases, fill::zeros);
  colvec fthetasave(grid.size(), fill::zeros);
  NumericMatrix weightssave(nmc, ncomps);
  NumericMatrix musave(nmc, ncomps);
  NumericMatrix varsave(nmc, ncomps);

  // Initialize MCMC
  int nsamples = nmc + nburn;
  mat PosteriorPrecision(nfeatures, nfeatures);
  mat PosteriorVariance(nfeatures, nfeatures);
  colvec beta(nfeatures, fill::zeros);
  colvec logodds(ncases, fill::zeros);
  colvec priorprob(ncases);
  colvec postprob(ncases, fill::randu);
  colvec omega(ncases);
  colvec PosteriorMean(nfeatures, fill::zeros);
  IntegerVector Gamma(ncases);

  // Mixture model parameters for alternative hypothesis
  NumericVector comp_alpha(ncomps,1.0); // Dirichlet prior for mixing weifhts
  NumericVector comp_means(ncomps, 0.0);
  NumericVector comp_marginalvariance(ncomps, 4.0); // = tau2[k] + 1, integrating out the random mean
  NumericVector comp_marginalsd(ncomps, 2.0);
  NumericVector comp_weights(ncomps, 1.0/ncomps);
  NumericVector M1 = dnorm(z,0.0,3.0);

  // Initialize intercept at something reasonable
  beta[0] = 0.0;
  for(int i=0; i<ncases; i++) {
    if( fabs(z[i])  > 2.0 ) {
      beta[0] += 1.0;
    }
  }
  beta[0] = beta[0]/ncases;
  beta[0] = log(beta[0]/(1-beta[0]));

  // Initialize component means
  comp_means = rnorm(5, 0.0, 3.0);

  // Main MCMC loop
  for(int i=0; i<nsamples; i++) {

    // Check for R user interrupt
    if(i % 100 == 0) {
      if (check_interrupt()) {
	interrupt = true;
      }
      // throw exception if interrupt occurred
      if (interrupt) {
	throw interrupt_exception("The FDRR model fit was interrupted.");
      }
    }


    // compute prior and posterior probability of alternative hypothesis
    logodds = X*beta;
    priorprob = InvLogit(logodds);
    postprob = PriorToPostprobFullBayes(priorprob, M1, M0);
    
    // Update mixture model for alternative hypothesis
    for(int j=0; j<ncases; j++) {
      Gamma[j] = Rf_rbinom(1.0, postprob[j]);
    }

    // Component-level sufficient statistics: just redeclare each pass through
    NumericVector thisprob(ncomps, 0.0);
    NumericVector sumsignals(ncomps, 0.0);
    NumericVector nsignals(ncomps,0.0);
    NumericVector totsumsq(ncomps, 0.0);
    NumericVector tval(ncomps, 0.0);
    int whichcomponent=ncomps;
    double muhat, muvar, tau2, resid;

    // for the signals, draw the mixture component and collect the sufficient statistics
    // Note that we don't actually need to save the model indicators
    for(int j=0; j < ncases; j++) {
      if(Gamma[j] != 0) {
	tval = (comp_means - z[j])/comp_marginalsd;
	thisprob = comp_weights*(dnorm(tval,0.0,1.0)/comp_marginalsd);
	whichcomponent = FDRreg::mysample(thisprob);
	nsignals[whichcomponent] += 1.0;
	sumsignals[whichcomponent] += z[j];
	resid = z[j] - comp_means[whichcomponent];
	totsumsq[whichcomponent] += resid*resid;
      }
    }

    // Update component-level parameters using the sufficient stats we've aggregated
    NumericVector eps = rnorm(ncomps, 0.0, 1.0);
    for(int k=0; k<ncomps; k++) {
      comp_marginalvariance[k] = 1.0/rtgamma_once( (1.0+sumsignals[k])/2.0, (1.0 + totsumsq[k])/2.0, 0.0, 1.0);
      comp_marginalsd[k] = sqrt(comp_marginalvariance[k]);
      tau2 =  comp_marginalvariance[k] - 1.0;
      muvar = tau2/(nsignals[k] + tau2*0.01);
      muhat = sumsignals[k]/(nsignals[k] + 0.01*tau2);
      comp_means[k] = muhat + sqrt(muvar)*eps[k];
    }
    comp_weights = rdirichlet_once(comp_alpha + nsignals);

    // // Recalculate marginal likelihood under mixture-of-normals alternative
    // for(int j=0; j < ncases; j++) {
    //   eps = (comp_means - z[j])/sqrt(comp_marginalvariance);
    //   M1[j] = sum(comp_weights * dnorm(eps, 0.0, 1.0)/sqrt(comp_marginalvariance));
    // }

    // Draw Polya-Gamma latents, given current regression estimate of log odds
    omega = rpg(pgscale, logodds);

    // Update regression parameters
    PosteriorPrecision = X.t() * diagmat(omega) * X + PriorPrecision;
    PosteriorVariance = inv(PosteriorPrecision);
    PosteriorMean = PosteriorVariance*(X.t() * (postprob - onehalf) + PriorPrecxMean);
    beta = PosteriorMean + rmvnormArma(PosteriorVariance);
    if(i >= nburn) {
      betasave.row(i-nburn) = PosteriorMean.t();
      priorprobsave = priorprobsave + priorprob;
      postprobsave = postprobsave + postprob;
      weightssave(i-nburn,_) = nsignals;
      musave(i-nburn,_) = sumsignals;
      varsave(i-nburn,_) = totsumsq;
    }
  }
  postprobsave = (1.0/nmc)*postprobsave;
  priorprobsave = (1.0/nmc)*priorprobsave;

  return Rcpp::List::create(Rcpp::Named("priorprob")=priorprobsave,
			    Rcpp::Named("postprob")=postprobsave,
			    Rcpp::Named("betasave")=betasave,
			    Rcpp::Named("musave")=musave,
			    Rcpp::Named("weightssave")=weightssave,
			    Rcpp::Named("varsave")=varsave
			    );
}




