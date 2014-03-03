#include "RcppArmadillo.h"
#include <R_ext/Utils.h>
#include <iostream>
#include <exception>
#include "RNG.h"
#include "PolyaGamma.h"
// Rcpp::depends(RcppArmadillo)

using namespace Rcpp;
using namespace arma;


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

// compute the inverse logit transform for an arma::vec
colvec InvLogit(const arma::vec & x) {
  int d = x.n_elem;
  colvec result(d);
  for(int i=0; i<d; i++) {
    result[i] = 1.0/(1.0+exp(-x(i)));
  }
  return result;
}

// multiplies each row of X by the corresponding entry in lambda
// equivalent to diag(lambda) %*% X
// just checking to see if this led to a speedup vis-a-vis native diagmat(lambda) * X
mat ProdOnRow(const arma::mat & X, const arma::vec & lambda) {
  int n = X.n_rows;
  int d = X.n_cols;
  mat result(n,d);
  for(int i=0; i<n; i++) {
    result.row(i) = lambda(i)*X.row(i);
  }
  return result;
}

// compute the inverse logit transform for an arma::vec
// mfull = vector of marginal likelihoods under global predictive (mixture of null + alternative)
// mnull = vector of likelihoods under the null predictive
colvec PriorToPostprob(const arma::vec & priorprob, NumericVector mfull, NumericVector mnull) {
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

// Draws random PG variances from arma::vectors of n's and psi's
colvec rpg(colvec shape, colvec scale) {
  // Set up the RNG
  RNG r;
  PolyaGamma pg;
#ifdef USE_R
  GetRNGstate();
#endif
  int d = shape.n_elem;
  colvec result(d);
  for(int i=0; i<d; i++) {
    result[i] = pg.draw(shape(i), scale(i), r);
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
    postprob = PriorToPostprob(priorprob, MTot, M0);
    
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

