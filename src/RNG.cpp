#include "RNG.h"

// #ifdef USE_R
// #include "RRNG.cpp"
// #else
// #include "GRNG.cpp"
// #endif

// Throw runtime exception or return.
#ifndef TREOR
#ifndef NTHROW
#define TREOR(MESS, VAL) throw std::runtime_error(MESS);
#else
#define TREOR(MESS, VAL) Rprintf( MESS); return VAL;
#endif
#endif

//////////////////////////////////////////////////////////////////////
	       // TRUNCATED NORMAL HELPER FUNCTIONS //
//////////////////////////////////////////////////////////////////////

double RNG::alphastar(double left)
{
  return 0.5 * (left + sqrt(left + 4));
} // alphastar

double RNG::lowerbound(double left)
{
  double astar  = alphastar(left);
  double lbound = left + exp(0.5 * left * (left - astar) + 0.5);
  return lbound;
} // lowerbound

//////////////////////////////////////////////////////////////////////
		     // DRAW TRUNCATED NORMAL //
//////////////////////////////////////////////////////////////////////

double RNG::tnorm(double left)
{
  double rho, ppsl;

  if (left < 0) { // Accept/Reject Normal
    while (true) {
      ppsl = norm(0.0, 1.0);
      if (ppsl > left) return ppsl;
    }
  }
  else { // Accept/Reject Exponential
    double astar = alphastar(left);
    while (true) {
      ppsl = expon_rate(astar) + left;
      rho  = exp( -0.5 * (ppsl - astar) * (ppsl - astar) );
      if (unif() < rho) return ppsl;
    }
  }
} // tnorm
//--------------------------------------------------------------------

double RNG::tnorm(double left, double right)
{
  // The most difficult part of this algorithm is figuring out all the
  // various cases.  An outline is summarized in the Appendix.

  // Check input 
  #ifdef USE_R
  if (ISNAN(right) || ISNAN(left))
  #else
  if (std::isnan(right) || std::isnan(left)) 
  #endif
  {
    Rprintf( "Warning: nan sent to RNG::tnorm: left=%g, right=%g\n", left, right);
    TREOR("RNG::tnorm: parameter problem.\n", 0);
    // throw std::runtime_error("RNG::tnorm: parameter problem.\n");
  }

  if (right < left) {
    Rprintf( "Warning: left: %g, right:%g.\n", left, right);
    TREOR("RNG::tnorm: parameter problem.\n", 0);
  }

  double rho, ppsl;

  if (left >= 0) {
    double lbound = lowerbound(left);
    if (right > lbound) { // Truncated Exponential.
      double astar = alphastar(left);
      while (true) {
	do
	  ppsl = expon_rate(astar) + left;
	while(ppsl > right);
	// REVIEW REVIEW - TAKE ANOTHER LOOK AT THIS.
	rho  = exp(-0.5*(ppsl - astar)*(ppsl-astar));
	if (unif() < rho) return ppsl;
	// if (ppsl < right) return ppsl;
      }
    }
    else {
      while (true) {
	ppsl = flat(left, right);
	rho  = exp(0.5 * (left*left - ppsl*ppsl));
	if (unif() < rho) return ppsl;
      }
    }
  }
  else if (right >= 0) {
    if ( (right - left) < SQRT2PI ){
      while (true) {
	ppsl = flat(left, right);
	rho  = exp(-0.5 * ppsl * ppsl);
	if (unif() < rho) return ppsl;
      }
    }
    else{
      while (true) {
	ppsl = norm(0, 1);
	if (left < ppsl && ppsl < right) return ppsl;
      }
    }
  }
  else {
    return -1. * tnorm(-1.0 * right, -1.0 * left);
  }
} // tnorm
//--------------------------------------------------------------------

double RNG::tnorm(double left, double mu, double sd)
{
  double newleft = (left - mu) / sd;
  return mu + tnorm(newleft) * sd;
} // tnorm
//--------------------------------------------------------------------

double RNG::tnorm(double left, double right, double mu, double sd)
{
  double newleft  = (left - mu) / sd;
  double newright = (right - mu) / sd;

  // I want to check this here as well so we can see what the input was.
  // It may be more elegant to try and catch tdraw.
  if (newright < newleft) {
    Rprintf( "left, right, mu, sd: %g, %g, %g, %g \n", left, right, mu, sd);
    Rprintf( "nleft, nright: %g, %g\n", newleft, newright);
    TREOR("RNG::tnorm: parameter problem.\n", 0);
  }

  double tdraw = tnorm(newleft, newright);
  double draw = mu + tdraw * sd;

  // It may be the case that there is some numerical error and that the draw
  // ends up out of bounds.
  if (draw < left || draw > right){
    Rprintf( "Error in tnorm: draw not in bounds.\n");
    Rprintf( "left, right, mu, sd: %g, %g, %g, %g\n", left, right, mu, sd);
    Rprintf( "nleft, nright, tdraw, draw: %g, %g, %g, %g\n", newleft, newright, tdraw, draw);
  }

  return draw;
} // tnorm
//--------------------------------------------------------------------

// Right tail of normal by Devroye
//------------------------------------------------------------------------------
double RNG::tnorm_tail(double t)
{
  double E1 = expon_rate(1.0);
  double E2 = expon_rate(1.0);
  while ( E1*E1 > 2 * E2 / t) {
    E1 = expon_rate(1.0);
    E2 = expon_rate(1.0);
  }
  return (1 + t * E1) / sqrt(t);
}

//------------------------------------------------------------------------------

// Truncatatoin at t = 1.
inline double RNG::right_tgamma_reject(double shape, double rate)
{
  double x = 2.0;
  while (x > 1.0) 
    x = gamma_rate(shape, rate);
  return x;
}

double RNG::omega_k(int k, double a, double b)
{
  double log_coef = -b + (a+k-1) * log(b) - Gamma(a+k, true) - p_gamma_rate(1.0, a, b, true);
  return exp(log_coef);
}

// Truncation at t = 1.
double RNG::right_tgamma_beta(double shape, double rate)
{
  double a = shape;
  double b = rate;

  double u = unif();

  int k = 1;
  double cdf = omega_k(1, a, b);
  while (u > cdf) {
    cdf += omega_k(++k, a, b);
    if (k % 100000 == 0) {
      Rprintf("right_tgamma_beta (itr k=%i): a=%g, b=%g, u=%g, cdf=%g\n", k, a, b, u, cdf);
      #ifdef USE_R
      R_CheckUserInterrupt();
      #endif
    }
  }

  return beta(a, k);
}

double RNG::rtgamma_rate(double shape, double rate, double right_t)
{
  // x \sim (a,b,t)
  // ty = x
  // y \sim (a, bt, 1);
  double a = shape;
  double b = rate * right_t;
  
  double p = p_gamma_rate(1, a, b);
  double y = 0.0;
  if (p > 0.95)
    y = right_tgamma_reject(a, b);
  else
    y = right_tgamma_beta(a,b);
  
  double x = right_t * y;
  return x;
}

//------------------------------------------------------------------------------
double RNG::ltgamma(double shape, double rate, double trunc)
{
  double a = shape;
  double b = rate * trunc;

  if (trunc <=0) {
    Rprintf( "ltgamma: trunc = %g < 0\n", trunc);
    return 0;
  }
  if (shape < 1) {
    Rprintf( "ltgamma: shape = %g < 1\n", shape);
    return 0;
  }

  if (shape ==1) return expon_rate(1) / rate + trunc;

  double d1 = b-a;
  double d3 = a-1;
  double c0 = 0.5 * (d1 + sqrt(d1*d1 + 4 * b)) / b;
   
  double x = 0.0;
  bool accept = false;

  while (!accept) {
    x = b + expon_rate(1) / c0;
    double u = unif();
    
    double l_rho = d3 * log(x) - x * (1-c0);
    double l_M   = d3 * log(d3 / (1-c0)) - d3;

    accept = log(u) <= (l_rho - l_M);
  }

  return trunc * (x/b);
}

//------------------------------------------------------------------------------
double RNG::igauss(double mu, double lambda)
{
  // See R code for specifics.
  double mu2 = mu * mu;
  double Y = norm(0.0, 1.0);
  Y *= Y;
  double W = mu + 0.5 * mu2 * Y / lambda;
  double X = W - sqrt(W*W - mu2);
  if (unif() > mu / (mu + X)) 
    X = mu2 / X;
  return X;
}

//------------------------------------------------------------------------------
double RNG::rtinvchi2(double scale, double trunc)
{
  double R = trunc / scale;
  // double X = 0.0;
  // // I need to consider using a different truncated normal sampler.
  // double E1 = r.expon_rate(1.0); double E2 = r.expon_rate(1.0);
  // while ( (E1*E1) > (2 * E2 / R)) {
  //   // Rprintf("E %g %g %g %g\n", E1, E2, E1*E1, 2*E2/R);
  //   E1 = r.expon_rate(1.0); E2 = r.expon_rate(1.0);
  // }
  // // Rprintf("E %g %g \n", E1, E2);
  // X = 1 + E1 * R;
  // X = R / (X * X);
  // X = scale * X;
  double E = tnorm(1/sqrt(R));
  double X = scale / (E*E);
  return X;
}

//------------------------------------------------------------------------------

double RNG::Beta(double a, double b, bool log)
{
  double out = Gamma(a, true) + Gamma(b, true) - Gamma(a+b,true);
  if (!log) out = exp(out);
  return out;
}

//------------------------------------------------------------------------------

double RNG::p_igauss(double x, double mu, double lambda)
{
  // z = 1 / mean
  double z = 1 / mu;
  double b = sqrt(lambda / x) * (x * z - 1);
  double a = sqrt(lambda / x) * (x * z + 1) * -1.0;
  double y = RNG::p_norm(b) + exp(2 * lambda * z) * RNG::p_norm(a);
  return y;
}
