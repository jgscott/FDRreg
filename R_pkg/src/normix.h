#ifndef __NORMIX__
#define __NORMIX__



using namespace Rcpp;

NumericVector dnormix(NumericVector y, NumericVector weights,  NumericVector mu, NumericVector tau2);
int mysample(NumericVector probs);
NumericVector marnormix(NumericVector y, NumericVector sigma2, NumericVector weights,  NumericVector mu, NumericVector tau2);
NumericVector rnormix(int n, NumericVector weights,  NumericVector mu, NumericVector tau2);
IntegerVector draw_mixture_component(NumericVector y, NumericVector sigma2, NumericVector weights,  NumericVector mu, NumericVector tau2);
inline double flag(double a, bool b);
NumericVector subsetter(NumericVector x, LogicalVector b);


#endif

