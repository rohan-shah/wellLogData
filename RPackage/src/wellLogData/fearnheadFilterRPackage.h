#ifndef FEARNHEAD_FILTER_RPACKAGE_HEADER_GUARD
#define FEARNHEAD_FILTER_RPACKAGE_HEADER_GUARD
#include <Rcpp.h>
SEXP fearnheadFilter(SEXP data, SEXP mu, SEXP sigma, SEXP nu, SEXP tau1, SEXP tau2, SEXP changeProbability, SEXP outlierProbability, SEXP outlierClusterProbability, SEXP nParticles, SEXP seed);
#endif
