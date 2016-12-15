#ifndef BASIC_FILTER_RPACKAGE_HEADER_GUARD
#define BASIC_FILTER_RPACKAGE_HEADER_GUARD
#include <Rcpp.h>
SEXP basicFilter(SEXP data, SEXP mu, SEXP sigma, SEXP nu, SEXP tau1, SEXP tau2, SEXP changeProbability, SEXP outlierProbability, SEXP outlierClusterProbability, SEXP nParticles, SEXP seed);
#endif
