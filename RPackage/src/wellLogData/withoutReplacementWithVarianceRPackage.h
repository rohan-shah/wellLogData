#ifndef WITHOUT_REPLACEMENT_WITH_VARIANCE_RPACKAGE_HEADER_GUARD
#define WITHOUT_REPLACEMENT_WITH_VARIANCE_RPACKAGE_HEADER_GUARD
#include <Rcpp.h>
SEXP withoutReplacementWithVariance(SEXP data, SEXP mu, SEXP sigma, SEXP nu, SEXP tau1, SEXP tau2, SEXP changeProbability, SEXP outlierProbability, SEXP outlierClusterProbability, SEXP nParticles, SEXP seed);
#endif
