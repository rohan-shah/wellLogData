#include "withoutReplacementWithVarianceRPackage.h"
#include "withoutReplacementWithVariance.h"
SEXP withoutReplacementWithVariance(SEXP data_sexp, SEXP mu_sexp, SEXP sigma_sexp, SEXP nu_sexp, SEXP tau1_sexp, SEXP tau2_sexp, SEXP changeProbability_sexp, SEXP outlierProbability_sexp, SEXP outlierClusterProbability_sexp, SEXP nParticles_sexp, SEXP seed_sexp)
{
BEGIN_RCPP
	std::vector<double> data;
	try
	{
		data = Rcpp::as<std::vector<double> >(data_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input data must be a numeric vector");
	}
	if(data.size() <= 1)
	{
		throw std::runtime_error("Input data must have at least two values");
	}
	
	double mu;
	try
	{
		mu = Rcpp::as<double>(mu_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input mu must be a number");
	}
	if(mu != mu) throw std::runtime_error("Input mu cannot be NA");

	double sigma;
	try
	{
		sigma = Rcpp::as<double>(sigma_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input sigma must be a number");
	}
	if(sigma != sigma) throw std::runtime_error("Input sigma cannot be NA");
	if(sigma <= 0) throw std::runtime_error("Input sigma must be positive");

	double nu;
	try
	{
		nu = Rcpp::as<double>(nu_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input nu must be a number");
	}
	if(nu != nu) throw std::runtime_error("Input nu cannot be NA");
	
	double tau1;
	try
	{
		tau1 = Rcpp::as<double>(tau1_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input tau1 must be a number");
	}
	if(tau1 != tau1) throw std::runtime_error("Input tau1 cannot be NA");
	if(tau1 <= 0) throw std::runtime_error("Input tau1 must be positive");

	double tau2;
	try
	{
		tau2 = Rcpp::as<double>(tau2_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input tau2 must be a number");
	}
	if(tau2 != tau2) throw std::runtime_error("Input tau2 cannot be NA");
	if(tau2 <= 0) throw std::runtime_error("Input tau2 must be positive");

	double changeProbability;
	try
	{
		changeProbability = Rcpp::as<double>(changeProbability_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input changeProbability must be a number");
	}
	if(changeProbability != changeProbability) throw std::runtime_error("Input changeProbability cannot be NA");
	if(changeProbability <= 0 || changeProbability >= 1) throw std::runtime_error("Input changeProbability must be in (0, 1)");

	double outlierProbability;
	try
	{
		outlierProbability = Rcpp::as<double>(outlierProbability_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input outlierProbability must be a number");
	}
	if(outlierProbability != outlierProbability) throw std::runtime_error("Input outlierProbability cannot be NA");
	if(outlierProbability <= 0 || outlierProbability >= 1) throw std::runtime_error("Input outlierProbability must be in (0, 1)");
	
	double outlierClusterProbability;
	try
	{
		outlierClusterProbability = Rcpp::as<double>(outlierClusterProbability_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input outlierClusterProbability must be a number");
	}
	if(outlierClusterProbability != outlierClusterProbability) throw std::runtime_error("Input outlierClusterProbability cannot be NA");
	if(outlierClusterProbability <= 0 || outlierClusterProbability >= 1) throw std::runtime_error("Input outlierClusterProbability must be in (0, 1)");

	int nParticles_int;
	try
	{
		nParticles_int = Rcpp::as<int>(nParticles_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input outlierClusterProbability must be an integer");
	}
	if(nParticles_int <= 1)
	{
		throw std::runtime_error("Input outlierClusterProbability must be positive");
	}
	std::size_t nParticles = (std::size_t)nParticles_int;

	int seed;
	try
	{
		seed = Rcpp::as<int>(seed_sexp);
	}
	catch(...)
	{
		throw std::runtime_error("Input seed must be an integer");
	}
	boost::mt19937 randomSource;
	randomSource.seed(seed);

	wellLogData::contextArgs args;
	args.mu = mu;
	args.sigma = sigma;
	args.tau1 = tau1;
	args.tau2 = tau2;
	args.nu = nu;
	args.changeProb = changeProbability;
	args.outlierProb = outlierProbability;
	args.outlierClusterProb = outlierClusterProbability;

	wellLogData::context contextObj(args, data);

	wellLogData::withoutReplacementWithVarianceArgs methodArgs(randomSource);
	methodArgs.nParticles = nParticles;
	wellLogData::withoutReplacementWithVariance(contextObj, methodArgs);
	return Rcpp::List::create(Rcpp::Named("outlierProbabilities") = Rcpp::wrap(methodArgs.outlierProbabilities), Rcpp::Named("changeProbabilities") = Rcpp::wrap(methodArgs.changeProbabilities));
END_RCPP
}
