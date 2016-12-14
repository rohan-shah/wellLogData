#include "context.h"
#include <stdexcept>
#define DOUBLE_NAN std::numeric_limits<double>::quiet_NaN()
namespace wellLogData
{
	context::context(const contextArgs& args, const std::vector<double>& data)
		: mu(args.mu), sigma(DOUBLE_NAN), sigmaSquared(DOUBLE_NAN), tau1Squared(DOUBLE_NAN), tau2(DOUBLE_NAN), tau2Squared(DOUBLE_NAN), changeProb(args.changeProb), outlierProb(args.outlierProb), outlierClusterProb(args.outlierClusterProb), data(data)
	{
		const double infty = std::numeric_limits<double>::infinity();
		if(mu != mu || mu == infty || mu == -infty)
		{
			throw std::runtime_error("Input mu must be a number");
		}
		if(args.sigma != args.sigma || args.sigma <= 0 || args.sigma == infty)
		{
			throw std::runtime_error("Input sigma must be a positive number");
		}
		sigma = args.sigma;
		sigmaSquared = args.sigma * args.sigma;
		if(args.tau1 != args.tau1 || args.tau1 <= 0 || args.tau1 == infty)
		{
			throw std::runtime_error("Input tau1 must be a positive number");
		}
		tau1Squared = args.tau1 * args.tau1;
		if(args.tau2 != args.tau2 || args.tau2 <= 0 || args.tau2 == infty)
		{
			throw std::runtime_error("Input tau2 must be a positive number");
		}
		tau2 = args.tau2;
		tau2Squared = args.tau2 * args.tau2;
		if(nu != nu || nu == infty || nu == -infty)
		{
			throw std::runtime_error("Input nu must be a number");
		}
		if(changeProb != changeProb || changeProb <= 0 || changeProb >= 1)
		{
			throw std::runtime_error("Input changeProb must be a probability");
		}
		if(outlierProb != outlierProb || outlierProb <= 0 || outlierProb >= 1)
		{
			throw std::runtime_error("Input outlierProb must be a probability");
		}
		if(outlierClusterProb != outlierClusterProb || outlierClusterProb <= 0 || outlierClusterProb >= 1)
		{
			throw std::runtime_error("Input outlierClusterProb must be a probability");
		}
		if(data.size() <= 1)
		{
			throw std::runtime_error("There must be at least two data values");
		}
	}
	double context::getNu() const
	{
		return mu;
	}
	double context::getOutlierProbability() const
	{
		return outlierProb;
	}
	double context::getOutlierClusterProbability() const
	{
		return outlierClusterProb;
	}
	double context::getMu() const
	{
		return mu;
	}
	double context::getSigma() const
	{
		return sigma;
	}
	double context::getSigmaSquared() const
	{
		return sigmaSquared;
	}
	double context::getTau1Squared() const
	{
		return tau1Squared;
	}
	double context::getTau2() const
	{
		return tau2;
	}
	double context::getTau2Squared() const
	{
		return tau2Squared;
	}
	const std::vector<double>& context::getData() const
	{
		return data;
	}
}
#undef DOUBLE_NAN
