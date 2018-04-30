#ifndef CONTEXT_HEADER_GUARD
#define CONTEXT_HEADER_GUARD
#include <boost/noncopyable.hpp>
#include <limits>
#include <vector>
#define DOUBLE_NAN std::numeric_limits<double>::quiet_NaN()
namespace wellLogData
{
	struct contextArgs
	{
	public:
		contextArgs()
			:mu(DOUBLE_NAN), sigma(DOUBLE_NAN), tau1(DOUBLE_NAN), tau2(DOUBLE_NAN), nu(DOUBLE_NAN), changeProb(DOUBLE_NAN), outlierProb(DOUBLE_NAN), outlierClusterProb(DOUBLE_NAN)
		{}
		double mu, sigma, tau1, tau2, nu;
		double changeProb, outlierProb, outlierClusterProb;
	};
	class context : public boost::noncopyable
	{
	public:
		context(const contextArgs&, const std::vector<double>& data);
		double getChangeProbability() const;
		//This is the probability of having an outlier
		double getOutlierProbability() const;
		//This is the probability of having consecutive outliers, so that they form a cluster. 
		double getOutlierClusterProbability() const;
		const std::vector<double>& getData() const;
		double getMu() const;
		double getNu() const;
		double getSigma() const;
		double getSigmaSquared() const;
		double getTau1Squared() const;
		double getTau2Squared() const;
		double getTau2() const;
	private:
		double mu, sigma, sigmaSquared, tau1Squared, tau2, tau2Squared, nu;
		double changeProb, outlierProb, outlierClusterProb;
		std::vector<double> data;
		context();
		context(const context&);
	};
}
#undef DOUBLE_NAN
#endif
