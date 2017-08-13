#ifndef WITHOUT_REPLACEMENT_WITH_VARIANCE_HEADER_GUARD
#define WITHOUT_REPLACEMENT_WITH_VARIANCE_HEADER_GUARD
#include "context.h"
#include <boost/random/mersenne_twister.hpp>
#include "conditionalPoissonSequential.h"
namespace wellLogData
{
	struct withoutReplacementWithVarianceArgs
	{
	public:
		withoutReplacementWithVarianceArgs(boost::mt19937& randomSource)
			: randomSource(randomSource)
		{}
		boost::mt19937& randomSource;
		std::size_t nParticles;
		std::vector<double> outlierProbabilities, changeProbabilities;
		std::vector<double> changeEstimateNumeratorVariances;
		double normalisingConstant;
	};
	void withoutReplacementWithVariance(const context& contextObj, withoutReplacementWithVarianceArgs& args);
}
#endif
