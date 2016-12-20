#ifndef WITHOUT_REPLACEMENT_HEADER_GUARD
#define WITHOUT_REPLACEMENT_HEADER_GUARD
#include "context.h"
#include <boost/random/mersenne_twister.hpp>
namespace wellLogData
{
	struct withoutReplacementArgs
	{
	public:
		withoutReplacementArgs(boost::mt19937& randomSource)
			: randomSource(randomSource)
		{}
		boost::mt19937& randomSource;
		std::size_t nParticles;
		std::vector<double> outlierProbabilities, changeProbabilities;
	};
	void withoutReplacement(const context& contextObj, withoutReplacementArgs& args);
}
#endif
