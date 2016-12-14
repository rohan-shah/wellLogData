#ifndef BASIC_FILTER_HEADER_GUARD
#define BASIC_FILTER_HEADER_GUARD
#include "context.h"
#include <boost/random/mersenne_twister.hpp>
namespace wellLogData
{
	struct basicFilterArgs
	{
	public:
		basicFilterArgs(boost::mt19937& randomSource)
			: randomSource(randomSource)
		{}
		boost::mt19937& randomSource;
		std::size_t nParticles;
		std::vector<double> outlierProbabilities, changeProbabilities;
	};
	void basicFilter(const context& contextObj, basicFilterArgs& args);
}
#endif
