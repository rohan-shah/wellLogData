#ifndef FEARNHEAD_FILTER_HEADER_GUARD
#define FEARNHEAD_FILTER_HEADER_GUARD
#include "context.h"
#include <boost/random/mersenne_twister.hpp>
namespace wellLogData
{
	struct fearnheadFilterArgs
	{
	public:
		fearnheadFilterArgs(boost::mt19937& randomSource)
			: randomSource(randomSource)
		{}
		boost::mt19937& randomSource;
		std::size_t nParticles;
		std::vector<double> outlierProbabilities, changeProbabilities;
	};
	void fearnheadFilter(const context& contextObj, fearnheadFilterArgs& args);
}
#endif
