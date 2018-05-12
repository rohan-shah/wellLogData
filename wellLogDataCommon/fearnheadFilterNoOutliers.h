#ifndef FEARNHEAD_FILTER_NO_OUTLIER_HEADER_GUARD
#define FEARNHEAD_FILTER_NO_OUTLIER_HEADER_GUARD
#include "context.h"
#include <boost/random/mersenne_twister.hpp>
namespace wellLogData
{
	struct fearnheadFilterNoOutliersArgs
	{
	public:
		fearnheadFilterNoOutliersArgs(boost::mt19937& randomSource)
			: randomSource(randomSource)
		{}
		boost::mt19937& randomSource;
		std::size_t nParticles;
		std::vector<double> changeProbabilities;
	};
	void fearnheadFilterNoOutliers(const context& contextObj, fearnheadFilterNoOutliersArgs& args);
}
#endif
