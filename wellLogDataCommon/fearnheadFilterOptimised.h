#ifndef FEARNHEAD_FILTER_OPTIMISED_HEADER_GUARD
#define FEARNHEAD_FILTER_OPTIMISED_HEADER_GUARD
#include "context.h"
#include <boost/random/mersenne_twister.hpp>
#include "fearnheadFilter.h"
namespace wellLogData
{
	void fearnheadFilterOptimised(const context& contextObj, fearnheadFilterArgs& args);
}
#endif
