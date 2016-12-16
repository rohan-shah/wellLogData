#ifndef FEARNHEAD_RESAMPLING_HEADER_GUARD
#define FEARNHEAD_RESAMPLING_HEADER_GUARD
#include <functional>
#include <boost/random/mersenne_twister.hpp>
namespace wellLogData
{
	double fearnheadGetKappa(std::size_t size, std::function<double(std::size_t)> getWeight, boost::mt19937& randomSource, int N, int& A, double& B);
}
#endif
