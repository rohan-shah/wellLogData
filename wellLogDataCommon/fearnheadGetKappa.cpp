#include "fearnheadGetKappa.h"
#include <boost/random/uniform_int_distribution.hpp>
namespace wellLogData
{
	void fearnheadGetKappa(std::size_t size, std::function<double(std::size_t)> getWeight, boost::mt19937& randomSource, int N, int& A, double& B)
	{
		int lowerBound = 0, upperBound = (int)size-1;
		double lowerBoundSum = 0;
		const double tolerance = 1e-8;
		while(true)
		{
			//Split the existing partition between currentIndex and currentIndex2
			int currentIndex = boost::random::uniform_int_distribution<>(lowerBound, upperBound-1)(randomSource);
			while(currentIndex <= upperBound - 1 && getWeight(currentIndex) == getWeight(currentIndex+1)) currentIndex++;
			if(currentIndex == (int)size - 1)
			{
				A = 0;
				B = lowerBoundSum;
				for(int i = lowerBound; i <= upperBound; i++) B += getWeight(i);
			}
			else
			{
				double partitionWeight1 = getWeight(currentIndex), partitionWeight2 = getWeight(currentIndex+1);
				double B1 = lowerBoundSum;
				int A1 = (int)size - lowerBound;
				int i = lowerBound;
				for(; getWeight(i) < partitionWeight1 && i < (int)size; i++)
				{
					B1 += getWeight(i);
					A1--;
				}
				double B2 = B1;
				int A2 = A1;
				for(; getWeight(i) < partitionWeight2 && i < (int)size; i++)
				{
					B2 += getWeight(i);
					A2--;
				}
				//We've found the minimum
				if(B1 / partitionWeight1 + A1 > N + tolerance && B2 / partitionWeight2 + A2 <= N + tolerance)
				{
					B = B2;
					A = A2;
					return;
				}
				//Take the smaller half
				else if(B1 / partitionWeight1 + A1 <= N + tolerance)
				{
					upperBound = currentIndex;
				}
				//Take the bigger half
				else
				{
					lowerBound = currentIndex + 1;
					lowerBoundSum = B2;
				}
			}
			if(getWeight(lowerBound) == getWeight(upperBound))
			{
				double partitionWeight = getWeight(lowerBound);
				B = lowerBoundSum;
				A = (int)size - lowerBound;
				int i = lowerBound;
				for(; getWeight(i) < partitionWeight; i++)
				{
					B += getWeight(i);
					A--;
				}
				if(B / partitionWeight + A <= N + tolerance)
				{
					return;
				}
				else
				{
					for(; i <= upperBound; i++)
					{
						B += getWeight(i);
						A--;
					}
					return;
				}
			}
		}
	}
}
