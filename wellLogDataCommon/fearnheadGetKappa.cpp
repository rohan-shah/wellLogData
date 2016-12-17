#include "fearnheadGetKappa.h"
#include <boost/random/uniform_int_distribution.hpp>
namespace wellLogData
{
	double fearnheadGetKappa(std::size_t size, std::function<double(std::size_t)> getWeight, boost::mt19937& randomSource, int N, int& A, double& B)
	{
		int lowerBound = 0, upperBound = (int)size-1;
		double lowerBoundSum = 0;
		const double tolerance = 1e-8;
		while(getWeight(lowerBound) != getWeight(upperBound))
		{
			//Split the existing partition between currentIndex and currentIndex2
			int currentIndex1 = boost::random::uniform_int_distribution<>(lowerBound, upperBound-1)(randomSource);
			int currentIndex2 = currentIndex1 + 1;
			double partitionWeight1 = getWeight(currentIndex1), partitionWeight2 = getWeight(currentIndex2);
			if(partitionWeight1 == partitionWeight2)
			{
				B = lowerBoundSum;
				A = (int)(size - lowerBound);
				for(int i = lowerBound; getWeight(i) <= partitionWeight1 && i < (int)size; i++)
				{
					B += getWeight(i);
					A--;
				}
				if(B / partitionWeight1 + A <= N + tolerance)
				{
					upperBound = currentIndex1;
				}
				else
				{
					lowerBoundSum = B;
					lowerBound = (int)size - A;
					if(A == 0) goto skipCheck;
				}
			}
			else
			{
				double B1 = lowerBoundSum;
				int A1 = (int)size - lowerBound;
				for(int i = lowerBound; getWeight(i) <= partitionWeight1 && i < (int)size; i++)
				{
					B1 += getWeight(i);
					A1--;
				}
				double B2 = B1;
				int A2 = A1;
				for(int i = currentIndex2; getWeight(i) <= partitionWeight2 && i < (int)size; i++)
				{
					B2 += getWeight(i);
					A2--;
				}
				//We've found the minimum
				if(B1 / partitionWeight1 + A1 > N + tolerance && B2 / partitionWeight2 + A2 <= N + tolerance)
				{
					upperBound = lowerBound = currentIndex2;
					B = B2;
					A = A2;
					goto skipCheck;
				}
				//Take the smaller half
				else if(B1 / partitionWeight1 + A1 <= N + tolerance)
				{
					upperBound = currentIndex1;
					A = A1;
					B = B1;
				}
				//Take the bigger half
				else
				{
					lowerBound = currentIndex2;
					lowerBoundSum = B1;
					B = B2;
					A = A2;
				}
			}
		}
		{
			double partitionWeight = getWeight(lowerBound);
			A = (int)size - lowerBound;
			B = lowerBoundSum;
			for(int i = lowerBound; getWeight(i) <= partitionWeight && i < (int)size; i++)
			{
				B += getWeight(i);
				A--;
			}
		}
skipCheck:
		return getWeight(lowerBound);
	}
}
