#include "fearnheadFilter.h"
#include "fearnheadGetKappa.h"
#include <boost/random/bernoulli_distribution.hpp>
#include "systematicSampling.h"
namespace wellLogData
{
	struct fearnheadFilterParticle
	{
		//mean and variance of X_t
		double mean, variance;
		//hidden state
		bool isOutlier, isChange;
		std::size_t timeLastChange;
		double weight;
	};
	bool particleWeightSorter(const fearnheadFilterParticle& first, const fearnheadFilterParticle& second)
	{
		return first.weight < second.weight;
	}
	void fearnheadFilter(const context& contextObj, fearnheadFilterArgs& args)
	{
		if(args.nParticles == 0)
		{
			throw std::runtime_error("Input nParticles cannot be zero");
		}

		const std::vector<double>& data = contextObj.getData();

		args.changeProbabilities.clear();
		args.outlierProbabilities.clear();
		boost::random::bernoulli_distribution<> isOutlierDist(contextObj.getOutlierProbability()), isChangeDist(contextObj.getChangeProbability());

		std::vector<double> systematicWeights;
		std::vector<int> systematicIndices;

		//Initial simulation
		std::vector<fearnheadFilterParticle> particles, childParticles;
		for(std::size_t i = 0; i < args.nParticles; i++)
		{
			fearnheadFilterParticle particle;
			particle.isChange = true;
			particle.isOutlier = isOutlierDist(args.randomSource);
			particle.timeLastChange = 0;
			particle.mean = contextObj.getMu();
			particle.variance = contextObj.getSigmaSquared();
			if(particle.isOutlier)
			{
				double tmp = (data[0] - contextObj.getNu());
				particle.weight = contextObj.getOutlierProbability() * (1/contextObj.getTau2()) * std::exp(-0.5 * tmp * tmp / contextObj.getTau2Squared()) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
			}
			else
			{
				double tmp = (data[0] - particle.mean);
				particle.weight = (1 - contextObj.getOutlierProbability()) * (1/contextObj.getSigma()) * std::exp(-0.5 * tmp * tmp / particle.variance) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
			}
			particles.push_back(particle);
		}
		std::vector<double> resamplingProbabilities;
		std::vector<std::ptrdiff_t> aliasTmp1, aliasTmp2;
		std::vector<std::pair<double, std::ptrdiff_t> > aliasTmp3;
		for(std::size_t time = 1; time < data.size(); time++)
		{
			//For each particle, make four successors
			childParticles.clear();
			for(std::size_t i = 0; i < args.nParticles; i++)
			{
				fearnheadFilterParticle& currentParticle = particles[i];
				for(int j = 0; j < 4; j++)
				{
					fearnheadFilterParticle childParticle;
					childParticle.weight = currentParticle.weight;
					childParticle.isChange = j % 2;
					childParticle.isOutlier = j / 2;
					if(currentParticle.isOutlier)
					{
						if(childParticle.isOutlier)
						{
							childParticle.weight *= contextObj.getOutlierClusterProbability();
						}
						else
						{
							childParticle.weight *= (1 - contextObj.getOutlierClusterProbability());
						}
					}
					else
					{
						if(childParticle.isOutlier)
						{
							childParticle.weight *= contextObj.getOutlierProbability();
						}
						else
						{
							childParticle.weight *= (1 - contextObj.getOutlierProbability());
						}
					}
					if(childParticle.isChange && childParticle.isOutlier)
					{
						childParticle.timeLastChange = time;
						childParticle.weight *= contextObj.getChangeProbability();
						childParticle.mean = contextObj.getMu();
						childParticle.variance = contextObj.getSigmaSquared();
						double tmp = (data[time] - contextObj.getNu());
						childParticle.weight *= (1/contextObj.getTau2()) * std::exp(-0.5 * tmp * tmp / contextObj.getTau2Squared()) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
					}
					else if(childParticle.isChange)
					{
						childParticle.timeLastChange = time;
						childParticle.weight *= contextObj.getChangeProbability();
						childParticle.mean = contextObj.getMu();
						childParticle.variance = contextObj.getSigmaSquared();
						double tmp = (data[time] - childParticle.mean);
						double obsSd = sqrt(contextObj.getSigmaSquared() + contextObj.getTau1Squared());
						childParticle.weight *= (1/obsSd) * std::exp(-0.5 * tmp * tmp / (childParticle.variance + contextObj.getTau1Squared())) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
					}
					else 
					{
						childParticle.timeLastChange = currentParticle.timeLastChange;
						childParticle.weight *= (1 - contextObj.getChangeProbability());
						if(currentParticle.isOutlier)
						{
							childParticle.mean = currentParticle.mean;
							childParticle.variance = currentParticle.variance;
						}
						else
						{
							childParticle.mean = ((currentParticle.mean / currentParticle.variance) + (data[time-1] / contextObj.getTau1Squared())) / ((1 / currentParticle.variance) + (1 / contextObj.getTau1Squared()));
							childParticle.variance = 1 / ((1 / currentParticle.variance) + (1 / contextObj.getTau1Squared()));
						}
						if(childParticle.isOutlier)
						{
							double tmp = data[time] - contextObj.getNu();
							childParticle.weight *= (1/contextObj.getTau2()) * std::exp(-0.5 * tmp * tmp / contextObj.getTau2Squared()) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
						}
						else
						{
							double tmp = data[time] - childParticle.mean;
							double obsSd = sqrt(childParticle.variance + contextObj.getTau1Squared());
							childParticle.weight *= (1/sqrt(obsSd)) * std::exp(-0.5 * tmp * tmp / (childParticle.variance + contextObj.getTau1Squared())) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
						}
					}
					childParticles.push_back(childParticle);
				}
			}
			resamplingProbabilities.clear();
			double sum = 0, outlierProbabilitiesSum = 0, changeProbabilitiesSum = 0;
			for(std::vector<fearnheadFilterParticle>::iterator i = childParticles.begin(); i != childParticles.end(); i++)
			{
				resamplingProbabilities.push_back(i->weight);
				sum += i->weight;
				if(i->isChange) changeProbabilitiesSum += i->weight;
				if(i->isOutlier) outlierProbabilitiesSum += i->weight;
			}
			args.changeProbabilities.push_back(changeProbabilitiesSum / sum);
			args.outlierProbabilities.push_back(outlierProbabilitiesSum / sum);

			//Resampling according to the method outlined in Fearnhead. 
			//First identify the value of c

			//Sort by increasing weighti
			std::sort(childParticles.begin(), childParticles.end(), particleWeightSorter);
			int A;
			double B;
			fearnheadGetKappa(childParticles.size(), [&childParticles](std::size_t i){ return childParticles[i].weight;}, args.randomSource, args.nParticles, A, B);
			double c = ((int)args.nParticles - A) / B;
#ifndef NDEBUG
			double checkSum = 0, checkSum2 = 0;
			for(std::vector<fearnheadFilterParticle>::iterator i = childParticles.begin(); i != childParticles.end(); i++)
			{
				checkSum += std::min(i->weight * c, 1.0);
				checkSum2 += i->weight;
			}
			if(fabs(checkSum - args.nParticles) > 1e-6)
			{
				throw std::runtime_error("Internal error");
			}
#endif
			double cInverse = 1/c;
			//Work out where "set 1" starts
			int startOfTakeAllStrata = 0;
			while(startOfTakeAllStrata < (int)childParticles.size() && childParticles[startOfTakeAllStrata].weight < cInverse) startOfTakeAllStrata++;
			int takeAllStrataSize = (int)childParticles.size() - startOfTakeAllStrata;

			//The rest of the particles form "set 2"
			systematicWeights.clear();
			double setTwoSum = 0;
			for(int i = 0; i < startOfTakeAllStrata; i++)
			{
				systematicWeights.push_back(childParticles[i].weight);
				setTwoSum += childParticles[i].weight;
			}
			systematicIndices.clear();
			sampling::systematicSamplingDouble(systematicWeights, setTwoSum / (args.nParticles - takeAllStrataSize), systematicIndices, args.randomSource);
			particles.clear();
			for(int i = 0; i < (int)systematicIndices.size(); i++)
			{
				particles.push_back(childParticles[systematicIndices[i]]);
				particles.back().weight = cInverse;
			}
			for(int i = startOfTakeAllStrata; i < (int)systematicIndices.size(); i++)
			{
				particles.push_back(childParticles[i]);
			}
		}
	}
}

