#include "basicFilter.h"
#include <boost/random/bernoulli_distribution.hpp>
#include "aliasMethod.h"
namespace wellLogData
{
	struct basicFilterParticle
	{
		//mean and variance of X_t
		double mean, variance;
		//hidden state
		bool isOutlier, isChange;
		std::size_t timeLastChange;
		double weight;
	};
	void basicFilter(const context& contextObj, basicFilterArgs& args)
	{
		if(args.nParticles == 0)
		{
			throw std::runtime_error("Input nParticles cannot be zero");
		}

		const std::vector<double>& data = contextObj.getData();

		args.changeProbabilities.clear();
		args.outlierProbabilities.clear();
		boost::random::bernoulli_distribution<> isOutlierDist(contextObj.getOutlierProbability()), isChangeDist(contextObj.getChangeProbability());

		//Initial simulation
		std::vector<basicFilterParticle> particles, childParticles;
		for(std::size_t i = 0; i < args.nParticles; i++)
		{
			basicFilterParticle particle;
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
				basicFilterParticle& currentParticle = particles[i];
				for(int j = 0; j < 4; j++)
				{
					basicFilterParticle childParticle;
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
					if(childParticle.isChange)
					{
						childParticle.timeLastChange = time;
						childParticle.weight *= contextObj.getChangeProbability();
						childParticle.mean = contextObj.getMu();
						childParticle.variance = contextObj.getSigmaSquared();
						double tmp = (data[time] - childParticle.mean);
						childParticle.weight *= (1/contextObj.getSigma()) * std::exp(-0.5 * tmp * tmp / childParticle.variance) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
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
							childParticle.mean = (currentParticle.mean / currentParticle.variance) + (data[time-1] / contextObj.getTau1Squared()) / ((1 / currentParticle.variance) + (1 / contextObj.getTau1Squared()));
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
							childParticle.weight *= (1/sqrt(childParticle.variance)) * std::exp(-0.5 * tmp * tmp / childParticle.variance) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
						}
					}
					childParticles.push_back(childParticle);
				}
			}
			resamplingProbabilities.clear();
			double sum = 0, outlierProbabilitiesSum = 0, changeProbabilitiesSum = 0;
			for(std::vector<basicFilterParticle>::iterator i = childParticles.begin(); i != childParticles.end(); i++)
			{
				resamplingProbabilities.push_back(i->weight);
				sum += i->weight;
				if(i->isChange) changeProbabilitiesSum += i->weight;
				if(i->isOutlier) outlierProbabilitiesSum += i->weight;
			}
			args.changeProbabilities.push_back(changeProbabilitiesSum / sum);
			args.outlierProbabilities.push_back(outlierProbabilitiesSum / sum);

			aliasMethod::aliasMethod resamplingObject(resamplingProbabilities, sum, aliasTmp1, aliasTmp2, aliasTmp3);
			particles.clear();
			for(std::size_t i = 0; i < args.nParticles; i++)
			{
				std::size_t index = resamplingObject(args.randomSource);
				particles.push_back(childParticles[index]);
				particles.back().weight = sum / args.nParticles;
			}
		}
	}
}

