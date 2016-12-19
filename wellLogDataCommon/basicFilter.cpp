#include "basicFilter.h"
#include <boost/random/bernoulli_distribution.hpp>
#include "aliasMethod.h"
namespace wellLogData
{
	struct basicFilterParticle
	{
	public:
		//mean and variance of X_t
		double mean, variance;
		std::size_t timeLastChange;
		double weight;
		//History of the hidden states
		std::vector<bool> isChange;
		std::vector<bool> isOutlier;
		basicFilterParticle()
			:mean(std::numeric_limits<double>::quiet_NaN()), variance(std::numeric_limits<double>::quiet_NaN()), timeLastChange(0), weight(std::numeric_limits<double>::quiet_NaN()) 
		{}
		basicFilterParticle(basicFilterParticle&& other)
			:mean(other.mean), variance(other.variance), timeLastChange(other.timeLastChange), weight(other.weight), isChange(std::move(other.isChange)), isOutlier(std::move(other.isOutlier))
		{}
		basicFilterParticle(const basicFilterParticle& other)
			:mean(other.mean), variance(other.variance), timeLastChange(other.timeLastChange), weight(other.weight), isChange(other.isChange), isOutlier(other.isOutlier)
		{}
		basicFilterParticle operator=(basicFilterParticle&& other)
		{
			mean = other.mean;
			variance = other.variance;
			timeLastChange = other.timeLastChange;
			weight = other.weight;
			isChange = std::move(other.isChange);
			isOutlier = std::move(other.isOutlier);
			return *this;
		}
	};
	void basicFilter(const context& contextObj, basicFilterArgs& args)
	{
		if(args.nParticles == 0)
		{
			throw std::runtime_error("Input nParticles cannot be zero");
		}

		const std::vector<double>& data = contextObj.getData();
		boost::random::bernoulli_distribution<> isOutlierDist(contextObj.getOutlierProbability()), isChangeDist(contextObj.getChangeProbability());

		//Initial simulation
		std::vector<basicFilterParticle> particles, childParticles;
		for(std::size_t i = 0; i < args.nParticles; i++)
		{
			basicFilterParticle particle;
			particle.isChange.resize(1, true);
			particle.isOutlier.resize(1, isOutlierDist(args.randomSource));
			particle.timeLastChange = 0;
			particle.mean = contextObj.getMu();
			particle.variance = contextObj.getSigmaSquared();
			if(particle.isOutlier[0])
			{
				double tmp = (data[0] - contextObj.getNu());
				particle.weight = contextObj.getOutlierProbability() * (1/contextObj.getTau2()) * std::exp(-0.5 * tmp * tmp / contextObj.getTau2Squared()) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
			}
			else
			{
				double tmp = (data[0] - particle.mean);
				particle.weight = (1 - contextObj.getOutlierProbability()) * (1/contextObj.getSigma()) * std::exp(-0.5 * tmp * tmp / particle.variance) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
			}
			particles.emplace_back(std::move(particle));
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
					childParticle.isChange = currentParticle.isChange;
					childParticle.isOutlier = currentParticle.isOutlier;
					childParticle.isChange.push_back(j % 2);
					childParticle.isOutlier.push_back(j / 2);
					if(currentParticle.isOutlier.back())
					{
						if(childParticle.isOutlier.back())
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
						if(childParticle.isOutlier.back())
						{
							childParticle.weight *= contextObj.getOutlierProbability();
						}
						else
						{
							childParticle.weight *= (1 - contextObj.getOutlierProbability());
						}
					}
					if(childParticle.isChange.back() && childParticle.isOutlier.back())
					{
						childParticle.timeLastChange = time;
						childParticle.weight *= contextObj.getChangeProbability();
						childParticle.mean = contextObj.getMu();
						childParticle.variance = contextObj.getSigmaSquared();
						double tmp = (data[time] - contextObj.getNu());
						childParticle.weight *= (1/contextObj.getTau2()) * std::exp(-0.5 * tmp * tmp / contextObj.getTau2Squared()) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
					}
					else if(childParticle.isChange.back())
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
						if(currentParticle.isOutlier.back())
						{
							childParticle.mean = currentParticle.mean;
							childParticle.variance = currentParticle.variance;
						}
						else
						{
							childParticle.mean = ((currentParticle.mean / currentParticle.variance) + (data[time-1] / contextObj.getTau1Squared())) / ((1 / currentParticle.variance) + (1 / contextObj.getTau1Squared()));
							childParticle.variance = 1 / ((1 / currentParticle.variance) + (1 / contextObj.getTau1Squared()));
						}
						if(childParticle.isOutlier.back())
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
					childParticles.emplace_back(std::move(childParticle));
				}
			}
			if(time == data.size()-1) 
			{
				particles.swap(childParticles);
				break;
			}
			resamplingProbabilities.clear();
			double sum = 0;
			for(std::vector<basicFilterParticle>::iterator i = childParticles.begin(); i != childParticles.end(); i++)
			{
				resamplingProbabilities.push_back(i->weight);
				sum += i->weight;
			}

			aliasMethod::aliasMethod resamplingObject(resamplingProbabilities, sum, aliasTmp1, aliasTmp2, aliasTmp3);
			particles.clear();
			for(std::size_t i = 0; i < args.nParticles; i++)
			{
				std::size_t index = resamplingObject(args.randomSource);
				particles.push_back(childParticles[index]);
				particles.back().weight = 1;
				//particles.back().weight = sum / args.nParticles;
			}
		}
		args.changeProbabilities.clear();
		args.outlierProbabilities.clear();
		double sumWeights = 0;
		for(std::size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++) sumWeights += particles[particleIndex].weight;
		for(std::size_t time = 0; time < data.size(); time++)
		{
			double sumOutlier = 0, sumChange = 0;
			for(std::size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
			{
				if(particles[particleIndex].isChange[time]) sumChange += particles[particleIndex].weight;
				if(particles[particleIndex].isOutlier[time]) sumOutlier += particles[particleIndex].weight;
			}
			args.changeProbabilities.push_back(sumChange / sumWeights);
			args.outlierProbabilities.push_back(sumOutlier / sumWeights);
		}
	}
}

