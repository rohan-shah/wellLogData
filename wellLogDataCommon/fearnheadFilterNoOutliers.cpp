#include "fearnheadFilterNoOutliers.h"
#include "fearnheadGetKappa.h"
#include <boost/random/bernoulli_distribution.hpp>
#include "fearnheadSampling.h"
namespace wellLogData
{
	struct fearnheadFilterNoOutliersParticle
	{
	public:
		//mean and variance of X_t
		double mean, variance;
		std::vector<double> isChange;
		std::size_t timeLastChange;
		double weight;
		fearnheadFilterNoOutliersParticle()
			: mean(std::numeric_limits<double>::quiet_NaN()), variance(std::numeric_limits<double>::quiet_NaN()), timeLastChange(0), weight(std::numeric_limits<double>::quiet_NaN())
		{}
		fearnheadFilterNoOutliersParticle& operator=(fearnheadFilterNoOutliersParticle&& other)
		{
			mean = other.mean;
			variance = other.variance;
			isChange = std::move(other.isChange);
			weight = other.weight;
			return *this;
		}
		fearnheadFilterNoOutliersParticle(fearnheadFilterNoOutliersParticle&& other)
			: mean(other.mean), variance(other.variance), isChange(std::move(other.isChange)), weight(other.weight)
		{}
		fearnheadFilterNoOutliersParticle(const fearnheadFilterNoOutliersParticle& other)
			: mean(other.mean), variance(other.variance), isChange(other.isChange), weight(other.weight)
		{}
	};
	bool particleWeightSorter(const fearnheadFilterNoOutliersParticle& first, const fearnheadFilterNoOutliersParticle& second)
	{
		return first.weight < second.weight;
	}
	void fearnheadFilterNoOutliers(const context& contextObj, fearnheadFilterNoOutliersArgs& args)
	{
		if(args.nParticles == 0)
		{
			throw std::runtime_error("Input nParticles cannot be zero");
		}

		const std::vector<double>& data = contextObj.getData();

		args.changeProbabilities.clear();
		boost::random::bernoulli_distribution<> isChangeDist(contextObj.getChangeProbability());

		samplingDouble::fearnheadSamplingArgs samplingArgs;
		//Initial simulation
		std::vector<fearnheadFilterNoOutliersParticle> particles, childParticles;
		{
			fearnheadFilterNoOutliersParticle particle;
			particle.isChange.resize(1, 1.0);
			particle.timeLastChange = 0;
			particle.mean = contextObj.getMu();
			particle.variance = contextObj.getSigmaSquared();
			double tmp = (data[0] - particle.mean);
			particle.weight = (1/contextObj.getSigma()) * std::exp(-0.5 * tmp * tmp / particle.variance) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
			particles.emplace_back(std::move(particle));
		}
		std::vector<double> changeSums;
		for(std::size_t time = 1; time < data.size(); time++)
		{
			//For each particle, make four successors
			childParticles.clear();
			changeSums.resize(time);
			std::fill(changeSums.begin(), changeSums.end(), 0);
			double sumWeights = 0;
			for(std::size_t i = 0; i < particles.size(); i++)
			{
				fearnheadFilterNoOutliersParticle& currentParticle = particles[i];
				for(std::size_t time2 = 0; time2 < time; time2++)
				{
					changeSums[time2] += currentParticle.isChange[time2] * currentParticle.weight;
				}
				sumWeights += currentParticle.weight;
			
				fearnheadFilterNoOutliersParticle childParticle;
				childParticle.weight = currentParticle.weight;

				childParticle.isChange = currentParticle.isChange;
				childParticle.isChange.push_back(0.0);

				childParticle.timeLastChange = currentParticle.timeLastChange;
				childParticle.weight *= (1 - contextObj.getChangeProbability());

				childParticle.mean = ((currentParticle.mean / currentParticle.variance) + (data[time-1] / contextObj.getTau1Squared())) / ((1 / currentParticle.variance) + (1 / contextObj.getTau1Squared()));
				childParticle.variance = 1 / ((1 / currentParticle.variance) + (1 / contextObj.getTau1Squared()));
				double tmp = data[time] - childParticle.mean;
				double obsSd = sqrt(childParticle.variance + contextObj.getTau1Squared());
				childParticle.weight *= (1/sqrt(obsSd)) * std::exp(-0.5 * tmp * tmp / (childParticle.variance + contextObj.getTau1Squared())) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
				childParticles.emplace_back(std::move(childParticle));
			}
			for(std::size_t time2 = 0; time2 < time; time2++)
			{
				changeSums[time2] /= sumWeights;
			}
			//The *single* particle in the next generation corresponding to a change
			{
				double tmp = (data[time] - contextObj.getMu());
				double obsSd = sqrt(contextObj.getSigmaSquared() + contextObj.getTau1Squared());
				double changeWeight = (1/obsSd) * std::exp(-0.5 * tmp * tmp / (contextObj.getSigmaSquared() + contextObj.getTau1Squared())) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
				fearnheadFilterNoOutliersParticle changeParticle;
				changeParticle.weight = sumWeights * contextObj.getChangeProbability() * changeWeight;

				changeParticle.isChange = changeSums;
				changeParticle.isChange.push_back(1.0);

				changeParticle.timeLastChange = time;
				changeParticle.mean = contextObj.getMu();
				changeParticle.variance = contextObj.getSigmaSquared();
				childParticles.emplace_back(std::move(changeParticle));
			}
			if(time == data.size()-1)
			{
				particles.swap(childParticles);
				break;
			}
			if(childParticles.size() <= args.nParticles)
			{
				particles.swap(childParticles);
			}
			else
			{
				//Resampling according to the method outlined in Fearnhead. 
				samplingArgs.weights.clear();
				for(int i = 0; i < (int)childParticles.size(); i++) samplingArgs.weights.push_back(childParticles[i].weight);

				samplingArgs.n = args.nParticles;
				samplingDouble::fearnheadSampling(samplingArgs, args.randomSource);

				particles.clear();
				for(int i = 0; i < (int)samplingArgs.indices.size(); i++)
				{
					particles.emplace_back(std::move(childParticles[samplingArgs.indices[i]]));
					//multiple all particle weights by c and divide by the inclusion probabilities
					if(samplingArgs.deterministicInclusion[samplingArgs.indices[i]])
					{
						particles.back().weight *= samplingArgs.c;
					}
					else
					{
						particles.back().weight = 1;
					}
				}
			}
		}
		args.changeProbabilities.clear();
		double sumWeights = 0;
		for(std::size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++) sumWeights += particles[particleIndex].weight;
		for(std::size_t time = 0; time < data.size(); time++)
		{
			double sumChange = 0;
			for(std::size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
			{
				sumChange += particles[particleIndex].isChange[time] * particles[particleIndex].weight;
			}
			args.changeProbabilities.push_back(sumChange / sumWeights);
		}
	}
}

