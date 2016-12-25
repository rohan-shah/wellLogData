#include "fearnheadFilter.h"
#include "fearnheadGetKappa.h"
#include <boost/random/bernoulli_distribution.hpp>
#include "systematicSampling.h"
namespace wellLogData
{
	struct fearnheadFilterParticle
	{
	public:
		//mean and variance of X_t
		double mean, variance;
		std::vector<double> isOutlier, isChange;
		std::size_t timeLastChange;
		double weight;
		fearnheadFilterParticle()
			: mean(std::numeric_limits<double>::quiet_NaN()), variance(std::numeric_limits<double>::quiet_NaN()), timeLastChange(0), weight(std::numeric_limits<double>::quiet_NaN())
		{}
		fearnheadFilterParticle& operator=(fearnheadFilterParticle&& other)
		{
			mean = other.mean;
			variance = other.variance;
			isOutlier = std::move(other.isOutlier);
			isChange = std::move(other.isChange);
			weight = other.weight;
			return *this;
		}
		fearnheadFilterParticle(fearnheadFilterParticle&& other)
			: mean(other.mean), variance(other.variance), isOutlier(std::move(other.isOutlier)), isChange(std::move(other.isChange)), weight(other.weight)
		{}
		fearnheadFilterParticle(const fearnheadFilterParticle& other)
			: mean(other.mean), variance(other.variance), isOutlier(other.isOutlier), isChange(other.isChange), weight(other.weight)
		{}
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
			particle.isChange.resize(1, 1.0);
			particle.isOutlier.resize(1, (double)isOutlierDist(args.randomSource));
			particle.timeLastChange = 0;
			particle.mean = contextObj.getMu();
			particle.variance = contextObj.getSigmaSquared();
			if(particle.isOutlier[0] != 0)
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
		std::vector<double> changeSums, outlierSums;
		for(std::size_t time = 1; time < data.size(); time++)
		{
			//For each particle, make four successors
			childParticles.clear();
			changeSums.resize(time);
			outlierSums.resize(time);
			std::fill(changeSums.begin(), changeSums.end(), 0);
			std::fill(outlierSums.begin(), outlierSums.end(), 0);
			double sumWeights = 0, outlierSumWeights = 0, notOutlierSumWeights = 0;
			for(std::size_t i = 0; i < particles.size(); i++)
			{
				fearnheadFilterParticle& currentParticle = particles[i];
#ifndef NDBEUG
				if(currentParticle.isOutlier[time-1] != 0 && currentParticle.isOutlier[time-1] != 1) throw std::runtime_error("Internal error");
#endif
				for(std::size_t time2 = 0; time2 < time; time2++)
				{
					outlierSums[time2] += currentParticle.isOutlier[time2] * currentParticle.weight;
					changeSums[time2] += currentParticle.isChange[time2] * currentParticle.weight;
				}
				sumWeights += currentParticle.weight;
				if(currentParticle.isOutlier[time-1] > 0) outlierSumWeights += currentParticle.weight;
				else notOutlierSumWeights += currentParticle.weight;
				for(int j = 0; j < 2; j++)
				{
					fearnheadFilterParticle childParticle;
					childParticle.weight = currentParticle.weight;

					childParticle.isChange = currentParticle.isChange;
					childParticle.isChange.push_back(0.0);

					childParticle.isOutlier = currentParticle.isOutlier;
					childParticle.isOutlier.push_back((double)(j % 2));

					if(currentParticle.isOutlier[time-1] > 0)
					{
						if(j % 2)
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
						if(j % 2)
						{
							childParticle.weight *= contextObj.getOutlierProbability();
						}
						else
						{
							childParticle.weight *= (1 - contextObj.getOutlierProbability());
						}
					}
					childParticle.timeLastChange = currentParticle.timeLastChange;
					childParticle.weight *= (1 - contextObj.getChangeProbability());
					if(j % 2)
					{
						childParticle.mean = currentParticle.mean;
						childParticle.variance = currentParticle.variance;
					}
					else
					{
						childParticle.mean = ((currentParticle.mean / currentParticle.variance) + (data[time-1] / contextObj.getTau1Squared())) / ((1 / currentParticle.variance) + (1 / contextObj.getTau1Squared()));
						childParticle.variance = 1 / ((1 / currentParticle.variance) + (1 / contextObj.getTau1Squared()));
					}
					if(childParticle.isOutlier[time])
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
					childParticles.emplace_back(std::move(childParticle));
				}
			}
			for(std::size_t time2 = 0; time2 < time; time2++)
			{
				outlierSums[time2] /= sumWeights;
				changeSums[time2] /= sumWeights;
			}
			//The *single* particle in the next generation corresponding to both an outlier *and* a change
			{
				double tmp = (data[time] - contextObj.getNu());
				double changeOutlierWeight = (1/contextObj.getTau2()) * std::exp(-0.5 * tmp * tmp / contextObj.getTau2Squared()) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
				fearnheadFilterParticle changeOutlierParticle;
				changeOutlierParticle.weight = (outlierSumWeights * contextObj.getOutlierClusterProbability() + notOutlierSumWeights * contextObj.getOutlierProbability()) * contextObj.getChangeProbability() * changeOutlierWeight;
				changeOutlierParticle.isChange = changeSums;
				changeOutlierParticle.isChange.push_back(1.0);

				changeOutlierParticle.isOutlier = outlierSums;
				changeOutlierParticle.isOutlier.push_back(1.0);

				changeOutlierParticle.timeLastChange = time;
				changeOutlierParticle.mean = contextObj.getMu();
				changeOutlierParticle.variance = contextObj.getSigmaSquared();
				childParticles.emplace_back(std::move(changeOutlierParticle));
			}
			//The *single* particle in the next generation corresponding to a change (but not an outlier)
			{
				double tmp = (data[time] - contextObj.getMu());
				double obsSd = sqrt(contextObj.getSigmaSquared() + contextObj.getTau1Squared());
				double changeWeight = (1/obsSd) * std::exp(-0.5 * tmp * tmp / (contextObj.getSigmaSquared() + contextObj.getTau1Squared())) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
				fearnheadFilterParticle changeParticle;
				changeParticle.weight = (outlierSumWeights * (1 - contextObj.getOutlierClusterProbability()) + notOutlierSumWeights * (1 - contextObj.getOutlierProbability())) * contextObj.getChangeProbability() * changeWeight;

				changeParticle.isChange = changeSums;
				changeParticle.isChange.push_back(1.0);

				changeParticle.isOutlier = outlierSums;
				changeParticle.isOutlier.push_back(0.0);

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
				samplingDouble::systematicSamplingDouble(systematicWeights, setTwoSum / (args.nParticles - takeAllStrataSize), systematicIndices, args.randomSource);
				particles.clear();
				for(int i = 0; i < (int)systematicIndices.size(); i++)
				{
					particles.emplace_back(std::move(childParticles[systematicIndices[i]]));
					//The particle weight is reset to 1 here, because we want to avoid the weights going to zero. To do this we divide through by the smallest particle weight, which is cInverse. 
					particles.back().weight = 1;
				}
				for(int i = startOfTakeAllStrata; i < (int)childParticles.size(); i++)
				{
					particles.emplace_back(std::move(childParticles[i]));
					particles.back().weight *= c;
				}
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
				sumChange += particles[particleIndex].isChange[time] * particles[particleIndex].weight;
				sumOutlier += particles[particleIndex].isOutlier[time] * particles[particleIndex].weight;
			}
			args.changeProbabilities.push_back(sumChange / sumWeights);
			args.outlierProbabilities.push_back(sumOutlier / sumWeights);
		}
	}
}

