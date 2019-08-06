#include "fearnheadFilterOptimised.h"
#include "fearnheadGetKappa.h"
#include <boost/random/bernoulli_distribution.hpp>
#include "fearnheadSampling.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
namespace wellLogData
{
	typedef float numericType;
	struct fearnheadFilterOptimisedParticle
	{
	public:
		//mean and variance of X_t
		numericType mean, variance;
		std::size_t timeLastChange;
		numericType weight;
		fearnheadFilterOptimisedParticle()
			: mean(std::numeric_limits<numericType>::quiet_NaN()), variance(std::numeric_limits<numericType>::quiet_NaN()), timeLastChange(0), weight(std::numeric_limits<numericType>::quiet_NaN())
		{}
	};
	struct customMatrix : public std::vector < std::vector<numericType> >
	{
	public:
		customMatrix(std::size_t nRow, std::size_t nCol)
		{
			for (std::size_t i = 0; i < nRow; i++)
			{
				this->emplace_back(std::vector<numericType>(nCol, 0));
			}
		}
		numericType& operator()(std::size_t row, std::size_t col)
		{
			return this->operator[](row)[col];
		}
	};
	bool particleWeightSorter(const fearnheadFilterOptimisedParticle& first, const fearnheadFilterOptimisedParticle& second)
	{
		return first.weight < second.weight;
	}
	void fearnheadFilterOptimised(const context& contextObj, fearnheadFilterArgs& args)
	{
		if(args.nParticles == 0)
		{
			throw std::runtime_error("Input nParticles cannot be zero");
		}

		const std::vector<double>& originalData = contextObj.getData();
		std::vector<float> data(originalData.begin(), originalData.end());

		typedef customMatrix probabilityMatrix;
		probabilityMatrix changeProbabilities(2*args.nParticles + 2, data.size()), outlierProbabilities(2 * args.nParticles + 2, data.size());
		probabilityMatrix changeProbabilitiesPossibilities(2 * args.nParticles + 2, data.size()), outlierProbabilitiesPossibilities(2 * args.nParticles + 2, data.size());
		args.changeProbabilities.clear();
		args.outlierProbabilities.clear();
		boost::random::bernoulli_distribution<> isOutlierDist(contextObj.getOutlierProbability()), isChangeDist(contextObj.getChangeProbability());

		samplingDouble::fearnheadSamplingArgs samplingArgs;
		//Initial simulation
		std::vector<fearnheadFilterOptimisedParticle> particles, childParticles;
		for(std::size_t i = 0; i < 2; i++)
		{
			fearnheadFilterOptimisedParticle particle;
			changeProbabilities(i, 0) = 1.0;
			outlierProbabilities(i, 0) = (numericType)i;
			particle.timeLastChange = 0;
			particle.mean = contextObj.getMu();
			particle.variance = contextObj.getSigmaSquared();
			if(outlierProbabilities(i, 0) != 0)
			{
				numericType tmp = (data[0] - contextObj.getNu());
				particle.weight = contextObj.getOutlierProbability() * (1/contextObj.getTau2()) * std::exp(-0.5 * tmp * tmp / contextObj.getTau2Squared()) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
			}
			else
			{
				numericType tmp = (data[0] - particle.mean);
				particle.weight = (1 - contextObj.getOutlierProbability()) * (1/contextObj.getSigma()) * std::exp(-0.5 * tmp * tmp / particle.variance) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
			}
			particles.emplace_back(std::move(particle));
		}
		std::vector<numericType> changeSums, outlierSums;
		for (std::size_t time = 1; time < data.size(); time++)
		{
			//For each particle, make four successors
			childParticles.clear();
			changeSums.resize(time);
			outlierSums.resize(time);
			std::fill(changeSums.begin(), changeSums.end(), 0);
			std::fill(outlierSums.begin(), outlierSums.end(), 0);
			numericType sumWeights = 0, outlierSumWeights = 0, notOutlierSumWeights = 0;
			for (std::size_t i = 0; i < particles.size(); i++)
			{
				fearnheadFilterOptimisedParticle& currentParticle = particles[i];
				std::vector<numericType>& currentOutlierRow = outlierProbabilities[i];
				std::vector<numericType>& currentChangeRow = changeProbabilities[i];
#ifndef NDBEUG
				if (currentOutlierRow[time - 1] != 0 && currentOutlierRow[time - 1] != 1) throw std::runtime_error("Internal error");
#endif
				for (std::size_t time2 = 0; time2 < time; time2++)
				{
					outlierSums[time2] += currentOutlierRow[time2] * currentParticle.weight;
					changeSums[time2] += currentChangeRow[time2] * currentParticle.weight;
				}
				sumWeights += currentParticle.weight;
				if (currentOutlierRow[time - 1] > 0) outlierSumWeights += currentParticle.weight;
				else notOutlierSumWeights += currentParticle.weight;

				fearnheadFilterOptimisedParticle childParticleOutlier, childParticleNotOutlier;
				{
					childParticleOutlier.weight = currentParticle.weight;

					std::vector<numericType>& destChangeRow = changeProbabilitiesPossibilities[2 * i + 0];
					std::vector<numericType>& destOutlierRow = outlierProbabilitiesPossibilities[2 * i + 0];

					std::copy(currentChangeRow.begin(), currentChangeRow.begin() + time, destChangeRow.begin());
					std::copy(currentOutlierRow.begin(), currentOutlierRow.begin() + time, destOutlierRow.begin());
					destOutlierRow[time] = 1.0;

					if (currentOutlierRow[time - 1] > 0)
					{
						childParticleOutlier.weight *= contextObj.getOutlierClusterProbability();
					}
					else
					{
						childParticleOutlier.weight *= contextObj.getOutlierProbability();
					}
					childParticleOutlier.timeLastChange = currentParticle.timeLastChange;
					childParticleOutlier.weight *= (1 - contextObj.getChangeProbability());
					if (destOutlierRow[time - 1] > 0.0)
					{
						childParticleOutlier.mean = currentParticle.mean;
						childParticleOutlier.variance = currentParticle.variance;
					}
					else
					{
						childParticleOutlier.mean = ((currentParticle.mean / currentParticle.variance) + (data[time - 1] / contextObj.getTau1Squared())) / ((1 / currentParticle.variance) + (1 / contextObj.getTau1Squared()));
						childParticleOutlier.variance = 1 / ((1 / currentParticle.variance) + (1 / contextObj.getTau1Squared()));
					}
					numericType tmp = data[time] - contextObj.getNu();
					childParticleOutlier.weight *= (1 / contextObj.getTau2()) * std::exp(-0.5 * tmp * tmp / contextObj.getTau2Squared()) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
				}
				{
					childParticleNotOutlier.weight = currentParticle.weight;

					std::vector<numericType>& destChangeRow = changeProbabilitiesPossibilities[2 * i + 1];
					std::vector<numericType>& destOutlierRow = outlierProbabilitiesPossibilities[2 * i + 1];

					destChangeRow.swap(currentChangeRow);
					destOutlierRow.swap(currentOutlierRow);

					if (currentOutlierRow[time - 1] > 0)
					{
						childParticleNotOutlier.weight *= (1 - contextObj.getOutlierClusterProbability());
					}
					else
					{
						childParticleNotOutlier.weight *= (1 - contextObj.getOutlierProbability());
					}
					childParticleNotOutlier.timeLastChange = currentParticle.timeLastChange;
					childParticleNotOutlier.weight *= (1 - contextObj.getChangeProbability());
					if (destOutlierRow[time - 1] > 0.0)
					{
						childParticleNotOutlier.mean = currentParticle.mean;
						childParticleNotOutlier.variance = currentParticle.variance;
					}
					else
					{
						childParticleNotOutlier.mean = ((currentParticle.mean / currentParticle.variance) + (data[time - 1] / contextObj.getTau1Squared())) / ((1 / currentParticle.variance) + (1 / contextObj.getTau1Squared()));
						childParticleNotOutlier.variance = 1 / ((1 / currentParticle.variance) + (1 / contextObj.getTau1Squared()));
					}
					numericType tmp = data[time] - childParticleNotOutlier.mean;
					numericType obsSd = sqrt(childParticleNotOutlier.variance + contextObj.getTau1Squared());
					childParticleNotOutlier.weight *= (1 / sqrt(obsSd)) * std::exp(-0.5 * tmp * tmp / (childParticleNotOutlier.variance + contextObj.getTau1Squared())) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
				}
				childParticles.emplace_back(std::move(childParticleOutlier));
				childParticles.emplace_back(std::move(childParticleNotOutlier));
			}
			if(sumWeights == 0)
			{
				throw std::runtime_error("Internal error");
			}
			for(std::size_t time2 = 0; time2 < time; time2++)
			{
				outlierSums[time2] /= sumWeights;
				changeSums[time2] /= sumWeights;
			}
			//The *single* particle in the next generation corresponding to both an outlier *and* a change
			{
				numericType tmp = (data[time] - contextObj.getNu());
				numericType changeOutlierWeight = (1/contextObj.getTau2()) * std::exp(-0.5 * tmp * tmp / contextObj.getTau2Squared()) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
				fearnheadFilterOptimisedParticle changeOutlierParticle;
				changeOutlierParticle.weight = (outlierSumWeights * contextObj.getOutlierClusterProbability() + notOutlierSumWeights * contextObj.getOutlierProbability()) * contextObj.getChangeProbability() * changeOutlierWeight;

				std::vector<numericType>& destChangeRow = changeProbabilitiesPossibilities[2 * particles.size()];
				std::vector<numericType>& destOutlierRow = outlierProbabilitiesPossibilities[2 * particles.size()];

				std::copy(changeSums.begin(), changeSums.end(), destChangeRow.begin());
				destChangeRow[time] = 1.0;

				std::copy(outlierSums.begin(), outlierSums.end(), destOutlierRow.begin());
				destOutlierRow[time] = 1.0;

				changeOutlierParticle.timeLastChange = time;
				changeOutlierParticle.mean = contextObj.getMu();
				changeOutlierParticle.variance = contextObj.getSigmaSquared();
				childParticles.emplace_back(std::move(changeOutlierParticle));
			}
			//The *single* particle in the next generation corresponding to a change (but not an outlier)
			{
				numericType tmp = (data[time] - contextObj.getMu());
				numericType obsSd = sqrt(contextObj.getSigmaSquared() + contextObj.getTau1Squared());
				numericType changeWeight = (1/obsSd) * std::exp(-0.5 * tmp * tmp / (contextObj.getSigmaSquared() + contextObj.getTau1Squared())) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
				fearnheadFilterOptimisedParticle changeParticle;
				changeParticle.weight = (outlierSumWeights * (1 - contextObj.getOutlierClusterProbability()) + notOutlierSumWeights * (1 - contextObj.getOutlierProbability())) * contextObj.getChangeProbability() * changeWeight;

				std::vector<numericType>& destChangeRow = changeProbabilitiesPossibilities[2 * particles.size() + 1];
				std::vector<numericType>& destOutlierRow = outlierProbabilitiesPossibilities[2 * particles.size() + 1];

				std::copy(changeSums.begin(), changeSums.end(), destChangeRow.begin());
				destChangeRow[time] = 1.0;

				std::copy(outlierSums.begin(), outlierSums.end(), destOutlierRow.begin());
				destOutlierRow[time] = 0.0;

				changeParticle.timeLastChange = time;
				changeParticle.mean = contextObj.getMu();
				changeParticle.variance = contextObj.getSigmaSquared();
				childParticles.emplace_back(std::move(changeParticle));
			}
			if(time == data.size()-1)
			{
				particles.swap(childParticles);
				changeProbabilities.swap(changeProbabilitiesPossibilities);
				outlierProbabilities.swap(outlierProbabilitiesPossibilities);
				break;
			}
			if(childParticles.size() <= args.nParticles)
			{
				particles.swap(childParticles);
				changeProbabilities.swap(changeProbabilitiesPossibilities);
				outlierProbabilities.swap(outlierProbabilitiesPossibilities);
			}
			else
			{
				//Resampling according to the method outlined in Fearnhead. 
				samplingArgs.weights.clear();
				for(int i = 0; i < (int)childParticles.size(); i++) samplingArgs.weights.push_back(childParticles[i].weight);

				samplingArgs.n = args.nParticles;
				bool successful = false;
				try
				{
					samplingDouble::fearnheadSampling(samplingArgs, args.randomSource);
					successful = true;
				}
				catch(...)
				{
				}
				if(samplingArgs.c > 1e10)
				{
					successful = false;
				}
				//the sampling can fail for numerical reasons. E.g. we're unable to get an accurate value for c / c^-1. In that case, assume we just sample the top 1000 particles. 
				if(!successful)
				{
					std::vector<std::pair<int, float> > sorted;
					for(std::size_t i = 0; i < samplingArgs.weights.size(); i++)
					{
						sorted.push_back(std::make_pair((int)i, samplingArgs.weights[i]));
					}
					std::sort(sorted.begin(), sorted.end(), [](const std::pair<int, float>& first, const std::pair<int, float>& second)
						{
							return first.second > second.second;
						});
					samplingArgs.indices.clear();
					for(int i = 0; i < (int)args.nParticles; i++)
					{
						samplingArgs.indices.push_back(sorted[i].first);
					}
					samplingArgs.c = 1 / sorted[0].second;
					samplingArgs.deterministicInclusion.resize(args.nParticles);
					std::fill(samplingArgs.deterministicInclusion.begin(), samplingArgs.deterministicInclusion.end(), false);
				}

				particles.clear();
				for(int i = 0; i < (int)samplingArgs.indices.size(); i++)
				{
					int childIndex = samplingArgs.indices[i];
					particles.emplace_back(std::move(childParticles[childIndex]));
					
					std::vector<numericType>& sourceChange = changeProbabilitiesPossibilities[childIndex];
					std::vector<numericType>& destChange = changeProbabilities[i];
					destChange.swap(sourceChange);

					std::vector<numericType>& sourceOutlier = outlierProbabilitiesPossibilities[childIndex];
					std::vector<numericType>& destOutlier = outlierProbabilities[i];
					destOutlier.swap(sourceOutlier);
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
		args.outlierProbabilities.clear();
		numericType sumWeights = 0;
		for(std::size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++) sumWeights += particles[particleIndex].weight;
		for(std::size_t time = 0; time < data.size(); time++)
		{
			numericType sumOutlier = 0, sumChange = 0;
			for(std::size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
			{
				sumChange += changeProbabilities(particleIndex, time) * particles[particleIndex].weight;
				sumOutlier += outlierProbabilities(particleIndex, time) * particles[particleIndex].weight;
			}
			args.changeProbabilities.push_back(sumChange / sumWeights);
			args.outlierProbabilities.push_back(sumOutlier / sumWeights);
		}
	}
}

