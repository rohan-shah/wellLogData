#include "fearnheadFilterOptimised.h"
#include "fearnheadGetKappa.h"
#include <boost/random/bernoulli_distribution.hpp>
#include "fearnheadSampling.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
namespace wellLogData
{
	struct fearnheadFilterOptimisedParticle
	{
	public:
		//mean and variance of X_t
		double mean, variance;
		std::size_t timeLastChange;
		double weight;
		fearnheadFilterOptimisedParticle()
			: mean(std::numeric_limits<double>::quiet_NaN()), variance(std::numeric_limits<double>::quiet_NaN()), timeLastChange(0), weight(std::numeric_limits<double>::quiet_NaN())
		{}
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

		const std::vector<double>& data = contextObj.getData();

		typedef boost::numeric::ublas::matrix<double, boost::numeric::ublas::row_major> probabilityMatrix;
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
			outlierProbabilities(i, 0) = (double)i;
			particle.timeLastChange = 0;
			particle.mean = contextObj.getMu();
			particle.variance = contextObj.getSigmaSquared();
			if(outlierProbabilities(i, 0) != 0)
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
		for (std::size_t time = 1; time < data.size(); time++)
		{
			//For each particle, make four successors
			childParticles.clear();
			changeSums.resize(time);
			outlierSums.resize(time);
			std::fill(changeSums.begin(), changeSums.end(), 0);
			std::fill(outlierSums.begin(), outlierSums.end(), 0);
			double sumWeights = 0, outlierSumWeights = 0, notOutlierSumWeights = 0;
			for (std::size_t i = 0; i < particles.size(); i++)
			{
				fearnheadFilterOptimisedParticle& currentParticle = particles[i];
				auto currentOutlierRow = boost::numeric::ublas::row(outlierProbabilities, i);
				auto currentChangeRow = boost::numeric::ublas::row(changeProbabilities, i);
#ifndef NDBEUG
				if (currentOutlierRow(time - 1) != 0 && currentOutlierRow(time - 1) != 1) throw std::runtime_error("Internal error");
#endif
				for (std::size_t time2 = 0; time2 < time; time2++)
				{
					outlierSums[time2] += currentOutlierRow(time2) * currentParticle.weight;
					changeSums[time2] += currentChangeRow(time2) * currentParticle.weight;
				}
				sumWeights += currentParticle.weight;
				if (currentOutlierRow(time - 1) > 0) outlierSumWeights += currentParticle.weight;
				else notOutlierSumWeights += currentParticle.weight;

				fearnheadFilterOptimisedParticle childParticleOutlier, childParticleNotOutlier;
				{
					childParticleNotOutlier.weight = currentParticle.weight;

					auto destChangeRow = boost::numeric::ublas::row(changeProbabilitiesPossibilities, 2 * i + 1);
					auto destOutlierRow = boost::numeric::ublas::row(outlierProbabilitiesPossibilities, 2 * i + 1);

					std::copy(currentChangeRow.begin(), currentChangeRow.begin() + time, destChangeRow.begin());
					std::copy(currentOutlierRow.begin(), currentOutlierRow.begin() + time, destOutlierRow.begin());

					if (currentOutlierRow(time - 1) > 0)
					{
						childParticleNotOutlier.weight *= (1 - contextObj.getOutlierClusterProbability());
					}
					else
					{
						childParticleNotOutlier.weight *= (1 - contextObj.getOutlierProbability());
					}
					childParticleNotOutlier.timeLastChange = currentParticle.timeLastChange;
					childParticleNotOutlier.weight *= (1 - contextObj.getChangeProbability());
					if (destOutlierRow(time - 1) > 0.0)
					{
						childParticleNotOutlier.mean = currentParticle.mean;
						childParticleNotOutlier.variance = currentParticle.variance;
					}
					else
					{
						childParticleNotOutlier.mean = ((currentParticle.mean / currentParticle.variance) + (data[time - 1] / contextObj.getTau1Squared())) / ((1 / currentParticle.variance) + (1 / contextObj.getTau1Squared()));
						childParticleNotOutlier.variance = 1 / ((1 / currentParticle.variance) + (1 / contextObj.getTau1Squared()));
					}
					double tmp = data[time] - childParticleNotOutlier.mean;
					double obsSd = sqrt(childParticleNotOutlier.variance + contextObj.getTau1Squared());
					childParticleNotOutlier.weight *= (1 / sqrt(obsSd)) * std::exp(-0.5 * tmp * tmp / (childParticleNotOutlier.variance + contextObj.getTau1Squared())) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
				}
				{
					childParticleOutlier.weight = currentParticle.weight;

					auto destChangeRow = boost::numeric::ublas::row(changeProbabilitiesPossibilities, 2 * i + 0);
					auto destOutlierRow = boost::numeric::ublas::row(outlierProbabilitiesPossibilities, 2 * i + 0);

					std::copy(currentChangeRow.begin(), currentChangeRow.begin() + time, destChangeRow.begin());
					std::copy(currentOutlierRow.begin(), currentOutlierRow.begin() + time, destOutlierRow.begin());
					destOutlierRow(time) = 1.0;

					if (currentOutlierRow(time - 1) > 0)
					{
						childParticleOutlier.weight *= contextObj.getOutlierClusterProbability();
					}
					else
					{
						childParticleOutlier.weight *= contextObj.getOutlierProbability();
					}
					childParticleOutlier.timeLastChange = currentParticle.timeLastChange;
					childParticleOutlier.weight *= (1 - contextObj.getChangeProbability());
					if (destOutlierRow(time - 1) > 0.0)
					{
						childParticleOutlier.mean = currentParticle.mean;
						childParticleOutlier.variance = currentParticle.variance;
					}
					else
					{
						childParticleOutlier.mean = ((currentParticle.mean / currentParticle.variance) + (data[time - 1] / contextObj.getTau1Squared())) / ((1 / currentParticle.variance) + (1 / contextObj.getTau1Squared()));
						childParticleOutlier.variance = 1 / ((1 / currentParticle.variance) + (1 / contextObj.getTau1Squared()));
					}
					double tmp = data[time] - contextObj.getNu();
					childParticleOutlier.weight *= (1 / contextObj.getTau2()) * std::exp(-0.5 * tmp * tmp / contextObj.getTau2Squared()) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
				}
				childParticles.emplace_back(std::move(childParticleOutlier));
				childParticles.emplace_back(std::move(childParticleNotOutlier));
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
				fearnheadFilterOptimisedParticle changeOutlierParticle;
				changeOutlierParticle.weight = (outlierSumWeights * contextObj.getOutlierClusterProbability() + notOutlierSumWeights * contextObj.getOutlierProbability()) * contextObj.getChangeProbability() * changeOutlierWeight;

				auto destChangeRow = boost::numeric::ublas::row(changeProbabilitiesPossibilities, 2 * particles.size());
				auto destOutlierRow = boost::numeric::ublas::row(outlierProbabilitiesPossibilities, 2 * particles.size());

				std::copy(changeSums.begin(), changeSums.end(), destChangeRow.begin());
				destChangeRow(time) = 1.0;

				std::copy(outlierSums.begin(), outlierSums.end(), destOutlierRow.begin());
				destOutlierRow(time) = 1.0;

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
				fearnheadFilterOptimisedParticle changeParticle;
				changeParticle.weight = (outlierSumWeights * (1 - contextObj.getOutlierClusterProbability()) + notOutlierSumWeights * (1 - contextObj.getOutlierProbability())) * contextObj.getChangeProbability() * changeWeight;

				auto destChangeRow = boost::numeric::ublas::row(changeProbabilitiesPossibilities, 2 * particles.size() + 1);
				auto destOutlierRow = boost::numeric::ublas::row(outlierProbabilitiesPossibilities, 2 * particles.size() + 1);

				std::copy(changeSums.begin(), changeSums.end(), destChangeRow.begin());
				destChangeRow(time) = 1.0;

				std::copy(outlierSums.begin(), outlierSums.end(), destOutlierRow.begin());
				destOutlierRow(time) = 0.0;

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
				samplingDouble::fearnheadSampling(samplingArgs, args.randomSource);

				particles.clear();
				for(int i = 0; i < (int)samplingArgs.indices.size(); i++)
				{
					int childIndex = samplingArgs.indices[i];
					particles.emplace_back(std::move(childParticles[childIndex]));
					
					auto sourceChange = boost::numeric::ublas::row(changeProbabilitiesPossibilities, childIndex);
					auto destChange = boost::numeric::ublas::row(changeProbabilities, i);
					std::copy(sourceChange.begin(), sourceChange.end(), destChange.begin());

					auto sourceOutlier = boost::numeric::ublas::row(outlierProbabilitiesPossibilities, childIndex);
					auto destOutlier = boost::numeric::ublas::row(outlierProbabilities, i);
					std::copy(sourceOutlier.begin(), sourceOutlier.end(), destOutlier.begin());
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
		double sumWeights = 0;
		for(std::size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++) sumWeights += particles[particleIndex].weight;
		for(std::size_t time = 0; time < data.size(); time++)
		{
			double sumOutlier = 0, sumChange = 0;
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

