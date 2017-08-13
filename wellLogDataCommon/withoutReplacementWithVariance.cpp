#include "withoutReplacementWithVariance.h"
#include <boost/random/bernoulli_distribution.hpp>
#include "conditionalPoissonSequential.h"
#include <boost/graph/graphml.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>
namespace wellLogData
{
	struct varianceGraphVertex
	{
	public:
		varianceGraphVertex()
			: indexWithinDesign(-1), indexWithinSelected(-1), samplingStage(-1), accumulatedMean(0), trueDensity(0), V(0)
		{}
		int indexWithinDesign, indexWithinSelected, samplingStage;
		mutable ::sampling::mpfr_class accumulatedMean, trueDensity, V;
	};
	typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS, boost::property<boost::vertex_name_t, varianceGraphVertex>, boost::property<boost::edge_name_t, double> > varianceGraph;
	struct withoutReplacementWithVarianceParticle
	{
	public:
		//mean and variance of X_t
		double mean, variance;
		std::vector<sampling::mpfr_class> isOutlier, isChange;
		std::size_t timeLastChange;
		sampling::mpfr_class weight;
		sampling::mpfr_class trueDensity;
		varianceGraph::vertex_descriptor vertex;
		withoutReplacementWithVarianceParticle()
			: mean(std::numeric_limits<double>::quiet_NaN()), variance(std::numeric_limits<double>::quiet_NaN()), timeLastChange(0), weight(std::numeric_limits<double>::quiet_NaN()), trueDensity(-1), vertex(-1)
		{}
		withoutReplacementWithVarianceParticle& operator=(withoutReplacementWithVarianceParticle&& other)
		{
			mean = other.mean;
			variance = other.variance;
			isOutlier = std::move(other.isOutlier);
			isChange = std::move(other.isChange);
			weight = other.weight;
			trueDensity = other.trueDensity;
			vertex = other.vertex;
			return *this;
		}
		withoutReplacementWithVarianceParticle(withoutReplacementWithVarianceParticle&& other)
			: mean(other.mean), variance(other.variance), isOutlier(std::move(other.isOutlier)), isChange(std::move(other.isChange)), weight(other.weight), trueDensity(other.trueDensity), vertex(other.vertex)
		{}
		withoutReplacementWithVarianceParticle(const withoutReplacementWithVarianceParticle& other)
			: mean(other.mean), variance(other.variance), isOutlier(other.isOutlier), isChange(other.isChange), weight(other.weight), trueDensity(other.trueDensity), vertex(other.vertex)
		{}
	};
	class vertexPropertyWriter
	{
	public:
		vertexPropertyWriter(const varianceGraph& g)
			:g(g)
		{}
		template<typename vertexDesc> void operator()(std::ostream& out, const vertexDesc& v)
		{
			const varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, g, v);
			out << std::endl << "[trueDensity=\"" << vertexInfo.trueDensity.convert_to<double>() << "\"";
			out << ",accumulatedMean=\"" << vertexInfo.accumulatedMean.convert_to<double>() << "\"";
			out << ",indexWithinDesign=\"" << vertexInfo.indexWithinDesign << "\"";
			out << ",samplingStage=\"" << vertexInfo.samplingStage << "\"";
			out << ",V=\"" << vertexInfo.V.convert_to<double>() << "\"";
			out << ",indexWithinSelected=\"" << vertexInfo.indexWithinSelected <<"\"]" << std::endl;
		}
		const varianceGraph& g;
	};
	bool particleWeightSorter(const withoutReplacementWithVarianceParticle& first, const withoutReplacementWithVarianceParticle& second)
	{
		return first.weight < second.weight;
	}
	void withoutReplacementWithVariance(const context& contextObj, withoutReplacementWithVarianceArgs& args)
	{
		if(args.nParticles == 0)
		{
			throw std::runtime_error("Input nParticles cannot be zero");
		}

		const std::vector<double>& data = contextObj.getData();

		sampling::conditionalPoissonSequentialArgs samplingArgs(true);
		args.changeProbabilities.clear();
		args.outlierProbabilities.clear();
		boost::random::bernoulli_distribution<> isOutlierDist(contextObj.getOutlierProbability()), isChangeDist(contextObj.getChangeProbability());

		std::vector<double> systematicWeights;
		std::vector<int> systematicIndices;

		//Initial simulation
		std::vector<withoutReplacementWithVarianceParticle> particles, childParticles;

		std::vector<std::vector<int> > graphVertices(data.size());
		std::vector<std::vector<::sampling::mpfr_class> > allInclusionProbabilities;
		std::vector<boost::numeric::ublas::matrix<sampling::mpfr_class> > allSecondOrderInclusionProbabilities;

		varianceGraph varianceEstimationGraph;
		varianceGraph::vertex_descriptor rootVertex = boost::add_vertex(varianceEstimationGraph);
		//Outier
		{
			withoutReplacementWithVarianceParticle particle;
			particle.isChange.resize(1, 1.0);
			particle.isOutlier.resize(1, 1.0);
			particle.timeLastChange = 0;
			particle.mean = contextObj.getMu();
			particle.variance = contextObj.getSigmaSquared();
			double tmp = (data[0] - contextObj.getNu());
			particle.weight = contextObj.getOutlierProbability() * (1/contextObj.getTau2()) * std::exp(-0.5 * tmp * tmp / contextObj.getTau2Squared()) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
			particle.trueDensity = particle.weight;
			particle.vertex = boost::add_vertex(varianceEstimationGraph);
			boost::add_edge(rootVertex, particle.vertex, particle.weight.convert_to<double>(), varianceEstimationGraph);
			particles.emplace_back(std::move(particle));
			graphVertices[0].push_back((int)particle.vertex);

			varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, particle.vertex);
			vertexInfo.trueDensity = particle.trueDensity;
			vertexInfo.indexWithinDesign = vertexInfo.indexWithinSelected = 0;
			vertexInfo.samplingStage = 0;
		}
		//Non outlier
		{
			withoutReplacementWithVarianceParticle particle;
			particle.isChange.resize(1, 1.0);
			particle.isOutlier.resize(1, 0.0);
			particle.timeLastChange = 0;
			particle.mean = contextObj.getMu();
			particle.variance = contextObj.getSigmaSquared();
			double tmp = (data[0] - particle.mean);
			particle.weight = (1 - contextObj.getOutlierProbability()) * (1/contextObj.getSigma()) * std::exp(-0.5 * tmp * tmp / particle.variance) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
			particle.trueDensity = particle.weight;
			particle.vertex = boost::add_vertex(varianceEstimationGraph);
			boost::add_edge(rootVertex, particle.vertex, particle.weight.convert_to<double>(), varianceEstimationGraph);
			particles.emplace_back(std::move(particle));
			graphVertices[0].push_back((int)particle.vertex);

			varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, particle.vertex);
			vertexInfo.trueDensity = particle.trueDensity;
			vertexInfo.indexWithinDesign = vertexInfo.indexWithinSelected = 1;
			vertexInfo.samplingStage = 0;
		}
		{
			std::vector<::sampling::mpfr_class> inclusion(2, 1.0);
			allInclusionProbabilities.emplace_back(std::move(inclusion));
			boost::numeric::ublas::matrix<sampling::mpfr_class> secondOrder(2, 2, 1);
			allSecondOrderInclusionProbabilities.emplace_back(std::move(secondOrder));
		}

		std::vector<sampling::mpfr_class> changeSums, outlierSums;
		for(std::size_t time = 1; time < data.size(); time++)
		{
			//For each particle, make four successors
			childParticles.clear();
			changeSums.resize(time);
			outlierSums.resize(time);
			std::fill(changeSums.begin(), changeSums.end(), 0);
			std::fill(outlierSums.begin(), outlierSums.end(), 0);
			sampling::mpfr_class sumWeights = 0, outlierSumWeights = 0, notOutlierSumWeights = 0;
			sampling::mpfr_class outlierSumTrueDensity = 0, notOutlierSumTrueDensity = 0;
			samplingArgs.weights.clear();
			for(std::size_t i = 0; i < particles.size(); i++)
			{
				withoutReplacementWithVarianceParticle& currentParticle = particles[i];
#ifndef NDBEUG
				if(currentParticle.isOutlier[time-1] != 0 && currentParticle.isOutlier[time-1] != 1) throw std::runtime_error("Internal error");
#endif
				for(std::size_t time2 = 0; time2 < time; time2++)
				{
					outlierSums[time2] += currentParticle.isOutlier[time2] * currentParticle.weight;
					changeSums[time2] += currentParticle.isChange[time2] * currentParticle.weight;
				}
				sumWeights += currentParticle.weight;
				if(currentParticle.isOutlier[time-1] > 0)
				{
					outlierSumWeights += currentParticle.weight;
					outlierSumTrueDensity += currentParticle.trueDensity;
				}
				else
				{
					notOutlierSumWeights += currentParticle.weight;
					notOutlierSumTrueDensity += currentParticle.trueDensity;
				}
				for(int j = 0; j < 2; j++)
				{
					withoutReplacementWithVarianceParticle childParticle;
					childParticle.vertex = boost::add_vertex(varianceEstimationGraph);
					varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, childParticle.vertex);
					vertexInfo.indexWithinDesign = childParticles.size();

					childParticle.isChange = currentParticle.isChange;
					childParticle.isChange.push_back(0.0);

					childParticle.isOutlier = currentParticle.isOutlier;
					childParticle.isOutlier.push_back((double)(j % 2));

					double extraPart = 1;
					if(currentParticle.isOutlier[time-1] > 0)
					{
						if(j % 2)
						{
							extraPart *= contextObj.getOutlierClusterProbability();
						}
						else
						{
							extraPart *= (1 - contextObj.getOutlierClusterProbability());
						}
					}
					else
					{
						if(j % 2)
						{
							extraPart *= contextObj.getOutlierProbability();
						}
						else
						{
							extraPart *= (1 - contextObj.getOutlierProbability());
						}
					}
					childParticle.timeLastChange = currentParticle.timeLastChange;
					extraPart *= (1 - contextObj.getChangeProbability());
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
						extraPart *= (1/contextObj.getTau2()) * std::exp(-0.5 * tmp * tmp / contextObj.getTau2Squared()) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
					}
					else
					{
						double tmp = data[time] - childParticle.mean;
						double obsSd = sqrt(childParticle.variance + contextObj.getTau1Squared());
						extraPart *= (1/sqrt(obsSd)) * std::exp(-0.5 * tmp * tmp / (childParticle.variance + contextObj.getTau1Squared())) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
					}
					boost::add_edge(currentParticle.vertex, childParticle.vertex, extraPart, varianceEstimationGraph);
					vertexInfo.trueDensity = childParticle.trueDensity = currentParticle.trueDensity * extraPart;
					vertexInfo.samplingStage = time;
					childParticle.weight = currentParticle.weight * extraPart;
					samplingArgs.weights.push_back(childParticle.weight);
					childParticles.emplace_back(std::move(childParticle));
				}
			}
			//The *single* particle in the next generation corresponding to both an outlier *and* a change
			withoutReplacementWithVarianceParticle changeOutlierParticle;
			{
				double tmp = (data[time] - contextObj.getNu());
				double changeOutlierWeight = (1/contextObj.getTau2()) * std::exp(-0.5 * tmp * tmp / contextObj.getTau2Squared()) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
				changeOutlierParticle.weight = (outlierSumWeights * contextObj.getOutlierClusterProbability() + notOutlierSumWeights * contextObj.getOutlierProbability()) * contextObj.getChangeProbability() * changeOutlierWeight;
				changeOutlierParticle.trueDensity = (outlierSumTrueDensity * contextObj.getOutlierClusterProbability() + notOutlierSumTrueDensity * contextObj.getOutlierProbability()) * contextObj.getChangeProbability() * changeOutlierWeight;
				changeOutlierParticle.isChange = changeSums;
				changeOutlierParticle.isChange.push_back(1.0);

				changeOutlierParticle.isOutlier = outlierSums;
				changeOutlierParticle.isOutlier.push_back(1.0);

				changeOutlierParticle.timeLastChange = time;
				changeOutlierParticle.mean = contextObj.getMu();
				changeOutlierParticle.variance = contextObj.getSigmaSquared();
				samplingArgs.weights.push_back(changeOutlierParticle.weight);
				
				changeOutlierParticle.vertex = boost::add_vertex(varianceEstimationGraph);
				for(std::size_t i = 0; i < particles.size(); i++)
				{
					if(particles[i].isOutlier[time-1] == 1.0)
					{
						boost::add_edge(particles[i].vertex, changeOutlierParticle.vertex, contextObj.getOutlierClusterProbability() * contextObj.getChangeProbability() * changeOutlierWeight, varianceEstimationGraph);
					}
					else
					{
						boost::add_edge(particles[i].vertex, changeOutlierParticle.vertex, contextObj.getOutlierProbability() * contextObj.getChangeProbability() * changeOutlierWeight, varianceEstimationGraph);
					}
				}
				varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, changeOutlierParticle.vertex);
				vertexInfo.trueDensity = changeOutlierParticle.trueDensity;
				vertexInfo.indexWithinDesign = (int)childParticles.size();
				vertexInfo.samplingStage = time;
				
				childParticles.emplace_back(std::move(changeOutlierParticle));
			}
			//The *single* particle in the next generation corresponding to a change (but not an outlier)
			withoutReplacementWithVarianceParticle changeParticle;
			{
				double tmp = (data[time] - contextObj.getMu());
				double obsSd = sqrt(contextObj.getSigmaSquared() + contextObj.getTau1Squared());
				double changeWeight = (1/obsSd) * std::exp(-0.5 * tmp * tmp / (contextObj.getSigmaSquared() + contextObj.getTau1Squared())) * /* 1/sqrt(2 * pi) */ M_SQRT1_2 * 0.5 * M_2_SQRTPI;
				changeParticle.weight = (outlierSumWeights * (1 - contextObj.getOutlierClusterProbability()) + notOutlierSumWeights * (1 - contextObj.getOutlierProbability())) * contextObj.getChangeProbability() * changeWeight;
				changeParticle.trueDensity = (outlierSumTrueDensity * (1 - contextObj.getOutlierClusterProbability()) + notOutlierSumTrueDensity * (1 - contextObj.getOutlierProbability())) * contextObj.getChangeProbability() * changeWeight;

				changeParticle.isChange = changeSums;
				changeParticle.isChange.push_back(1.0);

				changeParticle.isOutlier = outlierSums;
				changeParticle.isOutlier.push_back(0.0);

				changeParticle.timeLastChange = time;
				changeParticle.mean = contextObj.getMu();
				changeParticle.variance = contextObj.getSigmaSquared();
				samplingArgs.weights.push_back(changeParticle.weight);
				
				changeParticle.vertex = boost::add_vertex(varianceEstimationGraph);
				for(std::size_t i = 0; i < particles.size(); i++)
				{
					if(particles[i].isOutlier[time-1] == 1.0)
					{
						boost::add_edge(particles[i].vertex, changeParticle.vertex, (1 - contextObj.getOutlierClusterProbability()) * contextObj.getChangeProbability() * changeWeight, varianceEstimationGraph);
					}
					else
					{
						boost::add_edge(particles[i].vertex, changeParticle.vertex, (1 - contextObj.getOutlierProbability()) * contextObj.getChangeProbability() * changeWeight, varianceEstimationGraph);
					}
				}
				varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, changeParticle.vertex);
				vertexInfo.trueDensity = changeParticle.trueDensity;
				vertexInfo.indexWithinDesign = (int)childParticles.size();
				vertexInfo.samplingStage = time;

				childParticles.emplace_back(std::move(changeParticle));
			}
			if(childParticles.size() <= args.nParticles)
			{
				for(int i = 0; i < (int)childParticles.size(); i++)
				{
					graphVertices[time].push_back((int)childParticles[i].vertex);
					boost::get(boost::vertex_name, varianceEstimationGraph, childParticles[i].vertex).indexWithinSelected  = i;
				}
				{
					std::vector<::sampling::mpfr_class> inclusion(childParticles.size(), 1.0);
					allInclusionProbabilities.emplace_back(std::move(inclusion));
					boost::numeric::ublas::matrix<sampling::mpfr_class> secondOrder(childParticles.size(), childParticles.size(), 1);
					allSecondOrderInclusionProbabilities.emplace_back(std::move(secondOrder));
				}
				particles.swap(childParticles);
			}
			else
			{
				samplingArgs.n = args.nParticles;
				samplingArgs.indices.clear();
				sampling::mpfr_class maximumWeight = 0;
				for(int i = 0; i < (int)samplingArgs.weights.size(); i++)
				{
					maximumWeight = std::max(maximumWeight, samplingArgs.weights[i]);
				}
				for(int i = 0; i < (int)samplingArgs.weights.size(); i++)
				{
					samplingArgs.weights[i] /= maximumWeight;
					samplingArgs.weights[i] = boost::multiprecision::max(sampling::mpfr_class(1e-5), samplingArgs.weights[i]);
				}
				sampling::conditionalPoissonSequential(samplingArgs, args.randomSource);
				particles.clear();
				for(int i = 0; i < (int)samplingArgs.indices.size(); i++)
				{
					particles.emplace_back(childParticles[samplingArgs.indices[i]]);
					graphVertices[time].push_back((int)particles.back().vertex);
					particles.back().weight /= samplingArgs.inclusionProbabilities[samplingArgs.indices[i]].convert_to<double>();
					boost::get(boost::vertex_name, varianceEstimationGraph, particles.back().vertex).indexWithinSelected  = i;
				}
				{
					boost::numeric::ublas::matrix<sampling::mpfr_class> secondOrder;
					sampling::conditionalPoissonSecondOrderInclusionProbabilities(samplingArgs, samplingArgs.inclusionProbabilities, secondOrder);
					allSecondOrderInclusionProbabilities.emplace_back(std::move(secondOrder));
					allInclusionProbabilities.emplace_back(std::move(samplingArgs.inclusionProbabilities));
				}
			}
		}
		args.changeProbabilities.clear();
		args.outlierProbabilities.clear();
		sampling::mpfr_class normalisingConstant = 0;
		for(std::size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++) normalisingConstant += particles[particleIndex].weight;
		for(std::size_t time = 0; time < data.size(); time++)
		{
			sampling::mpfr_class sumOutlier = 0, sumChange = 0;
			for(std::size_t particleIndex = 0; particleIndex < particles.size(); particleIndex++)
			{
				sumChange += particles[particleIndex].isChange[time] * particles[particleIndex].weight;
				sumOutlier += particles[particleIndex].isOutlier[time] * particles[particleIndex].weight;
			}
			args.changeProbabilities.push_back(sampling::mpfr_class(sumChange / normalisingConstant).convert_to<double>());
			args.outlierProbabilities.push_back(sampling::mpfr_class(sumOutlier / normalisingConstant).convert_to<double>());
		}
		for(int timeForEstimate = 0; timeForEstimate < (int)data.size(); timeForEstimate++)
		{
			//Initialize the graph
			for(int i = 0; i < (int)particles.size(); i++)
			{
				varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, particles[i].vertex);
				vertexInfo.accumulatedMean = particles[i].isChange[timeForEstimate];
			}
			std::vector<boost::numeric::ublas::matrix<sampling::mpfr_class> > allCovariances(data.size());
			allCovariances[data.size()-1].resize(graphVertices[data.size()-1].size(), graphVertices[data.size()-1].size(), false);
			for(int time = data.size() - 2; time >= 0; time--)
			{
				const std::vector<sampling::mpfr_class>& currentInclusionProbabilities = allInclusionProbabilities[time+1];
				const boost::numeric::ublas::matrix<sampling::mpfr_class>& currentSecondOrderInclusionProbabilities = allSecondOrderInclusionProbabilities[time+1];

				allCovariances[time].resize(graphVertices[time].size(), graphVertices[time].size(), false);
				boost::numeric::ublas::matrix<sampling::mpfr_class>& currentCovariance = allCovariances[time];
				const boost::numeric::ublas::matrix<sampling::mpfr_class>& previousCovariance = allCovariances[time+1];
				//First estimate the accumulated means
				for(int particleCounter = 0; particleCounter < (int)graphVertices[time].size(); particleCounter++)
				{
					int currentVertex = graphVertices[time][particleCounter];
					varianceGraphVertex& vertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, currentVertex);
					vertexInfo.accumulatedMean = 0;
					varianceGraph::out_edge_iterator current, end;
					boost::tie(current, end) = boost::out_edges(currentVertex, varianceEstimationGraph);
					for(; current != end; current++)
					{
						int targetVertex = boost::target(*current, varianceEstimationGraph);
						double edgeTrueDensity = boost::get(boost::edge_name, varianceEstimationGraph, *current);
						varianceGraphVertex& targetVertexInfo = boost::get(boost::vertex_name, varianceEstimationGraph, targetVertex);
						vertexInfo.accumulatedMean += (edgeTrueDensity * targetVertexInfo.accumulatedMean / currentInclusionProbabilities[targetVertexInfo.indexWithinDesign]);
					}
				}
				//now actually do the covariance estimation. 
				for(int particleCounter1 = 0; particleCounter1 < (int)graphVertices[time].size(); particleCounter1++)
				{
					int graphVertex1 = graphVertices[time][particleCounter1];
					varianceGraphVertex& graphVertex1Info = boost::get(boost::vertex_name, varianceEstimationGraph, graphVertex1);
					for(int particleCounter2 = 0; particleCounter2 < (int)graphVertices[time].size(); particleCounter2++)
					{
						int graphVertex2 = graphVertices[time][particleCounter2];
						varianceGraphVertex& graphVertex2Info = boost::get(boost::vertex_name, varianceEstimationGraph, graphVertex2);
						sampling::mpfr_class& currentCovarianceValue = currentCovariance(particleCounter1, particleCounter2);
					
						varianceGraph::out_edge_iterator current1, end1, current2, end2;
						boost::tie(current1, end1) = boost::out_edges(graphVertex1, varianceEstimationGraph);
						for(; current1 != end1; current1++)
						{
							int targetVertex1 = boost::target(*current1, varianceEstimationGraph);
							varianceGraphVertex& targetVertexInfo1 = boost::get(boost::vertex_name, varianceEstimationGraph, targetVertex1);
							if(targetVertexInfo1.indexWithinSelected != -1)
							{
								double multiple1 = boost::get(boost::edge_name, varianceEstimationGraph, *current1);
								boost::tie(current2, end2) = boost::out_edges(graphVertex2, varianceEstimationGraph);
								for(; current2 != end2; current2++)
								{
									int targetVertex2 = boost::target(*current2, varianceEstimationGraph);
									varianceGraphVertex& targetVertexInfo2 = boost::get(boost::vertex_name, varianceEstimationGraph, targetVertex2);
									if(targetVertexInfo2.indexWithinSelected != -1)
									{
										sampling::mpfr_class inclusionProduct = currentInclusionProbabilities[targetVertexInfo1.indexWithinDesign] * currentInclusionProbabilities[targetVertexInfo2.indexWithinDesign];
										double multiple2 = boost::get(boost::edge_name, varianceEstimationGraph, *current2);
										if(targetVertex1 == targetVertex2)
										{
											currentCovarianceValue += ((currentInclusionProbabilities[targetVertexInfo1.indexWithinDesign] - inclusionProduct) * targetVertexInfo1.accumulatedMean * targetVertexInfo2.accumulatedMean * graphVertex1Info.trueDensity * graphVertex2Info.trueDensity * multiple1 * multiple2 / (currentInclusionProbabilities[targetVertexInfo1.indexWithinDesign]*inclusionProduct)) + (previousCovariance(targetVertexInfo1.indexWithinSelected, targetVertexInfo2.indexWithinSelected) / (currentInclusionProbabilities[targetVertexInfo1.indexWithinDesign])) * (multiple1 * graphVertex1Info.trueDensity / targetVertexInfo1.trueDensity) * (multiple2 * graphVertex2Info.trueDensity / targetVertexInfo2.trueDensity);
										}
										else
										{
											currentCovarianceValue += ((currentSecondOrderInclusionProbabilities(targetVertexInfo1.indexWithinDesign, targetVertexInfo2.indexWithinDesign) - inclusionProduct) * targetVertexInfo1.accumulatedMean * targetVertexInfo2.accumulatedMean * graphVertex1Info.trueDensity * graphVertex2Info.trueDensity * multiple1 * multiple2 / (currentSecondOrderInclusionProbabilities(targetVertexInfo1.indexWithinDesign, targetVertexInfo2.indexWithinDesign) * inclusionProduct)) + (previousCovariance(targetVertexInfo1.indexWithinSelected, targetVertexInfo2.indexWithinSelected) / (currentSecondOrderInclusionProbabilities(targetVertexInfo1.indexWithinDesign, targetVertexInfo2.indexWithinDesign))) * (multiple1 * graphVertex1Info.trueDensity / targetVertexInfo1.trueDensity) * (multiple2 * graphVertex2Info.trueDensity / targetVertexInfo2.trueDensity);
										}
									}
								}
							}
						}
						if(particleCounter1 == particleCounter2)
						{
							graphVertex1Info.V = currentCovarianceValue;
						}
					}
				}
				//If a variance is zero, the covariances *have* to be zero. 
				for(int particleCounter1 = 0; particleCounter1 < (int)graphVertices[time].size(); particleCounter1++)
				{
					if(currentCovariance(particleCounter1, particleCounter1) == 0)
					{
						for(int particleCounter2 = 0; particleCounter2 < (int)graphVertices[time].size(); particleCounter2++)
						{
							currentCovariance(particleCounter1, particleCounter2) = currentCovariance(particleCounter2, particleCounter1) = 0;
						}
					}
				}
			}
			sampling::mpfr_class totalCovarianceFromGraph = allCovariances[0](0, 0) + allCovariances[0](1, 0) + allCovariances[0](0, 1) + allCovariances[0](1, 1);
			args.changeEstimateNumeratorVariances.push_back(totalCovarianceFromGraph.convert_to<double>());
			args.normalisingConstant = normalisingConstant.convert_to<double>();
		}
		//sampling::mpfr_class graphEstimate = (boost::get(boost::vertex_name, varianceEstimationGraph, 1).accumulatedMean * boost::get(boost::vertex_name, varianceEstimationGraph, 1).trueDensity + boost::get(boost::vertex_name, varianceEstimationGraph, 2).accumulatedMean * boost::get(boost::vertex_name, varianceEstimationGraph, 2).trueDensity) / normalisingConstant;
/*		{
			std::ofstream file("graph.graphml");
			boost::dynamic_properties dp;
			dp.property("samplingStage", boost::make_transform_value_property_map(boost::bind<int&>(&approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex::samplingStage, _1), boost::get(boost::vertex_name_t(), varianceEstimationGraph)));
			dp.property("indexWithinDesign", boost::make_transform_value_property_map(boost::bind<int&>(&approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex::indexWithinDesign, _1), boost::get(boost::vertex_name_t(), varianceEstimationGraph)));
			dp.property("indexWithinSelected", boost::make_transform_value_property_map(boost::bind<int&>(&approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex::indexWithinSelected, _1), boost::get(boost::vertex_name_t(), varianceEstimationGraph)));
			auto vertexToAccumulatedMean = [](approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& x)
				{
					return x.accumulatedMean.convert_to<double>();
				};
			dp.property("accumulatedMean", boost::make_transform_value_property_map(vertexToAccumulatedMean, boost::get(boost::vertex_name_t(), varianceEstimationGraph)));
			auto vertexToTrueDensity = [](approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& x)
				{
					return x.trueDensity.convert_to<double>();
				};
			dp.property("trueDensity", boost::make_transform_value_property_map(vertexToTrueDensity, boost::get(boost::vertex_name_t(), varianceEstimationGraph)));
			auto vertexToV = [](approximateZeroVarianceWORMergeWithVarianceImpl::varianceGraphVertex& x)
				{
					return x.V.convert_to<double>();
				};
			dp.property("V", boost::make_transform_value_property_map(vertexToV, boost::get(boost::vertex_name_t(), varianceEstimationGraph)));
			boost::write_graphml(file, varianceEstimationGraph, dp);
		}*/
		/*{
			std::ofstream file("graph.dot");
			vertexPropertyWriter vp(varianceEstimationGraph);
			//edgePropertyWriter ep(varianceEstimationGraph);
			//boost::write_graphviz(file, varianceEstimationGraph, vp, ep);
			boost::write_graphviz(file, varianceEstimationGraph, vp);
		}*/
	}
}

