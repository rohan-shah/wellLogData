library(wellLogData)
load("./results.RData")

normalisingConstant <- mean(do.call(c, lapply(results, function(x) x$normalisingConstant)))
numerator <- mean(do.call(c, lapply(results, function(x) x$changeEstimateNumerators[16])))

numerator / normalisingConstant

groupSize <- 10
indices <- split(seq_along(results), ceiling(seq_along(results) / groupSize))
optionalTransform <- identity#as.numeric
changeEstimateNumerators <- optionalTransform(do.call(c, lapply(results, function(x) x$changeEstimateNumerators[16])))
normalisingConstant <- optionalTransform(do.call(c, lapply(results, function(x) x$normalisingConstant)))
numeratorGroups <- do.call(c, lapply(indices, function(i)
	{
		mean(changeEstimateNumerators[i])
	}))
normalisingConstantGroups <- do.call(c, lapply(indices, function(i)
	{
		mean(normalisingConstant[i])
	}))
mean(numeratorGroups / normalisingConstantGroups)
var(as.numeric(numeratorGroups / normalisingConstantGroups))

empiricalVariances <- horvitzThompsonVariances <- horvitzThompsonVariances2 <- vector(mode = "numeric", length = length(indices) - 1)

changeEstimateNumeratorSecondMoments <- optionalTransform(do.call(c, lapply(results, function(x) x$changeEstimateNumeratorSecondMoments)))
secondMomentNormalisingConstant <- optionalTransform(do.call(c, lapply(results, function(x) x$secondMomentNormalisingConstant)))
changeEstimateProductExpectations <- optionalTransform(do.call(c, lapply(results, function(x) x$changeEstimateProductExpectations)))
changeEstimateNumeratorVariances <- optionalTransform(do.call(c, lapply(results, function(x) x$changeEstimateNumeratorVariances)))
mean(changeEstimateNumeratorVariances / groupSize)
var(as.numeric(numeratorGroups))

for(sampleSize in c(2:10,length(indices)))#2:length(indices))
{
	numeratorGroups <- do.call(c, lapply(indices[1:sampleSize], function(i)
		{
			mean(changeEstimateNumerators[i])
		}))
	numeratorVarianceGroups <- do.call(c, lapply(indices[1:sampleSize], function(i)
		{
			return(mean(changeEstimateNumeratorVariances[i])/groupSize)
		}))
	secondMomentNumeratorGroups <- do.call(c, lapply(indices[1:sampleSize], function(i)
		{
			firstMoments <- changeEstimateNumerators[i]
			secondMoments <- changeEstimateNumeratorSecondMoments[i]
			m <- outer(firstMoments, firstMoments)
			diag(m) <- secondMoments
			return(sum(m) / (groupSize^2))
		}))
	normalisingConstantGroups <- do.call(c, lapply(indices[1:sampleSize], function(i)
		{
			mean(normalisingConstant[i])
		}))
	secondMomentNormalisingConstantGroups <- do.call(c, lapply(indices[1:sampleSize], function(i)
		{
			firstMoments <- normalisingConstant[i]
			secondMoments <- secondMomentNormalisingConstant[i]
			m <- outer(firstMoments, firstMoments)
			diag(m) <- secondMoments
			return(sum(m) / (groupSize^2))
		}))
	productExpectationGroups <- do.call(c, lapply(indices[1:sampleSize], function(i)
		{
			normalisingConstants <- normalisingConstant[i]
			numerators <- changeEstimateNumerators[i]
			productTerms <- changeEstimateProductExpectations[i]
			m <- outer(numerators, normalisingConstants)
			diag(m) <- productTerms
			return(sum(m) / (groupSize^2))
		}))
	#The approximation
	horvitzThompsonVariances[sampleSize-1] <- as.numeric((mean(secondMomentNumeratorGroups) / (mean(numeratorGroups)^2) - (2 * mean(productExpectationGroups)/(mean(numeratorGroups) * mean(normalisingConstantGroups))) + mean(secondMomentNormalisingConstantGroups) / (mean(normalisingConstantGroups)^2)) * (mean(numeratorGroups)^2 / mean(normalisingConstantGroups)^2))
	horvitzThompsonVariances2[sampleSize-1] <- as.numeric(((mean(numeratorVarianceGroups) / (mean(numeratorGroups)^2)) + 1 - (2 * mean(productExpectationGroups)/(mean(numeratorGroups) * mean(normalisingConstantGroups))) + mean(secondMomentNormalisingConstantGroups) / (mean(normalisingConstantGroups)^2)) * (mean(numeratorGroups)^2 / mean(normalisingConstantGroups)^2))
	#The purely empirical 
	empiricalVariances[sampleSize - 1] <- as.numeric((mean(numeratorGroups^2) / (mean(numeratorGroups)^2) - (2 * mean(numeratorGroups * normalisingConstantGroups)/(mean(numeratorGroups) * mean(normalisingConstantGroups))) + mean(normalisingConstantGroups^2) / (mean(normalisingConstantGroups)^2)) * (mean(numeratorGroups)^2 / mean(normalisingConstantGroups)^2))
	cat("Done ", sampleSize - 1, " / ", length(indices) - 1, "\n", sep="")
}
