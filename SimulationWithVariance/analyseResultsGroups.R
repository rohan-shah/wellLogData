library(wellLogData)
load("./results.RData")

normalisingConstant <- mean(do.call(c, lapply(results, function(x) x$normalisingConstant)))
numerator <- mean(do.call(c, lapply(results, function(x) x$changeEstimateNumerators[16])))

numerator / normalisingConstant
mean(numeratorGroups / normalisingConstantGroups)
var(as.numeric(numeratorGroups / normalisingConstantGroups))

#Ok, that didn't work. What about taking groups of 10 at a time?
groupSize <- 10
indices <- split(seq_along(results), ceiling(seq_along(results) / groupSize))
empiricalVariances <- horvitzThompsonVariances <- vector(mode = "numeric", length = length(indices) - 1)
for(sampleSize in 2:length(indices))
{
	numeratorGroups <- do.call(c, lapply(indices[1:sampleSize], function(i)
		{
			mean(do.call(c, lapply(results[i], function(x) x$changeEstimateNumerators[16])))
		}))
	secondMomentNumeratorGroups <- do.call(c, lapply(indices[1:sampleSize], function(i)
		{
			firstMoments <- do.call(c, lapply(results[i], function(x) x$changeEstimateNumerators[16]))
			secondMoments <- do.call(c, lapply(results[i], function(x) x$changeEstimateNumeratorSecondMoments))
			m <- outer(firstMoments, firstMoments)
			diag(m) <- secondMoments
			return(sum(m) / (groupSize^2))
		}))
	normalisingConstantGroups <- do.call(c, lapply(indices[1:sampleSize], function(i)
		{
			mean(do.call(c, lapply(results[i], function(x) x$normalisingConstant)))
		}))
	secondMomentNormalisingConstantGroups <- do.call(c, lapply(indices[1:sampleSize], function(i)
		{
			firstMoments <- do.call(c, lapply(results[i], function(x) x$normalisingConstant))
			secondMoments <- do.call(c, lapply(results[i], function(x) x$secondMomentNormalisingConstant))
			m <- outer(firstMoments, firstMoments)
			diag(m) <- secondMoments
			return(sum(m) / (groupSize^2))
		}))
	productExpectationGroups <- do.call(c, lapply(indices[1:sampleSize], function(i)
		{
			normalisingConstants <- do.call(c, lapply(results[i], function(x) x$normalisingConstant))
			numerators <- do.call(c, lapply(results[i], function(x) x$changeEstimateNumerators[16]))
			productTerms <- do.call(c, lapply(results[i], function(x) x$changeEstimateProductExpectations))
			m <- outer(numerators, normalisingConstants)
			diag(m) <- productTerms
			return(sum(m) / (groupSize^2))
		}))
	#The approximation
	horvitzThompsonVariances[sampleSize-1] <- as.numeric((mean(secondMomentNumeratorGroups) * mean(normalisingConstantGroups)^2 - 2 * mean(productExpectationGroups) * mean(normalisingConstantGroups) * mean(numeratorGroups) + mean(secondMomentNormalisingConstantGroups) * mean(numeratorGroups)^2) / (mean(normalisingConstantGroups)^4))
	#The purely empirical 
	empiricalVariances[sampleSize - 1] <- as.numeric((mean(numeratorGroups^2) * mean(normalisingConstantGroups)^2 - 2 * mean(numeratorGroups*normalisingConstantGroups) * mean(normalisingConstantGroups) * mean(numeratorGroups) + (mean(normalisingConstantGroups^2) * mean(numeratorGroups)^2))/ (mean(normalisingConstantGroups)^4))
	cat("Done ", sampleSize - 1, " / ", length(indices) - 1, "\n", sep="")
}
