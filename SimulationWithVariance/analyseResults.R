library(wellLogData)
load("./results.RData")
secondMomentNumerator <- mean(do.call(c, lapply(results, function(x) x$changeEstimateNumeratorSecondMoments)))
secondMomentNumerator
mean(do.call(c, lapply(results, function(x) x$changeEstimateNumerators[16]^2)))

productExpectation <- mean(do.call(c, lapply(results, function(x) x$changeEstimateProductExpectations)))
productExpectation
mean(do.call(c, lapply(results, function(x) x$changeEstimateNumerators[16] * x$normalisingConstant)))
covarianceTerm <- cov(as.numeric(do.call(c, lapply(results, function(x) x$changeEstimateNumerators[16]))), as.numeric(do.call(c, lapply(results, function(x)x$normalisingConstant))))

secondMomentNormalisingConstant <- mean(do.call(c, lapply(results, function(x) x$secondMomentNormalisingConstant)))
secondMomentNormalisingConstant
mean(do.call(c, lapply(results, function(x) x$normalisingConstant^2)))

normalisingConstant <- mean(do.call(c, lapply(results, function(x) x$normalisingConstant)))
varNormalisingConstant <- var(as.numeric(do.call(c, lapply(results, function(x) x$normalisingConstant))))

numerator <- mean(do.call(c, lapply(results, function(x) x$changeEstimateNumerators[16])))
varNumerator <- var(as.numeric(do.call(c, lapply(results, function(x) x$changeEstimateNumerators[16]))))

estimates <- do.call(c, lapply(results, function(x) x$changeEstimateNumerators[16] / x$normalisingConstant))

#Actual estimate of expectation
as.numeric(mean(estimates))
#Approximations to expectation
numerator / normalisingConstant
(numerator / normalisingConstant) - (covarianceTerm / (normalisingConstant^2)) + ((varNormalisingConstant * numerator) / (normalisingConstant^3))

#Actual variance
var(as.numeric(estimates))
#Approximations to variance
(secondMomentNumerator * normalisingConstant^2 - 2 * productExpectation * normalisingConstant * numerator + secondMomentNormalisingConstant * numerator^2) / (normalisingConstant^4)
((numerator / normalisingConstant)^2) * ((varNumerator / (numerator^2)) - 2 * (covarianceTerm / (normalisingConstant * numerator)) + (varNormalisingConstant / (normalisingConstant^2)))

stepSize <- 10
ends <- seq(stepSize, length(results), by = stepSize)
normalisingConstantEmpiricalVariance <- normalisingConstantHorvitzThompsonVariance <- numeratorHorvitzThompsonVariance <- numeratorEmpiricalVariance <- vector(mode = "numeric", length = length(ends))
secondMomentNormalisingConstants <- do.call(c, lapply(results, function(x) x$secondMomentNormalisingConstant))
numeratorVariances <- do.call(c, lapply(results, function(x) x$changeEstimateNumeratorVariance))
normalisingConstants <- do.call(c, lapply(results, function(x) x$normalisingConstant))
numerators <- do.call(c, lapply(results, function(x) x$changeEstimateNumerators[16]))
for(i in 1:length(ends))
{
	normalisingConstantEmpiricalVariance[i] <- var(as.numeric(normalisingConstants[1:ends[i]]))
	numeratorEmpiricalVariance[i] <- var(as.numeric(numerators[1:ends[i]]))
	normalisingConstantHorvitzThompsonVariance[i] <- as.numeric(mean(secondMomentNormalisingConstants[1:ends[i]]) - mean(normalisingConstants[1:ends[i]])^2)
	numeratorHorvitzThompsonVariance[i] <- as.numeric(mean(numeratorVariances[1:ends[i]]))
}
