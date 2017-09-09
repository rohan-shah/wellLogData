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
