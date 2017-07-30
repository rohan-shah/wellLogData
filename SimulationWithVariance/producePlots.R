load("./results.RData")
estimates <- lapply(results, function(x) x$estimatedVariance)
estimates <- do.call(rbind, estimates)

estimatedVariances <- lapply(results, function(x) x$estimatedVariance)
estimatedVariances <- do.call(rbind, estimatedVariances)

htEstimates <- apply(estimatedVariances, 2, mean)
empiricalEstimates <- apply(estimates, 2, var)
