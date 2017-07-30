library(wellLogData)
nreps <- 1000
wellData <- read.table("~/Software/wellLogData/well.dat", row.names=NULL, colClasses = "numeric")
if(file.exists("./tmp.RData"))
{
	load("./tmp.RData")
	start <- min(which(apply(withoutReplacementChangeResults, 2, function(x) all(x == 0))))
} else
{
	withoutReplacementChangeResults <- fearnheadChangeResults <- withoutReplacementOutlierResults <- fearnheadOutlierResults <- matrix(0, nrow = length(wellData[,1]), ncol = nreps)
	start <- 1
}
for(i in start:nreps)
{
	result <- fearnheadFilter(data = as.numeric(wellData[,1]), mu = 115000, sigma = 10000, nu = 85000, tau1 = 2500, tau2 = 12500, changeProb = 1/250, outlierProb = 0.004, outlierClusterProb = 0.75, nParticles = 100, seed = i)
	fearnheadOutlierResults[, i] <- result$changeProbabilities
	fearnheadChangeResults[, i] <- result$outlierProbabilities

	result <- withoutReplacement(data = as.numeric(wellData[,1]), mu = 115000, sigma = 10000, nu = 85000, tau1 = 2500, tau2 = 12500, changeProb = 1/250, outlierProb = 0.004, outlierClusterProb = 0.75, nParticles = 100, seed = i)
	withoutReplacementOutlierResults[, i] <- result$changeProbabilities
	withoutReplacementChangeResults[, i] <- result$outlierProbabilities
	cat(i, " / ", nreps, "\n", sep = "")
	if(i %% 10 == 0) save(withoutReplacementChangeResults, fearnheadChangeResults, withoutReplacementOutlierResults, fearnheadOutlierResults, file = "results.RData")
}
