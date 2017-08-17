withoutReplacementWithVariance <- function(data, mu, sigma, nu, tau1, tau2, changeProb, outlierProb, outlierClusterProb, nParticles, seed)
{
	result <- .Call("withoutReplacementWithVariance", data, mu, sigma, nu, tau1, tau2, changeProb, outlierProb, outlierClusterProb, nParticles, seed, PACKAGE="wellLogData")
	return(lapply(result, mpfr))
}
