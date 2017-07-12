withoutReplacementWithVariance <- function(data, mu, sigma, nu, tau1, tau2, changeProb, outlierProb, outlierClusterProb, nParticles, seed)
{
	return(.Call("withoutReplacementWithVariance", data, mu, sigma, nu, tau1, tau2, changeProb, outlierProb, outlierClusterProb, nParticles, seed, PACKAGE="wellLogData"))
}
