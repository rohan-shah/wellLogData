withoutReplacement <- function(data, mu, sigma, nu, tau1, tau2, changeProb, outlierProb, outlierClusterProb, nParticles, seed)
{
	return(.Call("withoutReplacement", data, mu, sigma, nu, tau1, tau2, changeProb, outlierProb, outlierClusterProb, nParticles, seed, PACKAGE="wellLogData"))
}
