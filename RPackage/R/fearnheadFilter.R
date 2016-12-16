fearnheadFilter <- function(data, mu, sigma, nu, tau1, tau2, changeProb, outlierProb, outlierClusterProb, nParticles, seed)
{
	return(.Call("fearnheadFilter", data, mu, sigma, nu, tau1, tau2, changeProb, outlierProb, outlierClusterProb, nParticles, seed, PACKAGE="wellLogData"))
}
