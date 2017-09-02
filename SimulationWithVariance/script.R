library(wellLogData)
data <- read.table("../well.dat", skip = 7)
n <- 1000
if(file.exists("./results.RData"))
{
	load("./results.RData")
	start <- length(results) + 1
} else
{
	results <- list()
	start <- 1
}
for(i in start:n)
{
	results[[i]] <- withoutReplacementWithVariance(as.numeric(data[1020:1050,1]), mu = 115000, sigma = 10000, nu = 85000, tau1 = 2500, tau2 = 12500, changeProb = 1/250, outlierProb = 0.004, outlierClusterProb = 0.75, nParticles = 45, seed = i)
	if(i %% 10 == 0) save(results, file = "./results.RData")
	cat(i, " / ", n, "\n", sep="")
}
save(results, file = "./results.RData")
