load("./results.RData")
estimates <- lapply(results, function(x) x$changeProbabilities)
estimates <- do.call(rbind, estimates)

estimatedVariances <- lapply(results, function(x) x$estimatedVariance)
estimatedVariances <- do.call(rbind, estimatedVariances)

htEstimates <- apply(estimatedVariances, 2, mean)
empiricalEstimates <- apply(estimates, 2, var)

numeratorVariances <- do.call(rbind, lapply(results, function(x) x$numeratorVariances))
normalisingConstants <- unlist(lapply(results, function(x) x$normalisingConstant))

scientific_10 <- function(x) 
{
	parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}
library(ggplot2)
plotData <- data.frame("HorvitzThompson" = htEstimates, "Empirical" = empiricalEstimates)
pdf("./wellLogPlot.pdf")
	print(ggplot(data = plotData, mapping = aes(HorvitzThompson, Empirical)) + geom_point(size = 2) + scale_x_log10(labels = scientific_10) + scale_y_log10(labels = scientific_10) + xlab("Horvitz Thompson variance") + ylab("Empirical variance") + theme_bw() + theme(axis.text.x = element_text(size = rel(1.5)), axis.text.y = element_text(size = rel(1.5)), axis.title.x = element_text(size = rel(1.5)), axis.title.y = element_text(size = rel(1.5))) + geom_abline(slope = 1, intercept = 0))
dev.off()

spacing <- seq(10, nrow(estimates), by = 10)
pdf("./tracePlot.pdf")
	for(time in 1:31)
	{
		runningEmpiricalVariance1 <- sapply(spacing, function(x) var(estimates[1:x,time]))
		runningEmpiricalVariance2 <- sapply(spacing, function(x) var(estimates[1:x,time] * normalisingConstants[1:x]) / mean(normalisingConstants[1:x])^2)
		runningHtEstimates <- sapply(spacing, function(x) mean(numeratorVariances[1:x,time]) / mean(normalisingConstants[1:x])^2)
		data <- rbind(data.frame(value = runningEmpiricalVariance1, Method = "Empirical1", count = spacing), data.frame(value = runningEmpiricalVariance2, Method = "Empirical2", count = spacing), data.frame(Method = "Horvitz Thompson", value = runningHtEstimates, count = spacing))
		print(ggplot(data, mapping = aes(count, value, color = Method)) + geom_line() + ggtitle(paste0("Time = ", time, ", estimated value = ", mean(estimates[,time]))))
	}
dev.off()
