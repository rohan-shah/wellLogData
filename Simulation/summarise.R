load("results.RData")

outlierVarianceRatio <- apply(withoutReplacementOutlierResults, 1, var) / apply(fearnheadOutlierResults, 1, var)
changeVarianceRatio <- apply(withoutReplacementChangeResults, 1, var) / apply(fearnheadChangeResults, 1, var)

sum(outlierVarianceRatio > 1, na.rm=TRUE)
sum(outlierVarianceRatio < 1, na.rm=TRUE)

sum(changeVarianceRatio > 1, na.rm=TRUE)
sum(changeVarianceRatio < 1, na.rm=TRUE)

outlierVarianceRatio[outlierVarianceRatio > 100 | outlierVarianceRatio < 0.01] <- NA
changeVarianceRatio[changeVarianceRatio > 100 | changeVarianceRatio < 0.01] <- NA

pdf("./outlierVarianceRatio.pdf")
hist(log10(outlierVarianceRatio), breaks = 100)
dev.off()

pdf("./changeVarianceRatio.pdf")
hist(log10(changeVarianceRatio), breaks = 100)
dev.off()

