#' xcms version 3.16.1 
library(xcms)

#' Subset of MT01 data
files <- c("MT_20210803_005-NIST.mzML", "MT_20210803_051.mzML", "MT_20210803_089-Qstd.mzML", "MT_20210803_091.mzML",     
 "MT_20210803_139.mzML", "MT_20210803_181.mzML")

pd <- data.frame(sample_name = files)
raw_data <- readMSData(files = paste('MT01', files, sep='/'), pdata = new("NAnnotatedDataFrame", pd), mode = "onDisk")

cwp <- CentWaveParam(peakwidth = c(1, 30), ppm = 5,  noise = 1000, prefilter = c(3, 1000))
xdata <- findChromPeaks(raw_data, param = cwp)

#' for merging peaks
mpp <- MergeNeighboringPeaksParam(expandRt = 5, ppm = 1, minProp = 0.5)
xdata <- refineChromPeaks(xdata, mpp)

#' RT alignment
xdata <- adjustRtime(xdata, param = ObiwarpParam())

pdp <- PeakDensityParam(sampleGroups = rep(1,length(files)), minFraction = 0.1, bw = 3, binSize=0.001)
xdata <- groupChromPeaks(xdata, param = pdp)

ftDef <- featureDefinitions(xdata)
ftValues <- featureValues(xdata, value = "into", method = "sum")

write.csv(cbind(ftDef[,1:7],ftValues), "subsetMT01_XCMS_featureTable.csv")
