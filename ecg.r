# libraries
library("data.table")
library("ggpmisc")
library("ggplot2")
library("photobiology")

# files
ecg1 = as.data.frame(data.table::fread("ECG1.csv"))
ecg2 = as.data.frame(data.table::fread("ECG2.csv"))

ecg1 = ecg1[,-13]

########################################################
                   # FUNCTIONS #
########################################################

# takes a dataframe of samples and separates them into a list of each sample
# also converts sample data to be numeric
# samples[[1]][[1]] - sample # (ex. '3560')
# samples[[1]][[2]] - sample data
#                     dataframe of two columns
#                     - 'Time (ms)' - 271.6425, 271.7675, 271.8925...
#                     - 'ECG (mV)' - -69.1224, -73.5474, -77.7436...
listSamples = function(df) {
  samples = list()
  g = 1
  for (i in seq(1,ncol(df),2)) {
    j = i + 1
    
    sample = df[3:nrow(df), i:j]
    colnames(sample) = df[2, i:j]
    sample[,c(1,2)] = apply(sample[, c(1,2),drop=F], 2,           
                       function(x) as.numeric(as.character(x)))
    
    sample_name = df[1, i]
    
    samples[[g]] = list()
    samples[[g]][[1]] = sample_name
    samples[[g]][[2]] = as.data.frame(sample)
    
    g = g + 1
  }
  return(samples)
}

# gets a list of samples and returns each samples' average peak height and
# time between peaks, and their standard deviations
# samples[[1]][[1]] - sample # (ex. '3560')
# samples[[1]][[2]] - sample data
#                     dataframe of two columns
#                     - 'mean.peak.time (ms)' - 779.1486
#                     - 'mean.peakECG (mV)' - 4207.517
#                     - 'sd.peak.time (ms)' - 779.1486
#                     - 'sd.peakECG (mV)' - 4207.517
getMeans = function(data.list, threshold=0.5) {
  meanlist = list()
  for (i in seq(1, length(data.list))) {

    # `peak.time (ms)` = data.list[[i]][[2]]$`Time (ms)`[ggpmisc:::find_peaks(x = data.list[[i]][[2]]$`Time (ms)`, ignore_threshold=threshold)]
    `peak.time (ms)` = data.list[[i]][[2]]$`Time (ms)`[ggpmisc:::find_peaks(x = data.list[[i]][[2]]$`ECG (mV)`, ignore_threshold=threshold)]
    `peak.ECG (mV)` = data.list[[i]][[2]]$`ECG (mV)`[ggpmisc:::find_peaks(x = data.list[[i]][[2]]$`ECG (mV)`, ignore_threshold=threshold)]
    `peak.ECG (mv) - outliers removed` = `peak.ECG (mV)`[!`peak.ECG (mV)` %in% boxplot.stats(x)$out]

    meanlist[[i]] = list()
    meanlist[[i]][[1]] = data.list[[i]][[1]]
    meanlist[[i]][[2]] = data.frame(mean(`peak.time (ms)`), mean(`peak.ECG (mV)`), sd(`peak.time (ms)`), sd(`peak.ECG (mV)`))
  }
  return(meanlist)
}

# from ggpmisc2
find_peaks <-
  function(x,
           ignore_threshold = 0,
           span = 3,
           strict = TRUE,
           na.rm = FALSE) {
    # find peaks
    if(is.null(span)) {
      pks <- x == max(x, na.rm = na.rm)
      if (strict && sum(pks) != 1L) {
        pks <- logical(length(x)) # all FALSE
      }
    } else {
      pks <- splus2R::peaks(x = x, span = span, strict = strict)
    }
    # apply threshold to found peaks
    if (abs(ignore_threshold) < 1e-5) {
      pks
    } else {
      range_x <- range(x, na.rm = na.rm, finite = TRUE)
      min_x <- range_x[1]
      max_x <- range_x[2]
      x <- ifelse(!is.finite(x), min_x, x)
      # this can cater for the case when max_x < 0, as with logs
      delta <- max_x - min_x
      top_flag <- ignore_threshold > 0.0
      scaled_threshold <- delta * abs(ignore_threshold)
      if (top_flag) {
        ifelse(x - min_x > scaled_threshold, pks , FALSE)
      } else {
        ifelse(max_x - x > scaled_threshold, pks , FALSE)
      }
    }
  }

# gets a list of stable windows from a list of peaks in the table
# return the stable window if conditions met
# return no window if conditions not met
# is any peak; so windowsize of 10 might be too low 
# windowsize is a timespan rather than a # number of peaks
# example stat peaks
#  time | ecg
#  1501 | 13

# TESTING
peaks <- peaks
minPeakStd <- 50
minTimeStd = 50
startWinSize = 200
nHighPeaks = 10

peaks[nrow(peaks), 1]

getStablePeaks = function(peaks, minPeakStd, minTimeStd, startWinSize = 100, nHighPeaks = 10) {

  peaks = na.omit(peaks)
  colnames(peaks) = c("time","ecg")
  
  stablePeaks = list()

  # change window size to include at least 10 peaks to start
  windowSize = startWinSize
  if ((peaks[10, ]$time - peaks[1, ]$time) > startWinSize) {
    windowSize = peaks[10, ]$time - peaks[1, ]$time
  }
  
  winStart = peaks[1, 1]
  winEnd = peaks[1, 1] + windowSize

  max_peaks = data.frame()

  if (!is.numeric(winEnd) || !is.numeric(peaks[nrow(peaks), 1])) {
    return("Error")
  }
  
  while (winEnd < as.numeric(peaks[nrow(peaks), 1])) {
    
    # get the nHighPeaks ECG peaks within the window
    window = na.omit(peaks[peaks$time >= winStart, ])
    colnames(window) = c("time","ecg")
    window = na.omit(window[window$time <= winEnd, ])
    
    max_peaks = window[with(window, order(-ecg)), ] 
    max_peaks = max_peaks[1:nHighPeaks, ]
    colnames(max_peaks) = c("time","ecg")
    max_peaks <- na.omit(max_peaks)
    
    # not enough max peaks
    if (nrow(max_peaks) < 2) {
      # increment window
      winStart = winEnd
      winEnd = winStart + windowSize
      next
    }

    # get the time difference between each
    time_periods = data.frame()
    for (i in 1:(nrow(max_peaks)-1)) {
      time_periods = rbind(time_periods, max_peaks[i + 1, ] - max_peaks[i, ])
    }


    # get the standard deviations
    timeSD = sd(time_periods$time) # get standard deviation of time periods
    ecgSD = sd(max_peaks$ecg) # get standard deviation of ECG


    # # check if it has a small enough stdev for both ECG and time periods
    if ((timeSD <= minTimeStd) && (ecgSD <= minPeakStd)) {
      # no issues, add this as a stable window
      # stablePeaks = rbind(stablePeaks, max_peaks)  
      stablePeaks[[1]] = ecgSD
      stablePeaks[[2]] = window
      
      # increment entire window because all good
      winStart = winEnd
      winEnd = winStart + windowSize
    }
    else {
      stablePeaks[[1]] = ecgSD
      stablePeaks[[2]] = window
      
      # increment window by 2nd element
      winStart = window[2, 1]
      winEnd = winStart + windowSize 
    }

  }

  return(stablePeaks)
}


########################################################

ecg1.peaks = list()
for (i in 1:(length(ecg1.samples))) {

  peaks = ecg1.samples[[i]][[2]][find_peaks(x=ecg1.samples[[i]][[2]]$`ECG (mV)`, ignore_threshold = 0.3, span = 20), ]
  ecg1.peaks[[i]] = list()
  ecg1.peaks[[i]] = getStablePeaks(peaks, minPeakStd=10, minTimeStd=10, startWinSize=50, nHighPeaks=10)

}

ecg2.peaks = list()
for (i in 1:(length(ecg2.samples))) {

  peaks <- ecg2.samples[[i]][[2]][find_peaks(x=ecg1.samples[[1]][[2]]$`ECG (mV)`, ignore_threshold = 0.5, span = 20), ]
  ecg2.peaks <- append(ecg1.peaks, getStablePeaks(peaks))

}


########################################################


ecg1.samples = listSamples(ecg1)
ecg2.samples = listSamples(ecg2)


`peak.time (ms)` = data.list[[i]][[2]]$`Time (ms)`[ggpmisc::find_peaks(x = data.list[[1]][[2]]$`ECG (mV)`, ignore_threshold=threshold)]
`peak.ECG (mV)` = data.list[[i]][[2]]$`ECG (mV)`[ggpmisc::find_peaks(x = data.list[[1]][[2]]$`ECG (mV)`, ignore_threshold=threshold)]

peaks1 = stat_peaks(ignore_threshold=0.5)


with(df,plot(range(ecg1.samples[[1]][[2]]),range(ecg1.samples[[1]][[2]]$`Time (ms)`, ecg1.samples[[1]][[2]]$`ECG (mV)`),type = "n"))
points(ggpmisc:::find_peaks(x = ecg1.samples[[1]][[2]]$`ECG (mV)`, ignore_threshold=0.5)~ecg1.samples[[1]][[2]], data = df, pch=20, col="blue")
points(stat_peaks(col = "red", ignore_threshold=0.5)~x, data = df, pch=20, col="red")

ggplot(data = ecg1.samples[[1]][[2]], aes(x = ecg1.samples[[1]][[2]]$`Time (ms)`, y = ecg1.samples[[1]][[2]]$`ECG (mV)`)) + geom_line() + ggpmisc:::find_peaks(col = "red", threshold=0.5)

# you'll probably have to end up discarding this sample
ggplot(data = ecg1.samples[[1]][[2]], aes(x = ecg1.samples[[4]][[2]]$`Time (ms)`, y = ecg1.samples[[4]][[2]]$`ECG (mV)`)) + geom_line() + stat_peaks(col = "red", span = 20, ignore_threshold=0.5)


lmm = local.min.max(ecg1.samples[[4]][[2]], dev=mean, add.points=TRUE, main="Local Minima and Maxima") 

#### NOTES ####
# 3 peaks for
# window size is 10 peaks

# determine if it's a good window:
# - choose the 10 highest points
#   - if it's stable, those 10 high peaks will be close
#   - if not, then there we a greater amount of variability
# report the # of good windows in the sample

# note: variability is determined via standard deviation
# - smaller standard deviation better
#   - we'll choose the standard deviation cut-off (do one so the last time period
#     aren't included as stable windows)

# expected 10-20 stable windows (less than 10 with cut-off above)

stat_peaks(col = "red", span = 20, ignore_threshold=0.5)

x <- stat_peaks(
  data = ecg1.samples[[1]][[2]],
  span = 20,
  ignore_threshold = 0.5,
)

x = ggplot(data = ecg1.samples[[1]][[2]], aes(x = ecg1.samples[[4]][[2]]$`Time (ms)`, y = ecg1.samples[[4]][[2]]$`ECG (mV)`)) + geom_line() + stat_peaks(col = "red", span = 20, ignore_threshold=0.5)

xf <- x$setup_layer
xm <- x$geom
xx <- x$geom_params

peaks <- ecg1.samples[[1]][[2]][find_peaks(x=ecg1.samples[[1]][[2]]$`ECG (mV)`, ignore_threshold = 0.5, span = 20), ]
ggplot(data = ecg1.samples[[1]][[2]], aes(x = ecg1.samples[[1]][[2]]$`Time (ms)`, y = ecg1.samples[[1]][[2]]$`ECG (mV)`)) + 
  geom_line() +
  geom_point(peaks, mapping=aes(x = peaks$`Time (ms)`, y = peaks$`ECG (mV)`), color="red")
pp <- getStablePeaks(peaks, minPeakStd=0.5, minTimeStd=0.5, startWinSize=50, nHighPeaks=10)

ggplot(data = peaks, aes(x = ecg1.samples[[4]][[2]]$`Time (ms)`, y = ecg1.samples[[4]][[2]]$`ECG (mV)`))

data <- ecg1.samples[[1]][[2]]
x <- data[find_peaks(x = data$`ECG (mV)`, ignore_threshold = 0.5, span = 20, strict = TRUE, na.rm = FALSE), ]

windows <- getStableWindows(x, 10)


df = ecg1.samples[[1]][[2]]
peaks = ecg1.samples[[1]][[2]]$`ECG (mV)`[ggpmisc::find_peaks(x = ecg1.samples[[1]][[2]]$`ECG (mV)`, ignore_threshold=0.5)]
mean(peaks)


for (i in seq(1, length(ecg1))) {
  `Time (ms)` = ecg1.samples[[i]][[2]]$`Time (ms)`
  `ECG (mV)` = ecg1.samples[[i]][[2]]$`ECG (mV)`
  ggplot(data = ecg1.samples[[i]][[2]], aes(x = `Time (ms)`, y = `ECG (mV)`)) + ggtitle(paste("ECG1 Sample ", ecg1.samples[[i]][[1]], "(", i, ") Peaks")) + geom_line() + stat_peaks(col = "red", ignore_threshold=0.5)
  # ggsave(filename=paste0("ecg1.", ecg1.samples[[i]][[1]], "_", i, ".png"), device=png, height=7, width=8)
}

for (i in seq(1, length(ecg2))) {
  `Time (ms)` = ecg2.samples[[i]][[2]]$`Time (ms)`
  `ECG (mV)` = ecg2.samples[[i]][[2]]$`ECG (mV)`
  ggplot(data = ecg2.samples[[i]][[2]], aes(x = `Time (ms)`, y = `ECG (mV)`)) + ggtitle(paste("ECG2 Sample ", ecg2.samples[[i]][[1]], "(", i, ") Peaks")) + geom_line() + stat_peaks(col = "red", ignore_threshold=0.5)
  # ggsave(filename=paste0("ecg2.", ecg2.samples[[i]][[1]], "_", i, ".png"), device=png, height=7, width=8)
}


ecg1.peaks = getMeans(ecg1.samples, 0.5)
ecg2.peaks = getMeans(ecg2.samples, 0.5)

# peaks should have relatively stable level (ecg1.samples[[1]][[2]] not really stable, maybe the middle)
# should be uniformally distributed
# regions don't need to be big
# each sample - find 5 candidate region

# you have multiple recordings from the same mouse on occassion

# easiest way to decide what is "uniform" or "stable" is via the standard deviation
