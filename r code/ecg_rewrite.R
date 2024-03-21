# libraries
library("data.table")
library("ggpmisc")
library("ggplot2")
library("openxlsx")

# files
ecg1 = as.data.frame(data.table::fread("data/ECG1.csv"))
ecg2 = as.data.frame(data.table::fread("data/ECG2.csv"))
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
    sample[c(1,2)] = sapply(sample[c(1, 2)],as.numeric)
    sample_name = df[1, i]
    samples[[g]] = list()
    samples[[g]][[1]] = sample_name
    samples[[g]][[2]] = as.data.frame(sample)
    g = g + 1
  }
  return(samples)
}

getGlobalMaxPeak <- function(peaks) {
  peaks = na.omit(peaks)
  colnames(peaks) = c("time","ecg")
  max_peaks = peaks[with(peaks, order(-ecg)), ]
  return(max_peaks[1, 2])
}

getGlobalMinPeak <- function(peaks) {
  peaks = na.omit(peaks)
  colnames(peaks) = c("time","ecg")
  max_peaks = peaks[with(peaks, order(-ecg)), ]
  return(max_peaks[nrow(peaks), 2])
}

# from ggpmisc2
# returns data.frame of the local maxima from a vector of numerics
find_peaks <- function(x,
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

# get baselines
getBaselines = function(windows) {

  baselines = data.frame()

  # no windows
  if (length(windows) <= 0) {
    return (baselines)
  }

  for (window_pair in windows) {
    window = window_pair[[2]]

    colnames(window) = c("time","ecg")

    # remove outliers and get mean
    z_scores = as.data.frame(sapply(window, function(window) (abs(window-mean(window))/sd(window))))
    no_outliers = z_scores[!rowSums(z_scores>3), ]
    mean = mean(no_outliers$ecg)

    time = window$time
    base = rep(mean, (length(window)))
    baseline = data.frame( time, base )

    baselines = rbind(baselines, baseline)

  }

  return (baselines)
}


# get intersection
getIntersection = function(ecg, baseline) {
  fit <-lm(y~poly(x,2))
  newx <-data.frame(x=seq(0,5,0.01))
  fitline = predict(fit, newdata=newx)

  est <-data.frame(newx,fitline)

  plot(df1,type="o",lwd=2)
  abline(h=1, col="red")
  lines(est, col="blue",lwd=2)

  cross <-est[which.min(abs(1-est$fitline)),] #find closest to 1

  plot(df1,type="o",lwd=2)
  abline(h=1)
  abline(v=cross[1], col="green")
  abline(v=cross[1]-k, col="purple")
  abline(v=cross[1]+k, col="purple")
}


########################################################
### GETTING THE STABLE WINDOWS USING THE ABOVE METHODS ###
ecg1.samples = listSamples(ecg1)
stdevs = c(60, 55, 20, 55, 20, 40, 20, 20, 55, 40)
ecg2.samples = listSamples(ecg2)

ecg1.peaks = list()
ecg1.windows = list()

for (i in 1:(length(ecg1.samples))) {
  # this sample instance
  sample = ecg1.samples[[i]][[2]]
  sample = na.omit(sample)
  sample_name = ecg1.samples[[i]][[1]]

  # remove parts of sample if too high
  tooHigh = sample[sample$`ECG (mV)` > 1000, ]
  if (nrow(tooHigh) > 0) {
    if (tooHigh[1, 1] < sample[nrow(sample) %/% 2, 1]) {
      sample = sample[sample$`Time (ms)` > tooHigh[nrow(tooHigh), 1], ]
    }
    else {
      sample = sample[sample$`Time (ms)` < tooHigh[1, 1], ]
    }
  }

  # remove parts of sample if too low
  tooLow = sample[sample$`ECG (mV)` < -400, ]
  if (nrow(tooLow) > 0) {
    if (tooLow[1, 1] < sample[nrow(sample) %/% 2, 1]) {
      sample = sample[sample$`Time (ms)` > tooLow[nrow(tooLow), 1], ]
    }
    else {
      sample = sample[sample$`Time (ms)` < tooLow[1, 1], ]
    }
  }
  ecg1.samples[[i]][[2]] = sample
  
  # get the local peaks for this sample
  peaks = sample[find_peaks(x=sample$`ECG (mV)`, ignore_threshold = 0.5, span = 20), ]
  
  # try to only get highest peaks
  if (getGlobalMaxPeak(peaks) - getGlobalMinPeak(peaks) > 200) {
    peaks = sample[find_peaks(x=sample$`ECG (mV)`, ignore_threshold = 0.7, span = 20), ]
  }
  ecg1.peaks[[i]] = peaks

  # configure the values for this sample
  noOut = na.omit(sample[!sample$`ECG (mV)` %in% peaks$`ECG (mV)`, ])
  model = lm(`ECG (mV)` ~ `Time (ms)`, data=peaks)
  expectedStdev = stdevs[i]
  peaksPerWindow = 6
  minTimeGap = 60
  windowSize = 1000
  
  # get the windows for this peak
  windows = getStableWindows(peaks = peaks,
  minPeakStd = expectedStdev,
  minTimeGap = minTimeGap,
  windowSize = windowSize,
  nHighPeaks = peaksPerWindow,
  slide = 100)
  ecg1.windows[[i]] = windows

  # configure highlighting of the peaks
  start = c()
  end = c()
  stable.peaks = data.frame()
  if (length(windows) >= 1) { # there are stable peaks in a window
    for (j in 1:(length(windows))) { # iterate through stable peaks
      stable.peaks = windows[[j]][[2]]
      start = append(start, stable.peaks$time[1])
      end = append(end, stable.peaks$time[nrow(stable.peaks)])
    }
  }
  rects = data.frame(start=start, end=end, group=seq_along(start))

  # make a label for the plots
  label0 = paste("window size =", windowSize)
  label1 = paste("min std for ecg =", expectedStdev)
  label2 = paste("time gap =", minTimeGap)
  label3 = paste("# peaks per window =", peaksPerWindow)
  label = paste(label0, label1, label2, label3, sep="\n")

  baseline = getBaselines(windows)
  print(baseline)


  # generate images with sections highlighted
  if (length(windows) > 0) {
    ggplot(data = sample, aes(x = `Time (ms)`, y = `ECG (mV)`)) +
    labs(caption = label) +
    theme(plot.caption = element_text(size=18)) +
    ggtitle(paste("ECG1 Sample ", sample_name, "(", i, ") Peaks")) +
    geom_line() +
    geom_line(baseline, mapping=aes(x = time, y = base), color="steelblue") +
    # geom_line() +
    # geom_point(x = c(2,3), y = rep(dnorm(2, mean = 3), 2), color = "green") +
    geom_line() +
    geom_point(peaks, mapping=aes(x = `Time (ms)`, y = `ECG (mV)`), color="red") +
    geom_rect(data=rects, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=min(peaks$`ECG (mV)`),
    ymax=max(peaks$`ECG (mV)`), group=group), color="transparent",
    fill="orange", alpha=0.3)

  }
  else {
    ggplot(data = sample, aes(x = `Time (ms)`, y = `ECG (mV)`)) +
    labs(caption = label) +
    theme(plot.caption = element_text(size=18)) +
    ggtitle(paste("ECG1 Sample ", sample_name, "(", i, ") Peaks")) +
    geom_line() +
    geom_point(peaks, mapping=aes(x = `Time (ms)`, y = `ECG (mV)`), color="red")

  }
  ggsave(filename=paste0("ecg1rerun3.", sample_name, "_", i, ".png"), device=png, height=10, width=17)

  # export each window to its own sheet page on an excel file
  results = list()
  for (w in 1:length(windows)) {
    result = cbind(windows[[w]][[1]], windows[[w]][[2]]) # bind the stats df and the window df 
    results[[w]] = result
  }
  names(results) = c(paste0("window - ", 1:length(windows)))
  
  write.xlsx(results, file = paste0("ecg1rerun3.", sample_name, "_", i, ".xlsx"), asTable = TRUE)
}