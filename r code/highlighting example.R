####################
### HIGHLIGHTING ###
####################

## Example data
set.seed(0)
dat <- data.frame(dates=seq.Date(Sys.Date(), Sys.Date()+99, 1),
                  value=cumsum(rnorm(100)))

## Determine highlighted regions
v <- rep(0, 100)
v[c(5:20, 30:35, 90:100)] <- 1

## Get the start and end points for highlighted regions
inds <- diff(c(0, v))
start <- dat$dates[inds == 1]
end <- dat$dates[inds == -1]
if (length(start) > length(end)) end <- c(end, tail(dat$dates, 1))

## highlight region data
rects <- data.frame(start=start, end=end, group=seq_along(start))

library(ggplot2)
ggplot(data=dat, aes(dates, value)) +
  theme_minimal() +
  geom_line(lty=2, color="steelblue", lwd=1.1) +
  geom_point() +
  geom_rect(data=rects, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=min(dat$value),
                                               ymax=max(dat$value), group=group), color="transparent", fill="orange", alpha=0.3)
