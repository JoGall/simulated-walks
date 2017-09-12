library(fitdistrplus)
library(CircStats)
library(dplyr)
library(zoo)

n <- 30000 #number of frames to simulate

# get real player positions from sample dataset used in publication http://home.ifi.uio.no/paalh/publications/files/mmsys2014-dataset.pdf
d <- read.csv(url("http://home.ifi.uio.no/paalh/dataset/alfheim/2013-11-03/zxy/2013-11-03_tromso_stromsgodset_raw_first.csv"), head=F)
names(d) <- c("ts", "id", "X", "Y", "heading", "direction", "energy", "speed", "total_distance")
# d <- d %>%
#   group_by(id) %>%
#   mutate(speed2 = sqrt((X - lag(X))^2 + (Y - lag(Y))^2))
d$speed <- d$speed / 20


# 1: Fitting distribution to frequency of velocity changes
#------------------------------------------

# get duration of acceleration / deceleration phases (present as peaks and troughs in first derivative of speed)
# e.g. plot speed for one player for one one minute
dd <- subset(d, id == 1)[1:1200,]
peaks <- which(diff(sign(diff(dd$speed)))==-2)+1
troughs <- which(diff(sign(diff(dd$speed)))==2)+1
plot(dd$speed, type = 'l')
points(peaks, dd$speed[peaks], col="red", pch=19)
points(troughs, dd$speed[troughs], col="blue", pch=19)

# get for all data
delta.t <- lapply(unique(d$id), function(x) {
  ss <- subset(d, id == x)
  peaks = which(diff(sign(diff(ss$speed)))==-2)+1
  troughs = which(diff(sign(diff(ss$speed)))==2)+1
  diff(c(peaks, troughs))
}) %>%
  unlist()

# fitting distribution to peak / trough duration data
#hist(delta.t, 20)
#descdist(delta.t[delta.t>0], discrete = F, boot=500)
#fit.gamma <- fitdist(x, "gamma")
#ests <- bootdist(fit.gamma, niter = 100)
#plot(ests)
#quantile(ests, probs=.05)

# sample movement intervals from fitted gamma distribution
go_intervals <- NULL
while(max(cumsum(go_intervals)) < n) {
  next_interval <- round(rgamma(1, fit.gamma$estimate[1], fit.gamma$estimate[2]))
  go_intervals <- c(go_intervals, next_interval)
}

# 2: Fitting distribution to non-zero velocity data
#------------------------------------------

# fit distribution to non-zero speed data
#descdist(d$speed[d$speed>0], discrete = F, boot = 10)

# fit weibull distribution to non-zero velocity data
fit.weibull <- fitdist(d$speed[d$speed>0], "weibull")

# compare fit to real data
#hist(d$speed[d$speed>0], freq = F, 100)
#lines(density(rweibull(n, fit.weibull$estimate[1], fit.weibull$estimate[2])), col="red")

# 3: Simulate periods of zero movement
#------------------------------------------

# get duration of stationary phases
stop_runs <- lapply(unique(d$id), function(x) {
  ss <- subset(d, id == x)
  x <- ifelse(ss$speed > 0, TRUE, FALSE)
  rle(x)$lengths[rle(x)$values==FALSE]
}) %>%
  unlist()
stop_runs <- stop_runs[order(stop_runs)]
stop_runs <- stop_runs[-length(stop_runs)]
#hist(stop_runs, 100)

# calculate percentage of overall time stationary
length(d$speed[d$speed==0]) / length(d$speed) #11%

# fit beta distribution to stationary intervals
stop_runs <- lapply(unique(d$id), function(x) {
  ss <- subset(d, id == x)
  x <- ss$speed
  rle(x)$lengths[rle(x)$values==FALSE]
}) %>%
  unlist()

x <- (stop_runs - min(stop_runs)) / max(stop_runs)
x <- x[x>0]
fit.beta <- fitdistr(x, "beta", list(shape1 = 0.1, shape2 = 0.1),lower=0.001)

# visualise fitted distribution
#hist(x, 100, density=F)
#lines(density(rbeta(1e5, fit.beta$estimate[1], fit.beta$estimate[2])), col="red")

# sample stationary intervals from beta distribution
stop_intervals <- NULL
while(max(cumsum(stop_intervals)) < n * 0.11) {
  next_interval <- round(rbeta(1, fit.beta$estimate[1], fit.beta$estimate[2]) * max(runs) + min(runs))
  stop_intervals <- c(stop_intervals, next_interval)
}

# combine stop / go intervals and randomise
stops <- data.frame(interval = stop_intervals, movement = rep(0, length(stop_intervals)))
gos <- data.frame(interval = go_intervals, movement = rep(1, length(go_intervals)))


sim <- rbind(stops, gos)
sim <- sim[sample(nrow(sim)),]
sim$frame <- c(1, cumsum(sim$interval)[-nrow(sim)])

# sample target velocity from distribution
sim$speed <- rweibull(nrow(sim), fit.weibull$estimate[1], fit.weibull$estimate[2])

# simulate final velocities for n frames
velocity <- data.frame(frame = seq(min(sim$frame), max(sim$frame), by = 1)) %>%
  full_join(sim, by = "frame") %>%
  mutate(movement = zoo::na.locf(movement, na.rm = FALSE)) %>%
  mutate(speed = zoo::na.approx(speed) * movement / 20) %>%
  top_n(n, frame) %>%
  .$speed
  
# 4: Generate x,y-coordinates
#------------------------------------------

# set parameters for walk
rho <- 0.95 #concentration parameter for angles
xlim <- c(0,105)
ylim <- c(0,68)

# generate angles from wrapped normal distribution
ang <- CircStats::rwrpnorm(n-2,0,rho)
ang <- cumsum(c(runif(1,0,2*pi),ang))
x <- rep(NA, n)
y <- rep(NA, n)

# starting point at 0,0
x[1] = xlim[2] / 2
y[1] = ylim[2] / 2

# generate x,y-coordinates within x,y-limits
for (i in 2:n) {
  repeat {
    newx <- x[i-1]+velocity[i-1]*cos(ang[i-1])
    newy <- y[i-1]+velocity[i-1]*sin(ang[i-1])
    
    if (newx > 0 && newx < pitch_length &&
        newy > 0 && newy < pitch_width) break
  }
  x[i] <- newx
  y[i] <- newy
}

# generate final simulated x,y-coords
dat <- data.frame(X = x, Y = y)