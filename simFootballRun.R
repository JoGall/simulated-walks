simRun <- function(n, speed_shape = 1.62, speed_scale = 2.51, mov_shape = 16.36, mov_rate = 14.83, stat_shape1 = 0.16, stat_shape2 = 4.02, rho = 0.98, xlim = c(0,105), ylim = c(0,68), stat_pc = 0) {
  
  # get duration of movement phases
  go_intervals <- NULL
  while(max(cumsum(go_intervals)) < n) {
    next_interval <- round(rgamma(1, mov_shape, mov_rate))
    go_intervals <- c(go_intervals, next_interval)
  }
  
  # get duration of stationary phases
  stop_intervals <- NULL
  while(max(cumsum(stop_intervals)) < n * stat_pc) {
    next_interval <- round(rbeta(1, stat_shape1, stat_shape2) * max(runs) + min(runs))
    stop_intervals <- c(stop_intervals, next_interval)
  }

  # combine stop / go intervals and randomise
  stops <- data.frame(interval = stop_intervals, movement = rep(0, length(stop_intervals)))
  gos <- data.frame(interval = go_intervals, movement = rep(1, length(go_intervals)))
  
  # simulate
  sim <- rbind(stops, gos)
  sim <- sim[sample(nrow(sim)),]
  sim$frame <- c(1, cumsum(sim$interval)[-nrow(sim)])
  
  # sample target velocity from distribution
  sim$speed <- rweibull(nrow(sim), speed_shape, speed_scale)
  
  # simulate final velocities for n frames
  velocity <- data.frame(frame = seq(min(sim$frame), max(sim$frame), by = 1)) %>%
    full_join(sim, by = "frame") %>%
    mutate(movement = zoo::na.locf(movement, na.rm = FALSE)) %>%
    mutate(speed = zoo::na.approx(speed) * movement / 20) %>%
    top_n(n, frame) %>%
    .$speed
  
  # generate angles from wrapped normal distribution
  ang <- CircStats::rwrpnorm(n-2,0,rho)
  ang <- cumsum(c(runif(1,0,2*pi),ang))
  x <- rep(NA, n)
  y <- rep(NA, n)
  
  # starting point at centre of pitch
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
  
  data.frame(X = x, Y = y)
}
