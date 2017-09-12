#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 04/01/2016 by Joe Gallagher 
# joedgallagher@gmail.com
#
# This script contains functions to generate X,Y-trajectories following 
# random walk and correlated walk rules, for the purpose of simulating 
# insect locomotion.

# Functions made by modifying code from the adehabitatLT package (Calenge 
# et al., 2009) and scripts kindly made available by Julien Colomb 
# (https://github.com/jcolomb/CeTrAn/tree/master/CeTrAn/other_codes).
#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

require(adehabitatLT)
require(CircStats)

##FUNCTION DEFINES
##================

##Generate binary movement vector
makeMovementBinary <- function(frames, stat_gamma_shape = 0.221, stat_gamma_rate = 0.00764, mob_gamma_shape = 0.773, mob_gamma_rate = 0.0557, FPS = 19.05) {
#	frames <- time * 60 * FPS
	phases <- frames / 10
	
	#simulate distributions of stationary/mobile times
	#--------------------------------------
	## gamma distribution
	#stationary
# 	stat <- rgamma(phases, 0.217, 0.00792)
	stat <- rgamma(phases, stat_gamma_shape, stat_gamma_rate)
#	stat = stat[stat<30]
	stat = round(stat*FPS)
	#mobile
# 	mob <- rgamma(phases, 0.765, 0.0582)
	mob <- rgamma(phases, mob_gamma_shape, mob_gamma_rate)
	#	mob = mob[mob<30]
	mob = round(mob*FPS)
	
# 	## Levy distribution
# 	#stationary
# 	stat <- (runif(phases)^(1/(1 - stat_mu)) - 1)
# 	stat = round(stat)
# #	stat = stat[stat<300/FPS]
# 	#mobile
# 	mob <- 5 * (runif(phases)^(1/(1 - mob_mu)) - 1)
# 	mob <- round(mob)
# #	mob = mob[mob<30/FPS]

# 	## Cauchy distribution
# 	#stationary
# 	stat <- abs(rcauchy(phases, 4.685, 4.129))
# 	stat = round(stat*FPS)
# 	#mobile
# 	mob <- abs(rcauchy(phases, 0.480, 0.351))
# 	mob <- round(mob*FPS)
	
	#create RLE of mobile/stationary phases
	#--------------------------------------
	#create indexes
	stat_i <- (1:frames)[c(T,F)]
	mob_i <- (1:frames)[c(F,T)]
	#add data
	stat_rle <- data.frame(idx = stat_i, run = stat[1:length(stat_i)], values = rep(0, length(stat_i)))
	mob_rle <- data.frame(idx = mob_i, run = mob[1:length(mob_i)], values = rep(1, length(mob_i)))
	#combine and order
	rle <- rbind(stat_rle, mob_rle)
	rle <- rle[order(rle$idx),]
	rle <- na.omit(rle)
	
	#convert RLE to binary stream (gives a TRUE/FALSE value per frame)
	movement_vector <- rle[rep(1:nrow(rle), times=rle$run),]

	return(movement_vector$values[1:frames])
}

##Correlated walk
simm.crw2 <- function(date=1:100, h = 1, r = 0, R=45, x0=c(0,0), chi_df=0.65, id="A1", burst=id, typeII=TRUE) {
#===========================================================================
     if (!require(CircStats))
         stop("package CircStats required")
     if (typeII)
         class(date) <- c("POSIX","POSIXct")
     n <- length(date)
     dt <- c(diff(unclass(date)))
     if (all(dt-dt[1]>1e-7))
         stop("the time lag between relocations should be constant")
 
	if (h>0) {
		v=sqrt(dt)*rchi(n-1, chi_df) * h
	} else {
		v=dt/dt*(-h)
	}
	
	ang <- rwrpnorm(n-2,0,r)
	ang=cumsum(c(runif(1,0,2*pi),ang))
	x=c(1:n)
	y=x
	x[1]= x0[1]
	y[1] =x0[2]
	
	for (i in c(2: n)){
		x[i]= x[i-1]+v[i-1]*cos(ang[i-1])
		y[i]= y[i-1]+v[i-1]*sin(ang[i-1])
		RR=(x[i]^2+y[i]^2)
		if ((RR> R^2) && (i< (n-1)) ){
			
			x[i]= R*x[i]/sqrt(RR)
			y[i]= R*y[i]/sqrt(RR)
			ang2<-rwrpnorm(n-2-i,0,r)
			ang2=cumsum(c(runif(1,0,2*pi),ang2))
			j=c(1:(n-i-1))
			ang[i+j]=ang2[j]
		}
	}
	
    res <- as.ltraj(data.frame(x,y), date, id, burst, typeII=typeII)
    
    return(res)
}
 
 
##Levy walk
simm.levy2 <- function (date = 1:500, mu = 2, l0 = 1, x0 = c(0, 0), id = "A1", R=45, burst = id, typeII = TRUE, r = 0.6) {
#===========================================================================
    if (typeII) {
		class(date) <- c("POSIX", "POSIXct")
    }
    n <- length(date)
    dt <- c(diff(unclass(date)))
	if (all(dt - dt[1] > 1e-07)) {
    stop("the time lag between relocations should be constant")
    }
    
	ang <- runif(n - 2, -pi, pi)
	v <- dt * (l0 * (runif(n - 1, 1, 30)^(1/(1 - mu))))
	ang = cumsum(c(runif(1, 0, 2 * pi), ang))
    
#    newcode     
	x=c(1:n)
	y=x
	x[1]= x0[1]
	y[1] =x0[2]
	for (i in c(2: n)){
		x[i]= x[i-1]+v[i-1]*cos(ang[i-1])
		y[i]= y[i-1]+v[i-1]*sin(ang[i-1])
		RR=(x[i]^2+y[i]^2)
		
		if ((RR> R^2) && (i< (n-1)) ){
			x[i]= R*x[i]/sqrt(RR)
			y[i]= R*y[i]/sqrt(RR)
			ang2 <- rwrpnorm(n-2-i,0,r)
			ang2 = cumsum(c(runif(1,0,2*pi), ang2))
			j=c(1:(n-i-1))
			ang[i+j]=ang2[j]
			x = c(x0[2], x0[2] + cumsum(v * sin(ang2)))
			y = c(x0[1], x0[1] + cumsum(v * cos(ang2)))
		}
	}
    
    res <- as.ltraj(data.frame(x, y), date, id, burst, typeII = typeII)
   
	return(res)
}
 
 
##Levy correlated walk
simm.levycorr <- function (date = 1:500, mu = 2, l0 = 1, x0 = c(0, 0), id = "A1", R=45, r=0.6, burst = id, typeII = TRUE) {
#===========================================================================
     if (typeII)
         class(date) <- c("POSIX", "POSIXct")
     n <- length(date)
     dt <- c(diff(unclass(date)))
     if (all(dt - dt[1] > 1e-07))
         stop("the time lag between relocations should be constant")
    ang<-rwrpnorm(n-2,0,r)
    ang=cumsum(c(runif(1,0,2*pi),ang))
    v = dt * (l0 * (runif(n - 1, 1, 30)^(1/(1 - mu))))
         
 #    newcode     
     x=c(1:n)
     y=x
     x[1]= x0[1]
     y[1] =x0[2]
     for (i in c(2: n)){
         x[i]= x[i-1]+v[i-1]*cos(ang[i-1])
     		y[i]= y[i-1]+v[i-1]*sin(ang[i-1])
     		RR=(x[i]^2+y[i]^2)
     		if ((RR> R^2) && (i< (n-1)) ){
     			
     			x[i]= R*x[i]/sqrt(RR)
     			y[i]= R*y[i]/sqrt(RR)
     			ang2<-rwrpnorm(n-2-i,0,r)
           ang2=cumsum(c(runif(1,0,2*pi),ang2))
     			j=c(1:(n-i-1))
     			ang[i+j]=ang2[j]
     
     			}
     		}
     res <- as.ltraj(data.frame(x, y), date, id, burst, typeII = typeII)
     return(res)
}


##Correlated walk with pauses
simm.crw2pauses <- function(date=1:100, h=1, r=0, R=45, x0=c(0,0), id="A1", burst=id, typeII=TRUE, act=0.5, chi_df=0.65) {
#===========================================================================
	if (!require(CircStats))
		stop("package CircStats required")
	if (typeII)
		class(date) <- c("POSIX","POSIXct")
	
	n <- length(date)
	dt <- c(diff(unclass(date)))
	
	if (all(dt-dt[1]>1e-7))
		stop("the time lag between relocations should be constant")
    #activity of probability 1-ttt
	t= runif(n - 1,0,1)
	t
	t=ifelse(t<act,1,0)
	t
			
	if (h>0) {
		v=sqrt(dt)*rchi(n-1, chi_df) * h*t
	} else {
		v=dt/dt*(-h)*t
	}
	
	ang <- rwrpnorm(n-2,0,r)
	ang = cumsum(c(runif(1,0,2*pi),ang))
	
	x=c(1:n)
	y=x
	x[1]= x0[1]
	y[1] =x0[2]
	
	for (i in c(2: n)){
		x[i]= x[i-1]+v[i-1]*cos(ang[i-1])
		y[i]= y[i-1]+v[i-1]*sin(ang[i-1])
		RR=(x[i]^2+y[i]^2)
		if ((RR> R^2) && (i< (n-1)) ){	
			x[i]= R*x[i]/sqrt(RR)
			y[i]= R*y[i]/sqrt(RR)
			ang2<-rwrpnorm(n-2-i,0,r)
			ang2=cumsum(c(runif(1,0,2*pi),ang2))
			j=c(1:(n-i-1))
			ang[i+j]=ang2[j]
		}
	}

	res <- as.ltraj(data.frame(x,y),date, id, burst,
					typeII=typeII)
	return(res)
}
   
##Correlated walk with more realistic pauses
simm.crw2pauses2 <- function(date=1:100, h=1, r=0, R=45, x0=c(0,0), id="A1", burst=id, typeII=TRUE, chi_df=0.65, stat_mu = 6, mob_mu = 7.6) {
#===========================================================================
	if (!require(CircStats))
		stop("package CircStats required")
	if (typeII)
		class(date) <- c("POSIX","POSIXct")
	
	n <- length(date)
	dt <- c(diff(unclass(date)))
	
	if (all(dt-dt[1]>1e-7))
		stop("the time lag between relocations should be constant")
    
    #activity of probability bin_mov
	bin_mov <- makeMovementBinary(n-1, stat_mu, mob_mu)
	
	if (h>0) {
		v = sqrt(dt)*rchi(n-1, chi_df) * h * bin_mov
	} else {
		v = dt/dt*(-h) * bin_mov
	}
	
	ang <- rwrpnorm(n-2,0,r)
	ang = cumsum(c(runif(1,0,2*pi),ang))
	
	x=c(1:n)
	y=x
	x[1]= x0[1]
	y[1] =x0[2]
	
	for (i in c(2: n)){
		x[i]= x[i-1]+v[i-1]*cos(ang[i-1])
		y[i]= y[i-1]+v[i-1]*sin(ang[i-1])
		RR=(x[i]^2+y[i]^2)
		if ((RR> R^2) && (i< (n-1)) ){	
			x[i]= R*x[i]/sqrt(RR)
			y[i]= R*y[i]/sqrt(RR)
			ang2<-rwrpnorm(n-2-i,0,r)
			ang2=cumsum(c(runif(1,0,2*pi),ang2))
			j=c(1:(n-i-1))
			ang[i+j]=ang2[j]
		}
	}

	res <- as.ltraj(data.frame(x,y),date, id, burst,
					typeII=typeII)
	return(res)
}
   
##Levy correlated walk with pauses
simm.levycorrpauses <- function (date=1:500, mu=2, l0=1, x0=c(0, 0), id="A1", R=45, r=0.6, act=0.5, burst=id, typeII = TRUE) {
#===========================================================================
	if (typeII)
		class(date) <- c("POSIX", "POSIXct")
		
	n <- length(date)
	dt <- c(diff(unclass(date)))
	
	if (all(dt - dt[1] > 1e-07))
		stop("the time lag between relocations should be constant")
	
	ang <- rwrpnorm(n-2,0,r)
	ang = cumsum(c(runif(1,0,2*pi),ang))
	
	#activity of probability 1-ttt
	t = runif(n-1,0,1)
	t = ifelse(t<act,1,0)
	v = t * dt * (l0 * (runif(n - 1, 1, 30)^(1/(1 - mu))))
          
	#generate x,y
	x = c(1:n)
	y = x
	x[1] = x0[1]
	y[1] = x0[2]
	for (i in c(2: n)){
		x[i] = x[i-1]+v[i-1]*cos(ang[i-1])
		y[i] = y[i-1]+v[i-1]*sin(ang[i-1])
		RR =(x[i]^2+y[i]^2)
		if ((RR> R^2) && (i< (n-1)) ){
			x[i] = R*x[i]/sqrt(RR)
			y[i] = R*y[i]/sqrt(RR)
			ang2 <- rwrpnorm(n-2-i,0,r)
			ang2 = cumsum(c(runif(1,0,2*pi),ang2))
			j = c(1:(n-i-1))
			ang[i+j] = ang2[j]
		}
	}
		
	res <- as.ltraj(data.frame(x, y), date, id, burst, typeII = typeII)
	return(res)
}

##Levy correlated walk with more realistic pauses
simm.levycorrpauses2 <- function (date=1:500, gamma_shape = 0.6891, gamma_rate = 0.1073, Levy_mu = 1.5, l0 = 30, x0=c(0, 0), id="A1", R=45, r=0.99, burst=id, typeII = TRUE, FPS = 19.05) {
#===========================================================================
	if (typeII)
		class(date) <- c("POSIX", "POSIXct")
		
	n <- length(date)
	dt <- c(diff(unclass(date)))
	
	if (all(dt - dt[1] > 1e-07))
		stop("the time lag between relocations should be constant")
	
	#generate angles
	ang <- rwrpnorm(n-2,0,r)
	ang = cumsum(c(runif(1,0,2*pi),ang))
	
	#generate activity with probability bin_mov
	bin_mov <- makeMovementBinary(n-1)
	
	#generate velocity
	#using Levy distribution
# # v <- l0 * (runif(n-1, 1, 7)^(1/(1 - mu)))
# 	v = v - min(v) + 0.1
# #	v = v[v < l0/FPS]
# 	v = v * bin_mov * dt
# 	v = v / FPS
	
	#using gamma distribution
	v <- rgamma(n-1, gamma_shape, gamma_rate)
	v = v - min(v) + 0.1
#	v = v[v < 100/FPS]
	v = v * bin_mov * dt
	v = v / FPS
	
	#generate x,y
	x = c(1:n)
	y = x
	x[1] = x0[1]
	y[1] = x0[2]
	for (i in c(2: n)){
		x[i] = x[i-1]+v[i-1]*cos(ang[i-1])
		y[i] = y[i-1]+v[i-1]*sin(ang[i-1])
		RR =(x[i]^2+y[i]^2)
		if ((RR> R^2) && (i< (n-1)) ){
			x[i] = R*x[i]/sqrt(RR)
			y[i] = R*y[i]/sqrt(RR)
			ang2 <- rwrpnorm(n-2-i,0,r)
			ang2 = cumsum(c(runif(1,0,2*pi),ang2))
			j = c(1:(n-i-1))
			ang[i+j] = ang2[j]
		}
	}
		
	res <- as.ltraj(data.frame(x, y), date, id, burst, typeII = typeII)
	return(res)
}


# corrWalkRect()
#---------------------------------------------------------------------------
# Generate X,Y-coordinates confined to a rectangular area, following a
# correlated walk rule with velocities sampled from a gamma distribution and 
# including pauses
#----------------------------------------------------------------------------
# n = number of samples
# xlim = min and max X-coordinates
# ylim = min and max Y-coordinates
# rho = concentration parameter for generating angles (0 <= rho <= 1)
# h = average velocity
#----------------------------------------------------------------------------
corrWalkRect <- function(n, rho = 0.1, h = 0.2, xlim = c(0,105), ylim = c(0,68)) {
  
  # velocities sampled from chi distribution
  if (h>0) {
    v = rchi(n-1, chi_df) * h
  } else {
    v = -h
  }
  v <- rep(1, n)
  
  # angles sampled from wrapped normal distribution
  ang <- CircStats::rwrpnorm(n-2,0,rho)
  ang <- cumsum(c(runif(1,0,2*pi),ang))
  x <- c(1:n)
  y <- c(1:n)
  x[1] = pitch_length / 2
  y[1] = pitch_width / 2
  
  # generate x,y-coorrdinates within x,y-limits
  for (i in 2:n) {
    repeat {
      newx = x[i-1]+v[i-1]*cos(ang[i-1])
      newy = y[i-1]+v[i-1]*sin(ang[i-1])
      
      ## IF new locations are within bounds, then
      ##    break out of the repeat{} loop (otherwise
      ##    try again)
      if (newx > xlim[1] && newx < xlim[2] &&
          newy > ylim[1] && newy < ylim[2]) break
    }
    x[i] <- newx
    y[i] <- newy
  }
  
  data.frame(X = x, Y = y)
}
