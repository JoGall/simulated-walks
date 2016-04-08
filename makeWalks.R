#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 04/01/2016 by Joe Gallagher 
# joedgallagher@gmail.com
#
# This script contains a function to produce a set of simulated walks 
# between a range of speed and turning parameters, for use with the
# functions defined in 'walkFuns.R'.

#–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

#source simulated walk script
source('./walkFuns.R')

#default parameters
# FPS <- 19.05
# CORR_SPEED <- 10 / FPS
# LEVY_SPEED <- 50 / FPS

##==================================================================
## FUNCTION DEFINE
##==================================================================

makeWalks <- function(reps, minutes=10, FPS=19.05) {
	#create storage
	full_df <- data.frame()
	
	#vary the parameters
	i_speed_gam_shape <- rnorm(reps, 0.695, 0.144)
	i_speed_gam_rate <- rnorm(reps, 0.119, 0.0335)

	i_r <- runif(reps, 0.950, 0.999)
	i_gam_s <- runif(reps, 0.61, 0.65)
	i_gam_r <- runif(reps, 0.05, 0.15)

	#make progress bar
	tmp <- 0
	pb <- txtProgressBar(min = 0, max = reps, style = 3)
	
	#loop	
	for(i in 1:reps){			
		#simulate walk
		frames <- minutes * 60 * FPS
		##continuous correlated walk
# 		a <- simm.crw2(1:frames, h = i_h[i], r = i_r[i], chi_df = 0.65)
		##continuous Levy walk
# 		a <- simm.levy2(1:frames, mu = 2.2, l0 = i_l0[i], r = i_r[i])
		##Levy correlated walk
#		simm.levycorr(1:frames, mu = 2.2, l0 = i_l0[i], r = i_r[i])
		##correlated walk with random pauses
#  		a <- simm.crw2pauses(1:frames, h = i_h[i], r = i_r[i], act=0.3, chi_df=0.65)
		##correlated walk with RLE pauses
# 		a <- simm.crw2pauses2(1:frames, h = i_h[i], r = i_r[i], chi_df=0.65, stat_mu = i_stat_mu[i], mob_mu = i_mob_mu[i])
		##Levy walk with random pauses
# 		a <- simm.levycorrpauses(1:frames, mu = 2.2, l0 = i_l0[i], r = i_r[i], act=0.3)
		##Levy walk with RLE pauses
# 		a <- simm.levycorrpauses2(1:frames, gamma_shape = i_gam_s[i], gamma_rate = i_gam_r[i], l0 = i_l0[i], r = i_r[i])
		
		a <- simm.levycorrpauses2(1:frames, r = 0.994)

		#convert to dataframe and add time variable
		a = as.data.frame(a)
		tt <- seq(0, minutes*60, 1/FPS) * 1000
		tt = tt[-length(tt)]
		
		#change format to work with RUbitrail
		ss_df <- data.frame(id = rep(i, nrow(a)), exp = rep(NA, nrow(a)), rep = rep(NA, nrow(a)), area = rep(NA, nrow(a)), X = a$x, Y = a$y, Distance = a$dist, time = tt)
		
		full_df <- rbind(full_df, ss_df)
		
		# update progress bar
		Sys.sleep(0.1)
		tmp <- tmp+1
		setTxtProgressBar(pb, tmp)
		}
	close(pb)

	return(full_df)
}

#==================================================================

# #Example usage

# write.csv(makeWalks(20, time=10), "CW_sims.csv")

#==================================================================