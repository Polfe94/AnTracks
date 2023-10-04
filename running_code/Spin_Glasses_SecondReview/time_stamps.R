load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/sto.RData')

time_stamps <- data.frame()

for(i in seq_along(det)){
	d <- det[[i]]@date
	time_start <- min(det[[i]]@data[['Frame']]) /2 # first ant to arena
	time_interval <- range(do.call('rbind', det[[i]]@food)[['t']]) /2 # time interval of food detection
	
	time_stamps <- rbind(time_stamps, c(d, 'Deterministic',
			  time_start, time_interval, time_interval[2]-time_interval[1]))
}

for(i in seq_along(sto)){
	d <- sto[[i]]@date
	time_start <- min(sto[[i]]@data[['Frame']]) /2 # first ant to arena
	time_interval <- range(do.call('rbind', sto[[i]]@food)[['t']]) /2 # time interval of food stoection
	
	time_stamps <- rbind(time_stamps, c(d, 'Stochastic',
			  time_start, time_interval, time_interval[2]-time_interval[1]))
}

colnames(time_stamps) <- c('Experiment', 'Condition', 'Start_time', 
		  'First_detection', 'Last_detection',
		  'Detection_interval')


write.table(time_stamps,
	    file = "/home/polfer/research/papers/SPIN_GLASSES/review_2/times/time_stamps.csv",
	    sep = ',', dec = '.', row.names = F)
