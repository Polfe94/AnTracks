source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')

# load('~/research/gits/AnTracks/data/det.RData')
# 
# rbindlist(lapply(det, function(i){
# 	r <- min(rbindlist(i@food)[['t']])
# 	x <- min(i@data[['Frame']])
# 	data.table(tp1 = r, tp1_norm = r - x)
# }))[, .(t = tp1_norm/ 120)]
# 
# det_tps_med <- rbindlist(lapply(det, function(i){
# 	r <- range(rbindlist(i@food)[['t']]) / 2
# 	data.table(tp1 = r[1], tp2 = (r[2] - r[1]))
# }))[, .(tp1 = 1/median(tp1), tp2 = 1/median(tp2))]

path_SR <- '~/research/gits/AutomatAnts/results/2024/SR_no_listen/'
path_LR <- '~/research/gits/AutomatAnts/results/2024/LR_no_listen/'
path_both <- '~/research/gits/AutomatAnts/results/2024/both_no_listen/'

f <- list.files(path_SR)
files <- f[grepl('food', f)]
l <- length(files)

result_SR <- c()
for(i in seq_along(files)){
	
	test_data <- data.table(read_parquet(paste0(path_SR, files[i])))
	if(nrow(test_data[is.finite(t)])){
		ns_data <- data.table(read_parquet(paste0(path_SR, gsub('_food','', files[i]))))
		f1 <- min(ns_data[N > 0, Frame])
		test_data <- test_data[is.finite(t), .(node = revert_node(node), t = t, origin = revert_node(origin))]
		splt <- strsplit(files[i], '_')[[1]]
		rho <- round(as.numeric(splt[2]), 2)
		eps <- round(as.numeric(splt[4]), 2)
		rng <- range(test_data[['t']]) - f1
		
		result_SR <- rbind(result_SR, data.table(rho = rho, eps = eps, tp1 = rng[1], tp2 = rng[2] - rng[1]))
	}
	cat(paste0('Iter ',i, ' from ', l, '\r'))
}


result_LR <- c()
for(i in seq_along(files)){
	
	test_data <- data.table(read_parquet(paste0(path_LR, files[i])))
	if(nrow(test_data[is.finite(t)])){
		ns_data <- data.table(read_parquet(paste0(path_LR, gsub('_food','', files[i]))))
		f1 <- min(ns_data[N > 0, Frame])
		test_data <- test_data[is.finite(t), .(node = revert_node(node), t = t, origin = revert_node(origin))]
		splt <- strsplit(files[i], '_')[[1]]
		rho <- round(as.numeric(splt[2]), 2)
		eps <- round(as.numeric(splt[4]), 2)
		rng <- range(test_data[['t']]) - f1
		
		result_LR <- rbind(result_LR, data.table(rho = rho, eps = eps, tp1 = rng[1], tp2 = rng[2] - rng[1]))
	}
	cat(paste0('Iter ',i, ' from ', l, '\r'))
}

result_both <- c()
for(i in seq_along(files)){
	
	test_data <- data.table(read_parquet(paste0(path_both, files[i])))
	if(nrow(test_data[is.finite(t)])){
		ns_data <- data.table(read_parquet(paste0(path_both, gsub('_food','', files[i]))))
		f1 <- min(ns_data[N > 0, Frame])
		test_data <- test_data[is.finite(t), .(node = revert_node(node), t = t, origin = revert_node(origin))]
		splt <- strsplit(files[i], '_')[[1]]
		rho <- round(as.numeric(splt[2]), 2)
		eps <- round(as.numeric(splt[4]), 2)
		rng <- range(test_data[['t']]) - f1
		
		result_both <- rbind(result_both, data.table(rho = rho, eps = eps, tp1 = rng[1], tp2 = rng[2] - rng[1]))
	}
	cat(paste0('Iter ',i, ' from ', l, '\r'))
}

