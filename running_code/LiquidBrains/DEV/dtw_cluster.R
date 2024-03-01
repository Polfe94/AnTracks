#### LIBRARIES, DATA AND GENERIC FUNCTIONS ####
library(dtw)
library(dtwclust)
library(latex2exp)
library(arrow)
source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')

load("/home/polfer/research/gits/AutomatAnts/data/zSeries.RData")
load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/sto.RData')


ndtw <- function(x, y, ...) {
	dtw::dtw(x, y, step.pattern = asymmetric,
		 distance.only = TRUE, ...)$normalizedDistance
}

process_data <- function(path, 
			 vars = c('Frame', 'N')){
	ref <- data.table(Frame = 1:21600)
	t <- movingAverage(ref[['Frame']], 60, 0)

	files <- list.files(path)
	files <- files[!grepl('food', files) & !grepl('position', files) &
		       	!grepl('data', files) & !grepl('keys', files)]

	lapply(files, function(i){
		y <- merge(ref, data.table(read_parquet(paste0(path, i)))[, ..vars], all = TRUE, by = 'Frame')
		y[['N']] <- fillVec(y[['N']])
		n <- movingAverage(y[['N']], 60L, 0L)
		z <- zscore(n)
		data.table(Frame = t, N = n, Z = z)
	})
}

plot_clusters <- function(data_list, clusters){
	kprop <- 100* (table(clusters) / length(clusters))
	n <- names(kprop)
	xPlot <- rbindlist(lapply(seq_along(data_list), function(i){
		data.table(Frame = data_list[[i]][['Frame']], N = data_list[[i]][['N']], k = clusters[i], exp = i)
	}))
	
	ggplot(data = xPlot, aes(Frame, N)) + geom_path(aes(group = exp), color = 'grey70', alpha = 0.75)+
		geom_path(data = xPlot[, .(N = mean(N), exp = exp), by = c('Frame', 'k')], linewidth = 1)+
		facet_wrap(~factor(k, levels = c('norm', 'flat', 'late', 'low'),
				   labels = c(paste0('Experimental pattern (', 
				   		  ifelse('norm' %in% n, kprop[['norm']], 0),
				   		  '%)'), 
				   	   paste0('Plateau (', 
				   	          ifelse('flat' %in% n, kprop[['flat']], 0),
				   	          '%)'), 
				   	   paste0('Late start (', 
				   	          ifelse('late' %in% n, kprop[['late']], 0),
				   	          '%)'), 
				   	   paste0('Low activity (', 
				   	          ifelse('low' %in% n, kprop[['low']], 0),
				   	          '%)'))), 
			   nrow = 1) +
		ylab('Activity (number of ants in arena)') + 
		scale_x_continuous('Time (min)', breaks = seq(0, 21600, 30 * 120), labels = seq(0, 180, 30))+
		theme(strip.text = element_text(size = 16, margin = margin(t = 5, b = 5, unit = 'pt')),
		      aspect.ratio = 0.75)
}

path <- '/home/polfer/research/gits/AutomatAnts/results/rec_noRec_nf/rec/uniform/'

knames <- c('norm', 'late', 'flat', 'low')
data <- process_data(path)
distances <- rbindlist(lapply(data, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
}))
colnames(distances) <- knames
k <- knames[apply(distances, 1, which.min)]

plot_clusters(data, k)

det_avg <- rbindlist(lapply(det, function(i){
	data.table(i@data)[, .N, by = 'Frame']
}))[order(Frame)][, .(N = mean(N)), by = 'Frame'][N == 1 & Frame < 500, N := 0]

det_avg <- rbindlist(lapply(det, function(i){
	setDT(i@data)
	x <- merge(i@data[, .N, by = 'Frame'], data.table(Frame = 1:21600), all = TRUE, by = 'Frame')
	if(i@date != det[[7]]@date){
		x[is.na(N), 'N'] <- 0		
	}
	x
}))[, .(N = mean(N, na.rm = T)), by = 'Frame']
# det_avg[['condition']] <- 'det'
# det_avg[['.id']] <- -1

data_peak <- rbindlist(data[lapply(k, function(i) i == 'norm') == TRUE], idcol = TRUE)#[, condition := 'sim']
data_avg <- data_peak[, .(N = mean(N, na.rm = TRUE)), by = 'Frame']# [, .id := 0][, condition := 'sim_avg']

ggplot(data = data_peak, aes(Frame, N)) +
	geom_path(aes(group = .id), color = 'grey80', alpha = 0.6)+
	geom_path(data = det_avg, color = 'gold3', linewidth = 2)+
	geom_path(data = data_avg, color = 'mediumpurple', linewidth = 2)


############################################################################
############################################################################
############################################################################
############################################################################
############################################################################


path <- '/home/polfer/research/gits/AutomatAnts/results/rec_noRec_nf/noRec/uniform/'

knames <- c('norm', 'late', 'flat', 'low')
data <- process_data(path)
distances <- rbindlist(lapply(data, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
}))
colnames(distances) <- knames
k <- knames[apply(distances, 1, which.min)]

plot_clusters(data, k)

det_avg <- rbindlist(lapply(det, function(i){
	data.table(i@data)[, .N, by = 'Frame']
}))[order(Frame)][, .(N = mean(N)), by = 'Frame'][N == 1 & Frame < 500, N := 0]

det_avg <- rbindlist(lapply(det, function(i){
	setDT(i@data)
	x <- merge(i@data[, .N, by = 'Frame'], data.table(Frame = 1:21600), all = TRUE, by = 'Frame')
	if(i@date != det[[7]]@date){
		x[is.na(N), 'N'] <- 0		
	}
	x
}))[, .(N = mean(N, na.rm = T)), by = 'Frame']
# det_avg[['condition']] <- 'det'
# det_avg[['.id']] <- -1

data_peak <- rbindlist(data[lapply(k, function(i) i == 'norm') == TRUE], idcol = TRUE)#[, condition := 'sim']
data_avg <- data_peak[, .(N = mean(N, na.rm = TRUE)), by = 'Frame']# [, .id := 0][, condition := 'sim_avg']

ggplot(data = data_peak, aes(Frame, N)) +
	geom_path(aes(group = .id), color = 'grey80', alpha = 0.6)+
	geom_path(data = det_avg, color = 'gold3', linewidth = 2)+
	geom_path(data = data_avg, color = 'mediumpurple', linewidth = 2)



############################################################################
############################################################################
############################################################################
############################################################################
############################################################################


path <- '/home/polfer/research/gits/AutomatAnts/results/2024/uniform_noRec/'

knames <- c('norm', 'late', 'flat', 'low')
data <- process_data(path)
distances <- rbindlist(lapply(data, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
}))
colnames(distances) <- knames
k <- knames[apply(distances, 1, which.min)]

plot_clusters(data, k)

det_avg <- rbindlist(lapply(det, function(i){
	data.table(i@data)[, .N, by = 'Frame']
}))[order(Frame)][, .(N = mean(N)), by = 'Frame'][N == 1 & Frame < 500, N := 0]

det_avg <- rbindlist(lapply(det, function(i){
	setDT(i@data)
	x <- merge(i@data[, .N, by = 'Frame'], data.table(Frame = 1:21600), all = TRUE, by = 'Frame')
	if(i@date != det[[7]]@date){
		x[is.na(N), 'N'] <- 0		
	}
	x
}))[, .(N = mean(N, na.rm = T)), by = 'Frame']
# det_avg[['condition']] <- 'det'
# det_avg[['.id']] <- -1

data_peak <- rbindlist(data[lapply(k, function(i) i == 'norm') == TRUE], idcol = TRUE)#[, condition := 'sim']
data_avg <- data_peak[, .(N = mean(N, na.rm = TRUE)), by = 'Frame']# [, .id := 0][, condition := 'sim_avg']

ggplot(data = data_peak, aes(Frame, N)) +
	geom_path(aes(group = .id), color = 'grey80', alpha = 0.6)+
	geom_path(data = det_avg, color = 'gold3', linewidth = 2)+
	geom_path(data = data_avg, color = 'mediumpurple', linewidth = 2)

