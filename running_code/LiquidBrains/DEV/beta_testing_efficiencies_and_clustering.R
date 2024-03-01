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
	
	a <- lapply(files, function(i){

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

# path <- '/home/polfer/research/gits/AutomatAnts/results/2024/rec/'
path <- '/home/polfer/research/gits/AutomatAnts/results/2024/default/'
path <- '/home/polfer/research/gits/AutomatAnts/results/2024/beta_0.66/'
path <- '/home/polfer/research/gits/AutomatAnts/results/2024/beta_1'

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

kmod <- k
kmod[c(3,5,15,35)] <- 'low'

data_peak <- rbindlist(data[lapply(k, function(i) i == 'norm') == TRUE], idcol = TRUE)#[, condition := 'sim']
data_peak <- rbindlist(data[lapply(k, function(i) i %in% c('norm', 'flat')) == TRUE], idcol = TRUE)#[, condition := 'sim']
data_peak <- rbindlist(data[lapply(kmod, function(i) i == 'norm') == TRUE], idcol = TRUE)#[, condition := 'sim']

data_avg <- data_peak[, .(N = mean(N, na.rm = TRUE)), by = 'Frame']# [, .id := 0][, condition := 'sim_avg']

ggplot(data = data_peak, aes(Frame, N)) +
	geom_path(aes(group = .id), color = 'grey80', alpha = 0.5)+
	geom_path(data = det_avg, color = 'gold3', linewidth = 2)+
	geom_path(data = data_avg, color = 'mediumpurple', linewidth = 2)

ggplot(data = rbindlist(data, idcol = TRUE), aes(Frame, N)) + geom_path() + facet_wrap(~factor(.id))




files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('wRec_\\d{1,2}', files))))
pls <- lapply(seq_along(files), function(i){
	
	do.call('gc', args = list(verbose = FALSE))
	food <- data.table(read_parquet(paste0(path, files[i] ,'_food.parquet')))
	food <- unlist(food[!is.na(t), .(t = range(t))], use.names = FALSE) *2
	print(paste(food[1], food[2],food[2]-food[1], collapse = ' '))
	pos <- data.table(read_parquet(paste0(path, files[i], '_positions.parquet')))
	
	print(draw_hexagons(edges = sim_edges, color = 'grey50', linewidth = 1.2)+
	      	geom_point(data = pos, aes(x, y, fill = z), shape = 21,size = 4)+
	      	scale_fill_viridis_c(option = 'C') + 
	      	theme(legend.position = 'none', aspect.ratio = 0.5))
	browser()
})






parse_ids <- function(x){
	as.numeric(strsplit(x, ',')[[1]])
}


get_time <- function(x, threshold = 10){
	if(is.unsorted(x)){
		x <- sort(x)
	}
	y <- diff(x)
	y <- y[y <= threshold]
	sum(y)
}


eff <- function(path, filename, minpieces = 0){
	food <- data.table(read_parquet(paste0(path, filename ,'_food.parquet')))
	food <- unlist(food[!is.na(t), .(t = range(t))], use.names = FALSE) *2
	if(length(food)){
		dt <- data.table(read_parquet(paste0(path, filename, '_data.parquet')))
		
		
		filterdt <- dt[Frame <= food[2], .(id = parse_ids(id_out)), by = 'Frame']
		M <- dcast(data = filterdt, Frame ~ id, value.var = 'id')
		.tmp <- rbindlist(lapply(2:ncol(M), function(i){
			idx <- which(!is.na(M[[i]]))
			data.frame(id = unique(M[[i]][idx]),Frame = M[[1]][idx], d=c(1, diff(idx)))
		}))
		.tmp[d > 1, 'd'] <- 0
		.tmp[, interval := cumsum(c(TRUE, diff(d) != 0)), by = id]
		.tmp[, interval := ifelse(interval %% 2 == 0, interval + 1, interval), by = id]
		.tmp_result <- .tmp[, .(start_frame = min(Frame), end_frame = max(Frame)), by = .(id, interval)]
		
		.tp1 <- .tmp_result[start_frame <= food[1]][, end_frame := ifelse(end_frame < food[1], end_frame, food[1])]
		tp1_value <- sum(.tp1[, end_frame] - .tp1[, start_frame])/2
		.tp2 <- .tmp_result[start_frame <= food[2] & end_frame >= food[1]][, start_frame := ifelse(start_frame < food[1], food[1], start_frame)]
		.tp2 <- .tp2[, end_frame := ifelse(end_frame < food[2], end_frame, food[2])]
		tp2_value <- sum(.tp2[, end_frame] - .tp2[, start_frame])/2
		data.table(tp1 = tp1_value, tp2 = tp2_value)
	}
}


t0 <- Sys.time()
result_list <- lapply(seq_along(files), function(x){
	do.call('gc', args = list(verbose = FALSE))
	eff(path, files[x])
})
Sys.time()- t0
names(result_list) <- files
L <- result_list[lapply(result_list,length) > 1]
result_dt <- rbindlist(L, idcol = TRUE)



det_tps <- lapply(det, get_eff)
sto_tps <- lapply(sto, get_eff)
tps <- cbind(condition = 'exp', rbindlist(l = list(rbindlist(det_tps),rbindlist(sto_tps))))

colnames(result_dt)[1] <- 'condition'
result_dt[['condition']] <- 'sim_wRec'




ggplot(data = melt(rbindlist(list(result_dt, tps)), id.vars = 'condition'), 
       aes(condition, 1/value, fill = condition))+
	geom_boxplot(alpha = 0.6)+ 
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')),
		   scales = 'free_y')+
	scale_x_discrete('', labels = c('Experiment', 'With recruitment')) + 
	ylab(TeX('$Efficiency (s^{-1})$'))+
	scale_fill_manual('', values = c('mediumpurple', 'gold3'),
			  labels = c('Experiment', 'With recruitment'))+
	theme(legend.key.size = unit(30, 'pt'))
