#### LIBRARIES, DATA AND GENERIC FUNCTIONS ####
library(dtw)
library(dtwclust)
library(latex2exp)
library(arrow)
library(infotheo)
source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
load("/home/polfer/research/gits/AutomatAnts/data/zSeries.RData")

load('~/research/gits/AnTracks/data/det.RData')

foods_det <- rbindlist(lapply(det, function(i){
	a <- rbindlist(i@food)[['t']]
	data.table(mint = min(a)/2, maxt = max(a)/2)
}))[, lapply(.SD, mean)]

det_avg <- rbindlist(lapply(det, function(i){
	setDT(i@data)
	x <- merge(i@data[, .N, by = 'Frame'], data.table(Frame = 1:21600), all = TRUE, by = 'Frame')
	if(i@date != det[[7]]@date){
		x[is.na(N), 'N'] <- 0		
	}
	x
}))[, .(N = mean(N, na.rm = T)), by = 'Frame']

det_N <- rbindlist(lapply(det, function(i){
	setDT(i@data)
	x <- merge(i@data[, .N, by = 'Frame'], data.table(Frame = 1:21600), all = TRUE, by = 'Frame')
	if(i@date != det[[7]]@date){
		x[is.na(N), 'N'] <- 0		
	}
	x[, .(N = movingAverage(N, t = 60)[1:695], Frame = seq(31, 21545, 31))]
}), idcol = TRUE) 

#### FUNCTIONS ####
revert_node <- function(n){
	n <- do.call('rbind', n)
	apply(n, 1, function(i) paste0('(', paste0(i, collapse = ', '), ')'))
}

parse_nodes <- function(nodes){
	unlist(strsplit(nodes, ';'))
}

parse_ids <- function(ids){
	as.integer(unlist(strsplit(ids, ',')))
}


global_eff <- function(path){
	files <- list.files(path)[grepl('food', list.files(path))]
	results <- lapply(files, function(i){
		food <- data.table(read_parquet(paste0(path, i)))
		1/data.table(mint = min(food[['t']]), maxt = max(food[['t']]))
	})
	rbindlist(results, idcol = TRUE)
}

eff <- function(path, filename, minpieces = 0){
	food <- data.table(read_parquet(paste0(path, filename ,'_food.parquet')))
	food <- unlist(food[!is.na(t), .(t = range(t))], use.names = FALSE) *2
	if(length(food)){
		dt <- data.table(read_parquet(paste0(path, filename, '_data.parquet')))
		
		
		filterdt <- unique(dt[Frame <= food[2], .(id = parse_ids(id_out)), by = 'Frame'])
		M <- dcast.data.table(data = filterdt, Frame ~ id, value.var = 'id', )
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
	names(a) <- files
	a
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


path <- '/home/polfer/research/gits/AutomatAnts/results/2024/gains/'
f <- list.files(path)
f <- f[grepl('food', f)]
files <- unlist(regmatches(f, gregexpr('g_\\d{1}\\.\\d{1,2}_\\d{1,2}', f)))

eff_results <- rbindlist(lapply(seq_along(f), function(i){
	g <- as.numeric(gsub('g_', '', unlist(regmatches(f[i], gregexpr('g_\\d{1}\\.\\d{1,2}', f[i])))))
	data.table(read_parquet(paste0(path, f[i])))[, .(mint = min(t), maxt = max(t), g = g)]
}), idcol = TRUE)
filtered_eff <- eff_results[!is.finite(mint), mint := 10**10]
filtered_eff[!is.finite(maxt), maxt := 10**10]

ggplot(data = eff_results, aes(g, 1/mint, group = g)) + geom_boxplot()
ggplot(data = eff_results, aes(g, 1/maxt, group = g)) + geom_boxplot()

ggplot(data = melt(eff_results, id.vars = c('g', '.id')), aes(g, 1/value)) +
	geom_boxplot(aes(group = g)) + 
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')), scales = 'free')+
	ylab(TeX('Efficiency $(s^{-1})$')) + xlab('Sensitivity')

ggplot(data = melt(eff_results, id.vars = c('g', '.id')), aes(g, 1/value)) +
	geom_boxplot(aes(group = g)) + 
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')), scales = 'free')+
	ylab(TeX('Efficiency $(s^{-1})$')) + xlab('Sensitivity')

# ggplot(data = melt(eff_results[is.finite(mint), .(mint = mean(mint), maxt = mean(maxt), exps = .N), by = 'g'],
# 		   id.vars = c('g', 'exps')),
#        aes(g, 1/(value/exps))) +
# 	geom_point(aes(group = g)) + 
# 	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')), scales = 'free')+
# 	ylab(TeX('Efficiency $(s^{-1})$')) + xlab('Sensitivity')

ggplot(data = melt(rbind(eff_results[is.finite(mint), .(mint = mean(mint), maxt = mean(maxt), exps = .N), by = 'g'],
			 eff_results[!is.finite(mint), .(mint = 100000, maxt = 100000, exps = 1), by = 'g'])[, .(mint = min(mint), maxt = min(maxt), exps = max(exps)), by = 'g'],
		   id.vars = c('g', 'exps')),
       aes(g, 1/(value/exps))) +
	geom_point(aes(group = g)) + 
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')), scales = 'free')+
	ylab(TeX('Efficiency $(s^{-1})$')) + xlab('Sensitivity')


load('/home/polfer/research/gits/AutomatAnts/data/collective_efficiency_gains.RData')
# t0 <- Sys.time()
# coll_eff <- lapply(seq_along(f), function(i){
# 	do.call('gc', args = list(verbose = FALSE))
# 	x <- eff(path, files[i])
# 	x[['g']] <- as.numeric(gsub('g_', '', unlist(regmatches(f[i], gregexpr('g_\\d{1}\\.\\d{1,2}', f[i])))))
# 	x
# })
# Sys.time()-t0

pol <- rbindlist(lapply(coll_eff, function(i){
	if(length(i) == 1){
		data.table(tp1 = 10**10, tp2 = 10**10, g = i[['g']])
	} else {
		i
	}
}), idcol = TRUE)

ggplot(data = melt(pol, id.vars = c('g', '.id')),
       aes(g, 1/(value))) +
	geom_boxplot(aes(group = g)) + 
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')), scales = 'free')+
	ylab(TeX('Efficiency $(s^{-1})$')) +
	scale_x_continuous('Sensitivity', breaks = seq(0, 1, .15))

ggplot(data = melt(pol[, .(tp1 = mean(tp1), tp2 = mean(tp2)), by ='g'], id.vars = c('g')),
       aes(g, 1/(value))) +
	geom_point(aes(group = g)) + 
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')), scales = 'free')+
	ylab(TeX('Efficiency $(s^{-1})$')) + xlab('Sensitivity')
