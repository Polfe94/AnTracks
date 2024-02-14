path_static <- '/home/polfer/research/gits/AutomatAnts/results/2024/hetero_model/rho/static/'
path_shutdown <- '/home/polfer/research/gits/AutomatAnts/results/2024/hetero_model/rho/shutdown/'
f <- list.files(path_static)
# rhos <- unique(unlist(regmatches(f, gregexpr('rho_\\d{1}\\.\\d{1,3}', f))))
# rhos <- round(as.numeric(gsub('rho_', '', rhos)), 2)
files <- unique(unlist(regmatches(f, gregexpr('rho_\\d{1}\\.\\d{1,3}_\\d{1,2}', f))))


#### LIBRARIES, DATA AND GENERIC FUNCTIONS ####
library(dtw)
library(dtwclust)
library(latex2exp)
library(arrow)
# library(infotheo)
source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')

load('~/research/gits/AnTracks/data/det.RData')


## EXPERIMENTS
# foods_det <- rbindlist(lapply(det, function(i){
# 	a <- rbindlist(i@food)[['t']]
# 	data.table(mint = min(a)/2, maxt = max(a)/2)
# }))[, lapply(.SD, mean)] # <--- mean !
foods_det <- rbindlist(lapply(det, function(i){
	a <- rbindlist(i@food)[['t']]
	data.table(mint = min(a)/2, maxt = max(a)/2)
}))[, lapply(.SD, median)] # <--- median !

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

static_ <- rbindlist(lapply(f[grepl('food', f)], function(i){
	rho <- round(as.numeric(strsplit(i, '_')[[1]][2]), 2)
	f <- gsub('_food', '', i)
	d <- data.table(read_parquet(paste0(path_static, f)))
	m1 <- min(d[N >0, Frame])/2
	food <- data.table(read_parquet(paste0(path_static, i)))
	mint <- min(food[['t']]) - m1
	data.table(mint = 1/mint, maxt = 1/(max(food[['t']])-min(food[['t']])), rho = rho)
}), idcol = TRUE)

shutdown_ <- rbindlist(lapply(f[grepl('food', f)], function(i){
	rho <- round(as.numeric(strsplit(i, '_')[[1]][2]), 2)
	f <- gsub('_food', '', i)
	d <- data.table(read_parquet(paste0(path_shutdown, f)))
	m1 <- min(d[N >0, Frame])/2
	food <- data.table(read_parquet(paste0(path_shutdown, i)))
	mint <- min(food[['t']]) - m1
	data.table(mint = 1/mint, maxt = 1/(max(food[['t']])-min(food[['t']])), rho = rho)
}), idcol = TRUE)


nbins <- 100
trange_x <- c(0, 0.0275) ## mint approximate range
trange_y <- c(0, 0.0075) ## maxt approximate range
sq_x <- seq(trange_x[1], trange_x[2], length.out = nbins)
sq_y <- seq(trange_y[1], trange_y[2], length.out = nbins)
static_[, bin_x := cut(mint, breaks = sq_x)]
shutdown_[, bin_x := cut(mint, breaks = sq_x)]
static_[, bin_y := cut(maxt, breaks = sq_y)]
shutdown_[, bin_y := cut(maxt, breaks = sq_y)]

joint_mint <- rbindlist(list(static = static_[, c('mint', 'bin_x', 'rho')],
	       shutdown = shutdown_[, c('mint', 'bin_x', 'rho')]), idcol = TRUE)
joint_maxt <- rbindlist(list(static = static_[, c('maxt', 'bin_y', 'rho')],
			     shutdown = shutdown_[, c('maxt', 'bin_y', 'rho')]), idcol = TRUE)


PX_mint_static <- as.numeric(table(static_[['bin_x']]) / nrow(static_))
PX_mint_static_rho0 <- as.numeric(table(static_[rho == 0, bin_x]) / nrow(static_[rho == 0]))


#### FORAGING EFFICIENCY ####
det_global_eff <- rbindlist(lapply(det, function(i){
	x <- rbindlist(i@food)[['t']]
	mint <- min(x) - min(i@data[['Frame']])
	1/(data.table(mint = mint, maxt = max(x)-min(x))/2)
}))[, condition := 'Experiments']
det_global_eff[, bin_x := cut(mint, breaks = sq_x)]
det_global_eff[, bin_y := cut(maxt, breaks = sq_y)]

QX_mint <- as.numeric(table(det_global_eff[['bin_x']]) / nrow(det_global_eff))
QX_maxt <- as.numeric(table(det_global_eff[['bin_y']]) / nrow(det_global_eff))


library(philentropy)

KL_df <- suppressMessages({
	rbindlist(list(static_mint = rbindlist(lapply(unique(static_[['rho']]), function(i){
	PX_mint <- as.numeric(table(static_[rho == i, bin_x]) / nrow(static_[rho == i]))
	data.frame(rho = i, kl = KL(rbind(PX_mint, QX_mint)), var = 'mint', lbl = 'static')
})),
static_maxt = rbindlist(lapply(unique(static_[['rho']]), function(i){
	PX_maxt <- as.numeric(table(static_[rho == i, bin_y]) / nrow(static_[rho == i]))
	data.frame(rho = i, kl = KL(rbind(PX_maxt, QX_maxt)), var = 'maxt', lbl = 'static')
})),

shutdown_mint = rbindlist(lapply(unique(shutdown_[['rho']]), function(i){
	PX_mint <- as.numeric(table(shutdown_[rho == i, bin_x]) / nrow(shutdown_[rho == i]))
	data.frame(rho = i, kl = KL(rbind(PX_mint, QX_mint)), var = 'mint', lbl = 'shutdown')
})),
shutdown_maxt = rbindlist(lapply(unique(shutdown_[['rho']]), function(i){
	PX_maxt <- as.numeric(table(shutdown_[rho == i, bin_y]) / nrow(shutdown_[rho == i]))
	data.frame(rho = i, kl = KL(rbind(PX_maxt, QX_maxt)), var = 'maxt', lbl = 'shutdown')
}))))})

ggplot(data = KL_df, aes(rho, kl, color = lbl)) + geom_point(size = 4, shape = 21, fill =NA)+
	geom_path() + 
	scale_y_continuous('Kullback-Leibler Divergence', breaks = seq(0, 12, 1.25))+
	xlab(TeX('Proportion of LR scouts ($\\rho$)'))+
	scale_color_manual('', labels = c('Shutdown', 'Static'),
			   values = c('mediumpurple','gold3')) +
	facet_wrap(~ factor(var, levels = c('mint', 'maxt'),
			    labels = c('Exploration time', 'Exploitation time')))






# collective efficiency
l <- length(files)
t0 <- Sys.time()
collective_static <- lapply(seq_along(files), function(i){
	cat(paste0('Iteration = ', formatC(i, width = 4, format = 'd', flag = '0'), ' from ', l, '\r'))
	do.call('gc', args = list(verbose = FALSE))
	df <- eff(path_static, files[i])
	df[['rho']] <- round(as.numeric(strsplit(files[i], '_')[[1]][2]), 2)
})
Sys.time()-t0

print('--- finished static ---')
print('')

t0 <- Sys.time()
collective_shutdown <- lapply(seq_along(files), function(i){
	cat(paste0('Iteration = ', formatC(i, width = 4, format = 'd', flag = '0'), ' from ', l, '\r'))
	do.call('gc', args = list(verbose = FALSE))
	df <- eff(path_shutdown, files[i])
	df[['rho']] <- round(as.numeric(strsplit(files[i], '_')[[1]][2]), 2)
})
Sys.time()-t0

print('--- finished shutdown ---')
print('')



coll_dt <- rbindlist(collective_sim_effs)
coll_dt[['condition']] <- 'Simulations'


load('~/research/gits/AnTracks/results/coll_shutdown.RData')
load('~/research/gits/AnTracks/results/coll_static.RData')



nbins <- 100
trange_x <- c(0, 70000) ## mint approximate range
trange_y <- c(0, 75500) ## maxt approximate range
sq_x <- seq(trange_x[1], trange_x[2], length.out = nbins)
sq_y <- seq(trange_y[1], trange_y[2], length.out = nbins)
coll_static[, bin_x := cut(tp1, breaks = sq_x)]
coll_shutdown[, bin_x := cut(tp1, breaks = sq_x)]
coll_static[, bin_y := cut(tp2, breaks = sq_y)]
coll_shutdown[, bin_y := cut(tp2, breaks = sq_y)]

det_tps <- rbindlist(lapply(det, get_eff))[, bin_x := cut(tp1, breaks = sq_x)][, bin_y := cut(tp2, breaks = sq_y)]
QX_tp1 <- as.numeric(table(det_tps[['bin_x']]) / nrow(det_tps))
QX_tp2 <- as.numeric(table(det_tps[['bin_y']]) / nrow(det_tps))
coll_shutdown


library(philentropy)

KL_collective <- suppressMessages({
	rbindlist(list(coll_statictp1 = rbindlist(lapply(unique(coll_static[['rho']]), function(i){
		PX_tp1 <- as.numeric(table(coll_static[rho == i, bin_x]) / nrow(coll_static[rho == i]))
		data.frame(rho = i, kl = KL(rbind(PX_tp1, QX_tp1)), var = 'tp1', lbl = 'static')
	})),
	coll_statictp2 = rbindlist(lapply(unique(coll_static[['rho']]), function(i){
		PX_tp2 <- as.numeric(table(coll_static[rho == i, bin_y]) / nrow(coll_static[rho == i]))
		data.frame(rho = i, kl = KL(rbind(PX_tp2, QX_tp2)), var = 'tp2', lbl = 'static')
	})),
	
	coll_shutdowntp1 = rbindlist(lapply(unique(coll_shutdown[['rho']]), function(i){
		PX_tp1 <- as.numeric(table(coll_shutdown[rho == i, bin_x]) / nrow(coll_shutdown[rho == i]))
		data.frame(rho = i, kl = KL(rbind(PX_tp1, QX_tp1)), var = 'tp1', lbl = 'shutdown')
	})),
	coll_shutdowntp2 = rbindlist(lapply(unique(coll_shutdown[['rho']]), function(i){
		PX_tp2 <- as.numeric(table(coll_shutdown[rho == i, bin_y]) / nrow(coll_shutdown[rho == i]))
		data.frame(rho = i, kl = KL(rbind(PX_tp2, QX_tp2)), var = 'tp2', lbl = 'shutdown')
	}))))})

ggplot(data = KL_collective, aes(rho, kl, color = lbl)) + geom_point(size = 4, shape = 21, fill =NA)+
	geom_path() + 
	scale_y_continuous('Kullback-Leibler Divergence', breaks = seq(0, 12, 1.25))+
	xlab(TeX('Proportion of LR scouts ($\\rho$)'))+
	scale_color_manual('', labels = c('Shutdown', 'Static'),
			   values = c('mediumpurple','gold3')) +
	facet_wrap(~ factor(var, levels = c('tp1', 'tp2'),
			    labels = c('Exploration efficiency', 'Exploitation efficiency')))
