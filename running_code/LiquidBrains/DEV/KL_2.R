path_static_DET <- '/home/polfer/research/gits/AutomatAnts/results/2024/hetero_model/rho/DET/static/'
path_shutdown_DET <- '/home/polfer/research/gits/AutomatAnts/results/2024/hetero_model/rho/DET/shutdown/'
f <- list.files(path_static_DET)
files <- unique(unlist(regmatches(f, gregexpr('rho_\\d{1}\\.\\d{1,3}_\\d{1,2}', f))))


static_DET <- rbindlist(lapply(f[grepl('food', f)], function(i){
	rho <- round(as.numeric(strsplit(i, '_')[[1]][2]), 2)
	f <- gsub('_food', '', i)
	d <- data.table(read_parquet(paste0(path_static_DET, f)))
	m1 <- min(d[N >0, Frame])/2
	food <- data.table(read_parquet(paste0(path_static_DET, i)))
	mint <- min(food[['t']]) - m1
	data.table(mint = 1/mint, maxt = 1/(max(food[['t']])-min(food[['t']])), rho = rho)
}), idcol = TRUE)

shutdown_DET <- rbindlist(lapply(f[grepl('food', f)], function(i){
	rho <- round(as.numeric(strsplit(i, '_')[[1]][2]), 2)
	f <- gsub('_food', '', i)
	d <- data.table(read_parquet(paste0(path_shutdown_DET, f)))
	m1 <- min(d[N >0, Frame])/2
	food <- data.table(read_parquet(paste0(path_shutdown_DET, i)))
	mint <- min(food[['t']]) - m1
	data.table(mint = 1/mint, maxt = 1/(max(food[['t']])-min(food[['t']])), rho = rho)
}), idcol = TRUE)


nbins <- 100
trange_x <- c(0, 0.0275) ## mint approximate range
trange_y <- c(0, 0.0075) ## maxt approximate range
sq_x <- seq(trange_x[1], trange_x[2], length.out = nbins)
sq_y <- seq(trange_y[1], trange_y[2], length.out = nbins)
static_DET[, bin_x := cut(mint, breaks = sq_x)]
shutdown_DET[, bin_x := cut(mint, breaks = sq_x)]
static_DET[, bin_y := cut(maxt, breaks = sq_y)]
shutdown_DET[, bin_y := cut(maxt, breaks = sq_y)]

joint_mint <- rbindlist(list(static = static_DET[, c('mint', 'bin_x', 'rho')],
			     shutdown = shutdown_DET[, c('mint', 'bin_x', 'rho')]), idcol = TRUE)
joint_maxt <- rbindlist(list(static = static_DET[, c('maxt', 'bin_y', 'rho')],
			     shutdown = shutdown_DET[, c('maxt', 'bin_y', 'rho')]), idcol = TRUE)


PX_mint_static <- as.numeric(table(static_DET[['bin_x']]) / nrow(static_DET))
PX_mint_static_DETrho0 <- as.numeric(table(static_DET[rho == 0, bin_x]) / nrow(static_DET[rho == 0]))


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
	rbindlist(list(static_DETmint = rbindlist(lapply(unique(static_DET[['rho']]), function(i){
		PX_mint <- as.numeric(table(static_DET[rho == i, bin_x]) / nrow(static_DET[rho == i]))
		data.frame(rho = i, kl = KL(rbind(PX_mint, QX_mint)), var = 'mint', lbl = 'static')
	})),
	static_DETmaxt = rbindlist(lapply(unique(static_DET[['rho']]), function(i){
		PX_maxt <- as.numeric(table(static_DET[rho == i, bin_y]) / nrow(static_DET[rho == i]))
		data.frame(rho = i, kl = KL(rbind(PX_maxt, QX_maxt)), var = 'maxt', lbl = 'static')
	})),
	
	shutdown_DETmint = rbindlist(lapply(unique(shutdown_DET[['rho']]), function(i){
		PX_mint <- as.numeric(table(shutdown_DET[rho == i, bin_x]) / nrow(shutdown_DET[rho == i]))
		data.frame(rho = i, kl = KL(rbind(PX_mint, QX_mint)), var = 'mint', lbl = 'shutdown')
	})),
	shutdown_DETmaxt = rbindlist(lapply(unique(shutdown_DET[['rho']]), function(i){
		PX_maxt <- as.numeric(table(shutdown_DET[rho == i, bin_y]) / nrow(shutdown_DET[rho == i]))
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






################




load('~/research/gits/AnTracks/results/collective_DET_shutdown.RData')
load('~/research/gits/AnTracks/results/collective_DET_static.RData')

coll_shutdown_det_matrix <- rbindlist(coll_shutdown_det_matrix, fill = TRUE)
coll_static_det_matrix <- rbindlist(coll_static_det_matrix, fill = TRUE)


nbins <- 100
trange_x <- c(0, 70000) ## mint approximate range
trange_y <- c(0, 75500) ## maxt approximate range
sq_x <- seq(trange_x[1], trange_x[2], length.out = nbins)
sq_y <- seq(trange_y[1], trange_y[2], length.out = nbins)
coll_static_det_matrix[, bin_x := cut(tp1, breaks = sq_x)]
coll_shutdown_det_matrix[, bin_x := cut(tp1, breaks = sq_x)]
coll_static_det_matrix[, bin_y := cut(tp2, breaks = sq_y)]
coll_shutdown_det_matrix[, bin_y := cut(tp2, breaks = sq_y)]

det_tps <- rbindlist(lapply(det, get_eff))[, bin_x := cut(tp1, breaks = sq_x)][, bin_y := cut(tp2, breaks = sq_y)]
QX_tp1 <- as.numeric(table(det_tps[['bin_x']]) / nrow(det_tps))
QX_tp2 <- as.numeric(table(det_tps[['bin_y']]) / nrow(det_tps))
coll_shutdown_det_matrix


library(philentropy)

KL_collective <- suppressMessages({
	rbindlist(list(coll_static_det_matrixtp1 = rbindlist(lapply(unique(coll_static_det_matrix[['rho']]), function(i){
		PX_tp1 <- as.numeric(table(coll_static_det_matrix[rho == i, bin_x]) / nrow(coll_static_det_matrix[rho == i]))
		data.frame(rho = i, kl = KL(rbind(PX_tp1, QX_tp1)), var = 'tp1', lbl = 'static')
	})),
	coll_static_det_matrixtp2 = rbindlist(lapply(unique(coll_static_det_matrix[['rho']]), function(i){
		PX_tp2 <- as.numeric(table(coll_static_det_matrix[rho == i, bin_y]) / nrow(coll_static_det_matrix[rho == i]))
		data.frame(rho = i, kl = KL(rbind(PX_tp2, QX_tp2)), var = 'tp2', lbl = 'static')
	})),
	
	coll_shutdown_det_matrixtp1 = rbindlist(lapply(unique(coll_shutdown_det_matrix[['rho']]), function(i){
		PX_tp1 <- as.numeric(table(coll_shutdown_det_matrix[rho == i, bin_x]) / nrow(coll_shutdown_det_matrix[rho == i]))
		data.frame(rho = i, kl = KL(rbind(PX_tp1, QX_tp1)), var = 'tp1', lbl = 'shutdown')
	})),
	coll_shutdown_det_matrixtp2 = rbindlist(lapply(unique(coll_shutdown_det_matrix[['rho']]), function(i){
		PX_tp2 <- as.numeric(table(coll_shutdown_det_matrix[rho == i, bin_y]) / nrow(coll_shutdown_det_matrix[rho == i]))
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
