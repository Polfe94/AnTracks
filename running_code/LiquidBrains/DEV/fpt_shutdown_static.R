path_static_DET <- '/home/polfer/research/gits/AutomatAnts/results/2024/hetero_model/rho/DET/static/'
path_shutdown_DET <- '/home/polfer/research/gits/AutomatAnts/results/2024/hetero_model/rho/DET/shutdown/'
f <- list.files(path_static_DET)
files <- unique(unlist(regmatches(f, gregexpr('rho_\\d{1}\\.\\d{1,3}_\\d{1,2}', f))))

library(latex2exp)
library(arrow)
library(parallel)
# library(infotheo)
source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')

load('~/research/gits/AnTracks/data/nf.RData')
load('~/research/gits/AnTracks/data/det.RData')

nCores <- 20L
## ++++ check distance relationship... roughly experimental distance = 50 * simulation distance ++++ ##
# d1 <- pdist(as.matrix(hex[hex$node == 634, 1:2]), as.matrix(hex[hex$y > 1000, 1:2]))
# d2 <- pdist(as.matrix(hex_sim[node == '(0, 22)', 1:2]), as.matrix(hex_sim[, 1:2]))
# plot(sort(d1) / sort(d2))

R_sims <- seq(2.5, 20, 2.5)
R_exp <- R_sims * 50


# +++ RUN IN CLUSTER !! +++
# t0 <- Sys.time()
# result_shutdown_DET <- mclapply(f[grepl('data', f)], function(x){
# 	do.call('gc', args = list(verbose = FALSE))
# 	data <- data.table(read_parquet(paste0(path_shutdown_DET, x)))[, .(node = parse_nodes(pos)),
# 								       by = 'Frame']
# 	rho <- round(as.numeric(strsplit(x, '_')[[1]][2]), 2)
# 	mdata <- merge(data, hex_sim[, c('x', 'y', 'node')])[order(Frame)]
# 	dists <- pdist(as.matrix(hex_sim[node == '(0, 22)', c('x', 'y')]), as.matrix(mdata[, c('x', 'y')]))
# 	indices <- vapply(R_sims, function(i){
# 		which.max(dists > i)
# 	}, numeric(1))
# 	times <- mdata[indices, Frame]
# 	data.table(R = R_sims, t = times, rho = rho)
# }, mc.cores = nCores)
# Sys.time()-t0
# 
# t0 <- Sys.time()
# result_static_DET <- mclapply(f[grepl('data', f)], function(x){
# 	do.call('gc', args = list(verbose = FALSE))
# 	data <- data.table(read_parquet(paste0(path_static_DET, x)))[, .(node = parse_nodes(pos)),
# 								       by = 'Frame']
# 	rho <- round(as.numeric(strsplit(x, '_')[[1]][2]), 2)
# 	mdata <- merge(data, hex_sim[, c('x', 'y', 'node')])[order(Frame)]
# 	dists <- pdist(as.matrix(hex_sim[node == '(0, 22)', c('x', 'y')]), as.matrix(mdata[, c('x', 'y')]))
# 	indices <- vapply(R_sims, function(i){
# 		which.max(dists > i)
# 	}, numeric(1))
# 	times <- mdata[indices, Frame]
# 	data.table(R = R_sims, t = times, rho = rho)
# }, mc.cores = nCores)
# Sys.time()-t0

# save(result_shutdown_DET, file = '/home/usuaris/pol.fernandez/research/AnTracks/results/first_passage_shutdown_DET.RData')
# save(result_static_DET, file = '/home/usuaris/pol.fernandez/research/AnTracks/results/first_passage_static_DET.RData')

load('~/research/gits/AnTracks/results/fpt_shutdown_DET.RData')
fpt_shutdown_DET <- rbindlist(fpt_shutdown_DET)

ggplot(data = fpt_shutdown_DET[, .(t = mean(t, na.rm = TRUE)), by = c('rho', 'R')],
       aes(factor(rho), factor(R), fill = t)) +
	geom_tile() + scale_fill_viridis(option = 'C')
ggplot(data = fpt_shutdown_DET[is.finite(t), .(t = mean(t)), by = c('rho', 'R')],
       aes(factor(rho), factor(R), fill = t)) +
	geom_tile()
ggplot(data = fpt_shutdown_DET[is.finite(t), .(t = sum(t)), by = c('rho', 'R')],
       aes(factor(rho), factor(R), fill = t)) +
	geom_tile()
ggplot(data = fpt_shutdown_DET,
       aes(factor(rho), factor(R), fill = t)) +
	geom_tile()

ggplot(data = fpt_shutdown_DET[, .(t = median(t, na.rm = TRUE)), by = c('rho', 'R')],
       aes(factor(rho), factor(R), fill = t)) +
	geom_tile()

ggplot(data = fpt_shutdown_DET, aes(factor(R), t)) + geom_boxplot()
ggplot(data = fpt_shutdown_DET, aes(factor(R), t)) + geom_boxplot()+
	facet_wrap(~ rho)

ggplot(data = fpt_shutdown_DET[is.finite(t), .(t = median(t)), by = c('R', 'rho')],
       aes(R, t)) + geom_point()+ geom_path()+
	facet_wrap(~ rho)
ggplot(data = fpt_shutdown_DET[is.finite(t), .(t = mean(t)), by = c('R', 'rho')],
       aes(R, ) + geom_point()+ geom_path()+
	facet_wrap(~ rho)
ggplot(data = fpt_shutdown_DET[rho == 0.45], aes(factor(R), t)) + geom_boxplot()
# ggplot(data = fpt_shutdown_DET[is.finite(t), .(t = mean(t)), by = 'R'], aes(1/(t/2), R*50)) +
# 	geom_point(size = 3)+
# 	geom_smooth(formula = y ~ log(x), se = FALSE) + 
# 	scale_y_continuous('Radius (mm)', breaks = seq(0, 1000, 250))+
# 	scale_x_continuous(TeX('First passage efficiency $(s^{-1})$'))

ggplot(data = fpt_shutdown_DET[is.finite(t), .(t = 1/mean((t/2))), by = 'R'], aes(R*50, t)) +
	geom_point(size = 3)+
	geom_smooth(formula = y ~ log(x), se = FALSE) + 
	scale_x_continuous('Radius (mm)', breaks = seq(0, 1000, 250))+
	scale_y_continuous(TeX('First passage efficiency $(s^{-1})$'))

ggplot(data = fpt_shutdown_DET, aes(factor(R*50), 1/(t/2))) +
	geom_boxplot()+
	scale_x_discrete('Radius (mm)', breaks = seq(0, 1000, 250))+
	scale_y_continuous(TeX('First passage efficiency $(s^{-1})$'))+
	coord_cartesian(ylim = c(0, 0.04))

# ggplot(data = fpt_shutdown_DET[is.finite(t), .(t = 1/mean((t/2)), 
# 					       sd = 1/sd((t/2))), by = 'R'], aes(R*50, t)) +
# 	geom_pointrange(size = 1, aes(ymax = t +sd, ymin = t-sd))+
# 	geom_smooth(formula = y ~ log(x), se = FALSE) + 
# 	scale_x_continuous('Radius (mm)', breaks = seq(0, 1000, 250))+
# 	scale_y_continuous(TeX('First passage efficiency $(s^{-1})$'))

ggplot(data = fpt_shutdown_DET[is.finite(t), .(t = mean(t)), by = c('rho', 'R')],
       aes( R*50, 1/(t/2), color = factor(rho))) +
	geom_point(size = 3)+
	geom_smooth(formula = y ~ log(x), se = FALSE) + 
	scale_x_continuous('Radius (mm)', breaks = seq(0, 1000, 250))+
	scale_y_continuous(TeX('First passage efficiency $(s^{-1})$'))

######## LOOKING EXPLOITATION EFFICIENCY ########
l <- length(f[grepl('food', f)])

static_tps <- rbindlist(lapply(f[grepl('food', f)], function(x){
	food <- data.table(read_parquet(paste0(path_static_DET, x)))[, .(node = revert_node(node),
									   t = t, origin = origin)]
	if(any(is.finite(food[['t']]))){
		food[['patch']] <- c(rep(1, 6), rep(2, 6))
		rho <- round(as.numeric(strsplit(x, '_')[[1]][2]), 2)
		idx <- which.min(food[['t']])
		p1 <- food[idx, patch]
		mint_p1 <- food[idx, t]
		mint_p2 <- food[patch != p1, min(t)]
		tp2_sincr <- mint_p2 - mint_p1
		tp2_min <- min(mint_p2 - food[patch == p1 & t < mint_p2, t], na.rm = TRUE)
		cat(paste0('Finished iteration ', formatC(which(x == f[grepl('food', f)]), flag = '0', digits = 3), ' from ', l, '\r'))
		data.table(tp2_sincr = tp2_sincr, tp2_min = tp2_min, 
			   tp_diff = mean(diff(food[, sort(t)]), na.rm = TRUE),rho = rho)
	}
}))

l <- length(f[grepl('food', f)])

shutdown_tps <- rbindlist(lapply(f[grepl('food', f)], function(x){
	food <- data.table(read_parquet(paste0(path_shutdown_DET, x)))[, .(node = revert_node(node),
									 t = t, origin = origin)]
	if(any(is.finite(food[['t']]))){
		food[['patch']] <- c(rep(1, 6), rep(2, 6))
		rho <- round(as.numeric(strsplit(x, '_')[[1]][2]), 2)
		idx <- which.min(food[['t']])
		p1 <- food[idx, patch]
		mint_p1 <- food[idx, t]
		mint_p2 <- food[patch != p1, min(t)]
		tp2_sincr <- mint_p2 - mint_p1
		tp2_min <- min(mint_p2 - food[patch == p1 & t < mint_p2, t], na.rm = TRUE)
		cat(paste0('Finished iteration ', formatC(which(x == f[grepl('food', f)]), flag = '0', digits = 3), ' from ', l, '\r'))
		data.table(tp2_sincr = tp2_sincr, tp2_min = tp2_min, 
			   tp_diff = mean(diff(food[, sort(t)]), na.rm = TRUE),rho = rho)
	}
}))

ggplot(data = static_tps,
       aes(factor(rho), tp2_sincr)) + geom_boxplot()
ggplot(data = shutdown_tps,
       aes(factor(rho), tp2_sincr)) + geom_boxplot()

ggplot(data = static_tps,
       aes(factor(rho), tp2_min)) + geom_boxplot()
ggplot(data = shutdown_tps,
       aes(factor(rho), tp2_min)) + geom_boxplot()


ggplot(data = static_tps[is.finite(tp2_min), .(tp = mean(tp2_min)), by = 'rho'],
       aes(rho, tp)) + geom_point() +
	geom_smooth(formula = y ~ x, method = 'lm')

ggplot(data = static_tps[is.finite(tp2_min), .(tp = mean(tp2_min)), by = 'rho'],
       aes(rho, tp)) + geom_point() + geom_path()
ggplot(data = static_tps[is.finite(tp2_sincr), .(tp = mean(tp2_sincr)), by = 'rho'],
       aes(rho, tp)) + geom_point() + geom_path()

ggplot(data = shutdown_tps[is.finite(tp2_min), .(tp = mean(tp2_min)), by = 'rho'],
       aes(rho, tp)) + geom_point() +
	geom_smooth(formula = y ~ x, method = 'lm')


ggplot(data = static_tps[is.finite(tp_diff), .(tp = median(tp_diff)), by = 'rho'],
       aes(rho, tp)) + geom_point() + geom_path()
ggplot(data = shutdown_tps[is.finite(tp_diff), .(tp = mean(tp_diff)), by = 'rho'],
       aes(rho, tp)) + geom_point() + geom_path()




det_result <- rbindlist(lapply(det, function(p){
	data <- setDT(p@data)
	mdata <- merge(data[, c('node', 'Frame')], hex[, c('x', 'y', 'node')])[order(Frame)]
	dists <- pdist(as.matrix(hex[hex$node == 634, c('x', 'y')]), as.matrix(mdata[, c('x', 'y')]))
	indices <- vapply(R_exp, function(i){
		which.max(dists > i)
	}, numeric(1))
	times <- mdata[indices, Frame]
	data.table(R = R_exp, t = times, rho = 'det')
}))

nf_result <- rbindlist(lapply(nf, function(p){
	data <- setDT(p@data)
	mdata <- merge(data[, c('node', 'Frame')], hex[, c('x', 'y', 'node')])[order(Frame)]
	dists <- pdist(as.matrix(hex[hex$node == 634, c('x', 'y')]), as.matrix(mdata[, c('x', 'y')]))
	indices <- vapply(R_exp, function(i){
		which.max(dists > i)
	}, numeric(1))
	times <- mdata[indices, Frame]
	data.table(R = R_exp, t = times, rho = 'nf')
}))


ggplot(data = rbindlist(list(det = det_result, nf = nf_result)), 
       aes(factor(R), t / 120, fill = factor(rho)))+
	geom_boxplot(outlier.shape = NA, show.legend = FALSE, alpha = 0.6)+
	geom_jitter(width = 0.15, alpha = 0.4, size = 3, show.legend = FALSE)+
	facet_wrap(~ factor(rho, labels = c('DET', 'NFD')))+
	scale_y_continuous('First passage time (min)', breaks = seq(0, 80, 15))+
	scale_x_discrete('Radius (mm)')+
	scale_fill_manual('', values = c('mediumpurple', 'gold3'))


#################################################
#################################################
#################################################

#### 
LR_movement <- f[grepl('data', f) & grepl('rho_1.001', f)]
SR_movement <- f[grepl('data', f) & grepl('rho_0.001', f)]


maxt <- 2400
maxd <- 10
data <- data.table(read_parquet(paste0(path_static_DET, LR_movement[1])))[Frame <= maxt, .(id = parse_ids(id_out), node = parse_nodes(pos))]
data[, d := shift(node, type = 'lead'), by = id]
filtered <- data[node != d]
mdata <- merge(hex_sim, filtered)
dmatrix <- 
sim_scouts <- rbindlist(lapply(det, function(p){
	data <- setDT(p@data)[Frame <= maxt]
	filter <- data[, .N, by = 'N_ind']
	data <- data[N_ind %in% filter[N >= 10, N_ind]]
	data[['d']] <- get_segment(data[, c('Xmm', 'Ymm')])[, 2]
	colnames(data)[colnames(data) == 'node'] <- 'o'
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	data[['dist']] <- as.numeric(dmatrix)
	LR <- unique(data[dist > maxd, N_ind])
	SR <- unique(data[!N_ind %in% LR, N_ind])
	data_LR <- data[N_ind %in% LR][, scout := 'LR']
	data_SR <- data[N_ind %in% SR][, scout := 'SR']
	rbindlist(list(LR = data_LR, SR = data_SR))
}), idcol = TRUE)
det_scouts[['id']] <- apply(det_scouts[, c('.id', 'N_ind')], 1, paste, collapse = '_', sep = '_')

# nf_scouts <- rbindlist(lapply(nf, function(p){
nf_scouts <- rbindlist(lapply(nf[-c(1, 2)], function(p){ ## filter first two exps
	data <- setDT(p@data)[Frame <= maxt]
	filter <- data[, .N, by = 'N_ind']
	data <- data[N_ind %in% filter[N >= 10, N_ind]]
	data[['d']] <- get_segment(data[, c('Xmm', 'Ymm')])[, 2]
	colnames(data)[colnames(data) == 'node'] <- 'o'
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	data[['dist']] <- as.numeric(dmatrix)
	LR <- unique(data[dist > maxd, N_ind])
	SR <- unique(data[!N_ind %in% LR, N_ind])
	data_LR <- data[N_ind %in% LR][, scout := 'LR']
	data_SR <- data[N_ind %in% SR][, scout := 'SR']
	rbindlist(list(LR = data_LR, SR = data_SR))
}), idcol = TRUE)
nf_scouts[['id']] <- apply(nf_scouts[, c('.id', 'N_ind')], 1, paste, collapse = '_', sep = '_')

det_density <- det_scouts[, .(N = .N), by = c('scout','o','d')][, condition := 'DET'][, z := rank(N)]
nf_density <- nf_scouts[, .(N = .N), by = c('scout','o','d')][, condition := 'NFD'][, z := rank(N)]
traffic_flow <- rbind(det_density, nf_density)
