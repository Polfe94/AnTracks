source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')

load('~/research/gits/AnTracks/data/det.RData')

det_tps <- rbindlist(lapply(det, function(i){
	r <- range(rbindlist(i@food)[['t']]) / 2
	data.table(tp1 = r[1], tp2 = (r[2] - r[1]))
}))[, .(tp1 = 1/mean(tp1), tp2 = 1/mean(tp2))]

det_tps_med <- rbindlist(lapply(det, function(i){
	r <- range(rbindlist(i@food)[['t']]) / 2
	data.table(tp1 = r[1], tp2 = (r[2] - r[1]))
}))[, .(tp1 = 1/median(tp1), tp2 = 1/median(tp2))]

path_SR <- '~/research/gits/AutomatAnts/results/2024/SR_no_listen/'
path_LR <- '~/research/gits/AutomatAnts/results/2024/LR_no_listen/'
path_both <- '~/research/gits/AutomatAnts/results/2024/both_no_listen/'

patch_index <- c(rep(1, 6), rep(2, 6))
patch_1 <- c('(6, 33)', '(6, 34)', '(7, 34)', '(7, 33)', '(7, 32)', '(6, 32)')
patch_2 <- c('(6, 11)', '(6, 12)', '(7, 12)', '(7, 11)', '(7, 10)', '(6, 10)')

fSR <- list.files(path_SR)
fLR <- list.files(path_LR)
fb <- list.files(path_both)

files_SR <- fSR[grepl('food', fSR)]
files_LR <- fLR[grepl('food', fLR)]
files_both <- fb[grepl('food', fb)]

result_SR <- c()
l <- length(files_SR)
t0 <- Sys.time()
for(i in seq_along(files_SR)){
	
	test_data <- data.table(read_parquet(paste0(path_SR, files_SR[i])))
	test_data[['patch']] <- patch_index
	if(nrow(test_data[is.finite(t)])){
		ns_data <- data.table(read_parquet(paste0(path_SR, gsub('_food','', files_SR[i]))))
		test_data <- test_data[is.finite(t), .(node = revert_node(node), t = t, origin = revert_node(origin))]
		f1 <- min(ns_data[N > 0, Frame])
		splt <- strsplit(files_SR[i], '_')[[1]]
		rho <- round(as.numeric(splt[2]), 2)
		eps <- round(as.numeric(splt[4]), 2)
		rng <- range(test_data[['t']]) - f1
		inest <- test_data[['origin']] == '(0, 22)'
		ip1 <- test_data[['origin']] %in% patch_1
		ip2 <- test_data[['origin']] %in% patch_2
		
		pinterpatch <- sum(cbind(test_data[['patch']] == 1 & ip2, test_data[['patch']] == 2 & ip1)) / 12
		
		pnotfound <- 12 - nrow(test_data)
		
		result_SR <- rbind(result_SR, data.table(rho = rho, eps = eps, tp1 = rng[1], tp2 = rng[2] - rng[1], 
							 pnest = sum(inest)/12, 
							 pfood = sum(c(ip1, ip2))/12, pnotfound = pnotfound,
							 pinterpatch = pinterpatch))
	}
	cat(paste0('Iter ',i, ' from ', l, '\r'))
}
Sys.time() - t0

ggplot(data = result_SR[is.finite(tp1), .(tp1 = 1/median(tp1)), by = c('rho', 'eps')],
       aes(rho, eps, fill = tp1)) +
	geom_tile() + scico::scale_fill_scico('FPT efficiency') +
	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+
	ggtitle('SR ignore social feedbacks')

ggplot(data = result_SR[is.finite(tp2), .(tp2 = 1/median(tp2)), by = c('rho', 'eps')],
       aes(rho, eps, fill = tp2)) +
	geom_tile() + scico::scale_fill_scico('Exploitation efficiency') +
	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+
	ggtitle('SR ignore social feedbacks')


result_LR <- c()
l <- length(files_LR)
t0 <- Sys.time()
for(i in seq_along(files_LR)){
	
	test_data <- data.table(read_parquet(paste0(path_LR, files_LR[i])))
	test_data[['patch']] <- patch_index
	if(nrow(test_data[is.finite(t)])){
		ns_data <- data.table(read_parquet(paste0(path_SR, gsub('_food','', files_LR[i]))))
		test_data <- test_data[is.finite(t), .(node = revert_node(node), t = t, origin = revert_node(origin))]
		f1 <- min(ns_data[N > 0, Frame])
		splt <- strsplit(files_LR[i], '_')[[1]]
		rho <- round(as.numeric(splt[2]), 2)
		eps <- round(as.numeric(splt[4]), 2)
		rng <- range(test_data[['t']]) - f1
		inest <- test_data[['origin']] == '(0, 22)'
		ip1 <- test_data[['origin']] %in% patch_1
		ip2 <- test_data[['origin']] %in% patch_2
		
		pinterpatch <- sum(cbind(test_data[['patch']] == 1 & ip2, test_data[['patch']] == 2 & ip1)) / 12
		
		pnotfound <- 12 - nrow(test_data)
		
		result_LR <- rbind(result_LR, data.table(rho = rho, eps = eps, tp1 = rng[1], tp2 = rng[2] - rng[1], 
							 pnest = sum(inest)/12, 
							 pfood = sum(c(ip1, ip2))/12, pnotfound = pnotfound,
							 pinterpatch = pinterpatch))
	}
	cat(paste0('Iter ',i, ' from ', l, '\r'))
}
Sys.time() - t0

ggplot(data = result_LR[is.finite(tp1), .(tp1 = 1/median(tp1)), by = c('rho', 'eps')],
       aes(rho, eps, fill = tp1)) +
	geom_tile() + scico::scale_fill_scico('FPT efficiency') +
	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+
	ggtitle('LR ignore social feedbacks')

ggplot(data = result_LR[is.finite(tp2), .(tp2 = 1/median(tp2)), by = c('rho', 'eps')],
       aes(rho, eps, fill = tp2)) +
	geom_tile() + scico::scale_fill_scico('Exploitation efficiency') +
	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+
	ggtitle('LR ignore social feedbacks')


result_both <- c()
l <- length(files_both)
t0 <- Sys.time()
for(i in seq_along(files_both)){
	
	test_data <- data.table(read_parquet(paste0(path_both, files_both[i])))
	test_data[['patch']] <- patch_index
	if(nrow(test_data[is.finite(t)])){
		ns_data <- data.table(read_parquet(paste0(path_SR, gsub('_food','', files_both[i]))))
		test_data <- test_data[is.finite(t), .(node = revert_node(node), t = t, origin = revert_node(origin))]
		f1 <- min(ns_data[N > 0, Frame])
		splt <- strsplit(files_both[i], '_')[[1]]
		rho <- round(as.numeric(splt[2]), 2)
		eps <- round(as.numeric(splt[4]), 2)
		rng <- range(test_data[['t']]) - f1
		inest <- test_data[['origin']] == '(0, 22)'
		ip1 <- test_data[['origin']] %in% patch_1
		ip2 <- test_data[['origin']] %in% patch_2
		
		pinterpatch <- sum(cbind(test_data[['patch']] == 1 & ip2, test_data[['patch']] == 2 & ip1)) / 12
		
		pnotfound <- 12 - nrow(test_data)
		
		result_both <- rbind(result_both, data.table(rho = rho, eps = eps, tp1 = rng[1], tp2 = rng[2] - rng[1], 
							 pnest = sum(inest)/12, 
							 pfood = sum(c(ip1, ip2))/12, pnotfound = pnotfound,
							 pinterpatch = pinterpatch))
	}
	cat(paste0('Iter ',i, ' from ', l, '\r'))
}
Sys.time() - t0

ggplot(data = result_both[is.finite(tp1), .(tp1 = 1/median(tp1)), by = c('rho', 'eps')],
       aes(rho, eps, fill = tp1)) +
	geom_tile() + scico::scale_fill_scico('FPT efficiency') +
	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+
	ggtitle('Both listen to social feedbacks')

ggplot(data = result_both[is.finite(tp2), .(tp2 = 1/median(tp2)), by = c('rho', 'eps')],
       aes(rho, eps, fill = tp2)) +
	geom_tile() + scico::scale_fill_scico('Exploitation efficiency') +
	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+
	ggtitle('Both listen to social feedbacks')

ggplot(data = result_both[is.finite(pnotfound), .(pnotfound = 1-mean(pnotfound)), by = c('rho', 'eps')],
       aes(rho, eps, fill = pnotfound)) +
	geom_tile() + scico::scale_fill_scico('Probability of finding food', limits = c(0, 1)) +
	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+
	ggtitle('Both listen to social feedbacks')

ggplot(data = result_both[is.finite(pnest), .(pnest = mean(pnest)), by = c('rho', 'eps')],
       aes(rho, eps, fill = pnest)) +
	geom_tile() + scico::scale_fill_scico('Probability of nest to patch', limits = c(0, 1)) +
	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+
	ggtitle('Both listen to social feedbacks')

ggplot(data = result_both[is.finite(pfood), .(pfood = mean(pfood)), by = c('rho', 'eps')],
       aes(rho, eps, fill = pfood)) +
	scico::scale_fill_scico('Probability of patch to patch') +
	geom_tile() + #scico::scale_fill_scico('Probability of patch to patch', limits = c(0, 1)) +
	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+
	ggtitle('Both listen to social feedbacks')


# ggplot(data = result_both[, .(a = length(pnotfound)), by = c('rho', 'eps')],
#        aes(rho, eps, fill = a)) +
# 	geom_tile() + scico::scale_fill_scico('Exploitation efficiency') +
# 	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+
# 	ggtitle('Both listen to social feedbacks')

median_tp1 <- rbindlist(list(LR = result_LR[is.finite(tp1), .(tp1 = 1/median(tp1)), by = c('rho', 'eps')],
			     SR = result_SR[is.finite(tp1), .(tp1 = 1/median(tp1)), by = c('rho', 'eps')],
			     both = result_both[is.finite(tp1), .(tp1 = 1/median(tp1)), by = c('rho', 'eps')]), idcol = TRUE)
ggplot(data = median_tp1,
       aes(rho, eps, fill = tp1)) +
	geom_raster(interpolate = TRUE) + scale_fill_viridis('Exploration efficiency',
							     limits = c(0.0009, 0.0028),
							     option = 'C') +
	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+ 
	facet_wrap(~.id)


median_tp2 <- rbindlist(list(LR = result_LR[is.finite(tp2), .(tp2 = 1/median(tp2)), by = c('rho', 'eps')],
	       SR = result_SR[is.finite(tp2), .(tp2 = 1/median(tp2)), by = c('rho', 'eps')],
	       both = result_both[is.finite(tp2), .(tp2 = 1/median(tp2)), by = c('rho', 'eps')]), idcol = TRUE)
ggplot(data = median_tp2,
       aes(rho, eps, fill = tp2)) +
	geom_raster(interpolate = TRUE) + scale_fill_viridis('Exploitation efficiency', 
					      limits = c(4.5e-4, 9e-4), breaks = seq(4e-4,8e-4, 2e-4),
					      option = 'C') +
	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+ 
	facet_wrap(~.id)

### ALL TOGETHER ####
ggarrange(
ggplot(data = result_both[is.finite(tp2), .(tp2 = 1/median(tp2)), by = c('rho', 'eps')],
       aes(rho, eps, fill = tp2)) +
	geom_tile() + scico::scale_fill_scico('Exploitation efficiency') +
	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+
	ggtitle('Both listen to social feedbacks'),
ggplot(data = result_SR[is.finite(tp2), .(tp2 = 1/median(tp2)), by = c('rho', 'eps')],
       aes(rho, eps, fill = tp2)) +
	geom_tile() + scico::scale_fill_scico('Exploitation efficiency') +
	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+
	ggtitle('SR ignore social feedbacks'),
ggplot(data = result_LR[is.finite(tp2), .(tp2 = 1/median(tp2)), by = c('rho', 'eps')],
       aes(rho, eps, fill = tp2)) +
	geom_tile() + scico::scale_fill_scico('Exploitation efficiency') +
	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+
	ggtitle('LR ignore social feedbacks')
)

ggplot(data = result_LR[is.finite(tp2), .(tp2 = 1/median(tp2)), by = c('rho', 'eps')],
       aes(rho, eps, fill = tp2)) +
	geom_raster(interpolate = TRUE) + scico::scale_fill_scico('Exploitation efficiency') +
	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+
	ggtitle('LR ignore social feedbacks')

# ggarrange(
# 	ggplot(data = result_both[is.finite(tp2), .(tp2 = 1/mean(tp2)), by = c('rho', 'eps')],
# 	       aes(rho, eps, fill = tp2)) +
# 		geom_tile() + scico::scale_fill_scico('Exploitation efficiency') +
# 		xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+
# 		ggtitle('Both listen to social feedbacks'),
# 	ggplot(data = result_SR[is.finite(tp2), .(tp2 = 1/mean(tp2)), by = c('rho', 'eps')],
# 	       aes(rho, eps, fill = tp2)) +
# 		geom_tile() + scico::scale_fill_scico('Exploitation efficiency') +
# 		xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+
# 		ggtitle('SR ignore social feedbacks'),
# 	ggplot(data = result_LR[is.finite(tp2), .(tp2 = 1/mean(tp2)), by = c('rho', 'eps')],
# 	       aes(rho, eps, fill = tp2)) +
# 		geom_tile() + scico::scale_fill_scico('Exploitation efficiency') +
# 		xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+
# 		ggtitle('LR ignore social feedbacks')
# )

# ggarrange(
# 	ggplot(data = result_both[is.finite(tp1), .(tp1 = 1/median(tp1)), by = c('rho', 'eps')],
# 	       aes(eps, tp1, color = factor(rho), group = factor(rho))) +
# 		geom_path(linewidth = 1.5) +
# 		geom_point(size = 5, shape = 21, fill = NA)+
# 		scale_color_viridis_d('Proportion of LR', end = 0.9)+
# 		xlab(TeX('Proportion of listeners'))+ylab('Exploitation efficiency')+
# 		ggtitle('Both listen to social feedbacks')
# 		# geom_hline(yintercept = det_tps_med[['tp1']], linetype = 2, linewidth = 1.2, color = 'grey80')
# 	,
# 	ggplot(data = result_SR[is.finite(tp1), .(tp1 = 1/median(tp1)), by = c('rho', 'eps')],
# 	       aes(eps, tp1, color = factor(rho), group = factor(rho))) +
# 		geom_path(linewidth = 1.5) +
# 		geom_point(size = 5, shape = 21, fill = NA)+
# 		scale_color_viridis_d('Proportion of LR', end = 0.9)+
# 		xlab(TeX('Proportion of liste ners'))+ylab('Exploitation efficiency')+
# 		ggtitle('SR ignore social feedbacks')
# 		# geom_hline(yintercept = det_tps_med[['tp1']], linetype = 2, linewidth = 1.2, color = 'grey80')
# 	,
# 	ggplot(data = result_LR[is.finite(tp1), .(tp1 = 1/median(tp1)), by = c('rho', 'eps')],
# 	       aes(eps, tp1, color = factor(rho), group = factor(rho))) +
# 		geom_path(linewidth = 1.5) +
# 		geom_point(size = 5, shape = 21, fill = NA)+
# 		scale_color_viridis_d('Proportion of LR', end = 0.9)+
# 		xlab(TeX('Proportion of listeners'))+ylab('Exploitation efficiency')+
# 		ggtitle('LR ignore social feedbacks'),
# 		# geom_hline(yintercept = det_tps_med[['tp1']], linetype = 2, linewidth = 1.2, color = 'grey80'),
# 	ncol = 3, nrow = 1, common.legend = TRUE)

ggarrange(
ggplot(data = result_both[is.finite(tp2), .(tp2 = 1/mean(tp2)), by = c('rho', 'eps')],
       aes(eps, tp2, color = factor(rho), group = factor(rho))) +
	geom_path(linewidth = 1.5) +
	geom_point(size = 5, shape = 21, fill = NA)+
	scale_color_viridis_d('Proportion of LR', end = 0.9)+
	xlab(TeX('Proportion of listeners'))+ylab('Exploitation efficiency')+
	ggtitle('Both listen to social feedbacks')+
	geom_hline(yintercept = det_tps[['tp2']], linetype = 2, linewidth = 1.2, color = 'grey80')
,
ggplot(data = result_SR[is.finite(tp2), .(tp2 = 1/mean(tp2)), by = c('rho', 'eps')],
       aes(eps, tp2, color = factor(rho), group = factor(rho))) +
	geom_path(linewidth = 1.5) +
	geom_point(size = 5, shape = 21, fill = NA)+
	scale_color_viridis_d('Proportion of LR', end = 0.9)+
	xlab(TeX('Proportion of listeners'))+ylab('Exploitation efficiency')+
	ggtitle('SR ignore social feedbacks')+
	geom_hline(yintercept = det_tps[['tp2']], linetype = 2, linewidth = 1.2, color = 'grey80')
,
ggplot(data = result_LR[is.finite(tp2), .(tp2 = 1/mean(tp2)), by = c('rho', 'eps')],
       aes(eps, tp2, color = factor(rho), group = factor(rho))) +
	geom_path(linewidth = 1.5) +
	geom_point(size = 5, shape = 21, fill = NA)+
	scale_color_viridis_d('Proportion of LR', end = 0.9)+
	xlab(TeX('Proportion of listeners'))+ylab('Exploitation efficiency')+
	ggtitle('LR ignore social feedbacks')+
	geom_hline(yintercept = det_tps[['tp2']], linetype = 2, linewidth = 1.2, color = 'grey80'),
ncol = 3, nrow = 1, common.legend = TRUE)

ggarrange(
	ggplot(data = result_both[is.finite(pnest), .(pnest = mean(pnest)), by = c('rho', 'eps')],
	       aes(eps, pnest, color = factor(rho), group = factor(rho))) +
		geom_path(linewidth = 1.5) +
		geom_point(size = 5, shape = 21, fill = NA)+
		scale_color_viridis_d('Proportion of LR', end = 0.9)+
		xlab(TeX('Proportion of listeners'))+ylab('Probability of nest patch')+
		ggtitle('Both listen to social feedbacks')
	,
	ggplot(data = result_SR[is.finite(pnest), .(pnest = mean(pnest)), by = c('rho', 'eps')],
	       aes(eps, pnest, color = factor(rho), group = factor(rho))) +
		geom_path(linewidth = 1.5) +
		geom_point(size = 5, shape = 21, fill = NA)+
		scale_color_viridis_d('Proportion of LR', end = 0.9)+
		xlab(TeX('Proportion of listeners'))+ylab('Probability of nest patch')+
		ggtitle('SR ignore social feedbacks')
	,
	ggplot(data = result_LR[is.finite(pnest), .(pnest = mean(pnest)), by = c('rho', 'eps')],
	       aes(eps, pnest, color = factor(rho), group = factor(rho))) +
		geom_path(linewidth = 1.5) +
		geom_point(size = 5, shape = 21, fill = NA)+
		scale_color_viridis_d('Proportion of LR', end = 0.9)+
		xlab(TeX('Proportion of listeners'))+ylab('Probability of nest patch')+
		ggtitle('LR ignore social feedbacks'),
	ncol = 3, nrow = 1, common.legend = TRUE)

ggarrange(
	ggplot(data = result_both[is.finite(pfood), .(pfood = mean(pfood)), by = c('rho', 'eps')],
	       aes(eps, pfood, color = factor(rho), group = factor(rho))) +
		geom_path(linewidth = 1.5) +
		geom_point(size = 5, shape = 21, fill = NA)+
		scale_color_viridis_d('Proportion of LR', end = 0.9)+
		xlab(TeX('Proportion of listeners'))+ylab('Probability of nest patch')+
		ggtitle('Both listen to social feedbacks')
	,
	ggplot(data = result_SR[is.finite(pfood), .(pfood = mean(pfood)), by = c('rho', 'eps')],
	       aes(eps, pfood, color = factor(rho), group = factor(rho))) +
		geom_path(linewidth = 1.5) +
		geom_point(size = 5, shape = 21, fill = NA)+
		scale_color_viridis_d('Proportion of LR', end = 0.9)+
		xlab(TeX('Proportion of listeners'))+ylab('Probability of nest patch')+
		ggtitle('SR ignore social feedbacks')
	,
	ggplot(data = result_LR[is.finite(pfood), .(pfood = mean(pfood)), by = c('rho', 'eps')],
	       aes(eps, pfood, color = factor(rho), group = factor(rho))) +
		geom_path(linewidth = 1.5) +
		geom_point(size = 5, shape = 21, fill = NA)+
		scale_color_viridis_d('Proportion of LR', end = 0.9)+
		xlab(TeX('Proportion of listeners'))+ylab('Probability of nest patch')+
		ggtitle('LR ignore social feedbacks'),
	ncol = 3, nrow = 1, common.legend = TRUE)

ggarrange(
	ggplot(data = result_both[is.finite(pnotfound), .(pnotfound = mean(1-pnotfound)), by = c('rho', 'eps')],
	       aes(eps, pnotfound, color = factor(rho), group = factor(rho))) +
		geom_path(linewidth = 1.5) +
		geom_point(size = 5, shape = 21, fill = NA)+
		scale_color_viridis_d('Proportion of LR', end = 0.9)+
		xlab(TeX('Proportion of listeners'))+ylab('Probability of nest patch')+
		ggtitle('Both listen to social feedbacks')
	,
	ggplot(data = result_SR[is.finite(pnotfound), .(pnotfound = mean(1-pnotfound)), by = c('rho', 'eps')],
	       aes(eps, pnotfound, color = factor(rho), group = factor(rho))) +
		geom_path(linewidth = 1.5) +
		geom_point(size = 5, shape = 21, fill = NA)+
		scale_color_viridis_d('Proportion of LR', end = 0.9)+
		xlab(TeX('Proportion of listeners'))+ylab('Probability of nest patch')+
		ggtitle('SR ignore social feedbacks')
	,
	ggplot(data = result_LR[is.finite(pnotfound), .(pnotfound = mean(1-pnotfound)), by = c('rho', 'eps')],
	       aes(eps, pnotfound, color = factor(rho), group = factor(rho))) +
		geom_path(linewidth = 1.5) +
		geom_point(size = 5, shape = 21, fill = NA)+
		scale_color_viridis_d('Proportion of LR', end = 0.9)+
		xlab(TeX('Proportion of listeners'))+ylab('Probability of nest patch')+
		ggtitle('LR ignore social feedbacks'),
	ncol = 3, nrow = 1, common.legend = TRUE)

## medians
ggarrange(
	ggplot(data = result_both[is.finite(tp2), .(tp2 = 1/median(tp2)), by = c('rho', 'eps')],
	       aes(eps, tp2, color = factor(rho), group = factor(rho))) +
		geom_path(linewidth = 1.5) +
		geom_point(size = 5, shape = 21, fill = NA)+
		scale_color_viridis_d('Proportion of LR', end = 0.9)+
		xlab(TeX('Proportion of listeners'))+ylab('Exploitation efficiency')+
		ggtitle('Both listen to social feedbacks')+
		geom_hline(yintercept = det_tps_med[['tp2']], linetype = 2, linewidth = 1.2, color = 'grey80')
	,
	ggplot(data = result_SR[is.finite(tp2), .(tp2 = 1/median(tp2)), by = c('rho', 'eps')],
	       aes(eps, tp2, color = factor(rho), group = factor(rho))) +
		geom_path(linewidth = 1.5) +
		geom_point(size = 5, shape = 21, fill = NA)+
		scale_color_viridis_d('Proportion of LR', end = 0.9)+
		xlab(TeX('Proportion of listeners'))+ylab('Exploitation efficiency')+
		ggtitle('SR ignore social feedbacks')+
		geom_hline(yintercept = det_tps_med[['tp2']], linetype = 2, linewidth = 1.2, color = 'grey80'),
	ggplot(data = result_LR[is.finite(tp2), .(tp2 = 1/median(tp2)), by = c('rho', 'eps')],
	       aes(eps, tp2, color = factor(rho), group = factor(rho))) +
		geom_path(linewidth = 1.5) +
		geom_point(size = 5, shape = 21, fill = NA)+
		scale_color_viridis_d('Proportion of LR', end = 0.9)+
		xlab(TeX('Proportion of listeners'))+ylab('Exploitation efficiency')+
		ggtitle('LR ignore social feedbacks')+
		geom_hline(yintercept = det_tps_med[['tp2']], linetype = 2, linewidth = 1.2, color = 'grey80')
, common.legend = TRUE)

## SOME SUBSET OF RHOs
# ggplot(data = result_both[is.finite(tp2) & rho %in% seq(0, 1, 0.2), .(tp2 = 1/mean(tp2)), by = c('rho', 'eps')],
#        aes(eps, tp2, color = factor(rho), group = factor(rho))) +
# 	geom_path(linewidth = 1.5) +
# 	geom_point(size = 5, shape = 21, fill = NA)+
# 	scale_color_viridis_d('Proportion of LR', end = 0.9)+
# 	xlab(TeX('Proportion of listeners'))+ylab('Exploitation efficiency')+
# 	ggtitle('Both listen to social feedbacks')
# 
# ggplot(data = result_SR[is.finite(tp2) & rho %in% seq(0, 1, 0.2), .(tp2 = 1/mean(tp2)), by = c('rho', 'eps')],
#        aes(eps, tp2, color = factor(rho), group = factor(rho))) +
# 	geom_path(linewidth = 1.5) +
# 	geom_point(size = 5, shape = 21, fill = NA)+
# 	scale_color_viridis_d('Proportion of LR', end = 0.9)+
# 	xlab(TeX('Proportion of listeners'))+ylab('Exploitation efficiency')+
# 	ggtitle('SR ignore social feedbacks')
