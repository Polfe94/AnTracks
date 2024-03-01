source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
load('~/research/gits/AnTracks/results/mov_rho_R.RData')

library(arrow)

# path <- '/home/polfer/research/gits/AutomatAnts/results/2024/movement_results/'
# f <- list.files(path)
# files <- f[grepl('data', f)]
# l <- length(files)
# N <- 10000L

# result_list <- vector('list', l)
# for(i in seq_along(files)[130:(l)]){
# 	gc()
# 	n <- strsplit(files[i], '_')[[1]]
# 	rho <- round(as.numeric(n[2]), 2)
# 	R <- round(as.numeric(n[4]), 2)
# 	dt <- data.table(read_parquet(paste0(path, files[i])))[, .N, by = 'id']
# 	type <- rep('LR', N * rho)
# 	type <- c(type, rep('SR', N - length(type)))
# 	dt[['rho']] <- rho
# 	dt[['R']] <- R
# 	dt[['scout']] <- type
# 	cat(paste0('Iteration ', formatC(i, digits = 2, flag = '0'), ' from ', l, '\r'))
# 	result_list[[i]] <- dt
# }

# mov_rho_R <- rbindlist(result_list)[order(rho, R)]
		
ggplot(data = mov_rho_R[, .(N = mean(N)), by = c('rho', 'R')], aes(rho, R, fill = N)) + geom_tile()+
	ylab('Radius') +
	scico::scale_fill_scico('FPT', palette = 'davos', end = 0.8)
ggplot(data = mov_rho_R[, .(N = median(N)), by = c('rho', 'R')], aes(rho, R, fill = N)) + geom_tile()
ggplot(data = mov_rho_R[, .(N = quantile(N, q = 0.9)), by = c('rho', 'R')], aes(rho, R, fill = N)) + geom_tile()

ggplot(data = mov_rho_R[rho %in% c(0, 1)], aes(factor(scout), N, fill = factor(scout)))+
	geom_boxplot(alpha = 0.6, outlier.shape = NA) +
	ylab('Mean first passage (Iterations)') + xlab('')+
	facet_wrap(~ factor(R), scales = 'free', 
		   labeller = as_labeller(function(i){paste0('R = ',i)}))+
	scale_fill_manual('', values = c('mediumpurple', 'gold3'))

sapply(unique(sort(mov_rho_R[['R']])), function(i) wilcox.test(N ~ rho, data = mov_rho_R[rho %in% c(0, 1) & R == i])$p.value)
sapply(unique(sort(mov_rho_R[['R']])), function(i) t.test(N ~ rho, data = mov_rho_R[rho %in% c(0, 1) & R == i])$p.value)

# ggplot(data = mov_rho_R[rho %in% seq(0, 1, 0.25), .(N = mean(N)), by = c('rho', 'R')], aes(R, N, color = factor(rho))) +
# 	geom_point(size = 3) + geom_path(show.legend = FALSE)+
# 	ylab('Mean first passage (Iterations)') + xlab('Radius')+
# 	scale_color_viridis_d('Proportion of LR')

ggplot(data = mov_rho_R[rho %in% seq(0, 1, 0.25), .(N = mean(N)), by = c('rho', 'R')],
       aes(R, N, color = factor(rho))) +
	geom_point(size = 5) + geom_path(show.legend = FALSE, linewimov_rho_Rh = 1.25)+
	ylab('Mean first passage (Iterations)') + xlab('Radius')+
	scale_color_viridis_d('Proportion of LR')+
	theme(legend.position = c(0.2, 2/3), 
	      legend.background = element_rect(fill =NA, color = 'black'))


ggplot(data = mov_rho_R[, .(N = mean(N), sd = sd(N)), by = c('scout', 'R')],
       aes(R, N, color = factor(scout))) +
	geom_pointrange(aes(ymin = N - sd, ymax = N + sd), size = 2, alpha = 0.75,
			linewimov_rho_Rh = 3) +
	geom_path(show.legend = FALSE)+
	ylab('Mean first passage (Iterations)') + xlab('Radius')+
	scale_color_viridis_d('Scout behaviour', direction = -1, end = 0.85)+
	theme(legend.position = c(0.2, 2/3), 
	      legend.background = element_rect(fill =NA, color = 'black'))

# ggplot(data = mov_rho_R[rho %in% seq(0, 1, 0.25)],
#        aes(R, N, fill = factor(scout))) +
# 	geom_boxplot() + 
# 	ylab('Mean first passage (Iterations)') + xlab('Radius')+
# 	scale_fill_viridis_d('Proportion of LR')+
# 	facet_wrap(~ factor(rho))

ggplot(data = mov_rho_R[, .(N = mean(N)), by = c('rho', 'R', 'scout')], aes(rho, R, fill = N)) + geom_tile()+
	facet_wrap(~scout)

ggplot(data = mov_rho_R[rho %in% c(0, 1)], aes(factor(rho), N, fill = factor(scout)))+
	geom_boxplot() + facet_wrap(~ factor(R), scales = 'free')
