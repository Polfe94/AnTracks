source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
path <- '~/research/gits/AutomatAnts/results/2024/LR_det_rho/'
load('~/research/gits/AnTracks/data/det.RData')

det_sync <- vapply(det, function(i){
	x <- sapply(i@food, function(t) min(t[['t']]))
	(max(x) - min(x)) / 2
}, numeric(1))
det_tp2 <- vapply(det, function(i){
	x <- range(rbindlist(i@food)[['t']])
	(max(x) - min(x)) / 2
}, numeric(1))
det_tp1 <- vapply(det, function(i){
	x <- range(rbindlist(i@food)[['t']])
	(min(x)) / 2
}, numeric(1))
det_sync <- data.table(tp2_sincr = det_sync, tp1 = det_tp1,tp2 = det_tp2, rho = 'exp')

f <- list.files(path)

library(arrow)
library(latex2exp)
l <- length(f[grepl('food', f)])
static_tps <- rbindlist(lapply(f[grepl('food', f)], function(x){
	food <- data.table(read_parquet(paste0(path, x)))[, .(node = revert_node(node),
									 t = t, origin = origin)]
	if(any(is.finite(food[['t']]))){
		food[['patch']] <- c(rep(1, 6), rep(2, 6))
		p <- food[['origin']]
		rho <- round(as.numeric(strsplit(x, '_')[[1]][2]), 2)
		idx <- which.min(food[['t']])
		p1 <- food[idx, patch]
		mint_p1 <- food[idx, t]
		mint_p2 <- food[patch != p1, min(t)]
		tp2_sincr <- mint_p2 - mint_p1
		tp2_min <- min(mint_p2 - food[patch == p1 & t < mint_p2, t], na.rm = TRUE)
		cat(paste0('Finished iteration ', formatC(which(x == f[grepl('food', f)]), flag = '0', digits = 3), ' from ', l, '\r'))
		data.table(tp1 = mint_p1, 
			   tp2 = max(food[['t']], na.rm = TRUE) - mint_p1,
			   tp2_sincr = tp2_sincr, tp2_min = tp2_min, 
			   tp_diff = mean(diff(food[, sort(t)]), na.rm = TRUE),
			   pnest = sum(p == 'nest'), pfood = sum(p == 'food'), pnotfound = sum(!is.finite(food[['t']])),
			   rho = rho)
	}
}))

ggplot(data = static_tps[, .(tp1 = 1/mean(tp1, na.rm = TRUE)), by = 'rho'],
       aes(factor(rho), tp1)) + geom_point()
ggplot(data = static_tps[, .(tp2 = 1/mean(tp2, na.rm = TRUE)), by = 'rho'],
       aes(factor(rho), tp2)) + geom_point()

# ggplot(data = static_tps,
#        aes(factor(rho), tp1)) + geom_boxplot()

## 

ggplot(data = rbind(det_sync[, c('rho', 'tp1')], static_tps[rho %in% c(0, 1), c('rho', 'tp1')]),
       aes(factor(rho, levels = c('exp', '0', '1'),
       	   labels = c('Experiments', '0', '1')), 1/tp1, 
           fill = factor(rho, levels = c('exp', '0', '1')))) + 
	geom_boxplot(show.legend = FALSE, alpha = 0.6, outlier.shape = NA)+
	geom_jitter(size = 3, alpha = 0.25, show.legend = FALSE)+
	ylab(TeX('Exploration time ($s^{-1}$)'))+
	xlab(TeX('Proportion of LR $(\\rho)$'))+
	scale_fill_manual('',values = c('mediumpurple', 'gold3', 'brown4'))

# ggplot(data = static_tps[rho %in% c(0, 1)],
#        aes(factor(rho), 1/tp1, fill = factor(rho))) + 
# 	geom_boxplot(show.legend = FALSE, alpha = 0.6, outlier.shape = NA)+
# 	geom_jitter(size = 3, alpha = 0.25, show.legend = FALSE)+
# 	ylab(TeX('Exploration time ($s^{-1}$)'))+
# 	xlab(TeX('Proportion of LR $(\\rho)$'))+
# 	scale_fill_manual('', values = c('mediumpurple', 'gold3'))

ggplot(data = static_tps,
       aes(factor(rho), tp2)) + geom_boxplot()

ggplot(data = rbind(det_sync[, c('rho', 'tp2')], static_tps[rho %in% c(0, 1), c('rho', 'tp2')]),
       aes(factor(rho, levels = c('exp', '0', '1'),
       	   labels = c('Experiments', '0', '1')), 1/tp2, 
           fill = factor(rho, levels = c('exp', '0', '1')))) + 
	geom_boxplot(show.legend = FALSE, alpha = 0.6, outlier.shape = NA)+
	geom_jitter(size = 3, alpha = 0.25, show.legend = FALSE)+
	ylab(TeX('Exploitation time ($s^{-1}$)'))+
	xlab(TeX('Proportion of LR $(\\rho)$'))+
	scale_fill_manual('',values = c('mediumpurple', 'gold3', 'brown4'))

# ggplot(data = static_tps[rho %in% c(0, 1)],
#        aes(factor(rho), 1/tp2, fill = factor(rho))) + 
# 	geom_boxplot(show.legend = FALSE, alpha = 0.6, outlier.shape = NA)+
# 	geom_jitter(size = 3, alpha = 0.25, show.legend = FALSE)+
# 	ylab(TeX('Exploitation time ($s^{-1}$)'))+
# 	xlab(latex2exp::TeX('Proportion of LR $(\\rho)$'))+
# 	scale_fill_manual('', values = c('mediumpurple', 'gold3'))

ggplot(data = static_tps[is.finite(tp1), .(tp = 1/mean(tp1)), by = 'rho'],
       aes(rho, tp)) + geom_point() + 
	# geom_hline(yintercept = 1/median(det_sync[['tp1']]), linetype = 'dashed')+
	geom_smooth(formula = y ~ x, method = 'lm')+
	ylab(TeX('Patch discovery ($s^{-1}$)'))+
	xlab(latex2exp::TeX('Proportion of LR $(\\rho)$'))

ggplot(data = static_tps[is.finite(tp2_sincr), .(tp = 1/mean(tp2_sincr)), by = 'rho'],
       aes(rho, tp)) + geom_point() + 
	geom_hline(yintercept = 1/mean(det_sync[['tp2_sincr']]), linetype = 'dashed')+
	geom_smooth(formula = y ~ x, method = 'lm')+
	ylab(TeX('Patch syncronization ($s^{-1}$)'))+
	xlab(latex2exp::TeX('Proportion of LR $(\\rho)$'))

ggplot(data = static_tps[is.finite(tp2), .(tp = 1/mean(tp2)), by = 'rho'],
       aes(rho, tp)) + geom_point() + 
	geom_hline(yintercept = 1/mean(det_sync[['tp2']]), linetype = 'dashed')+
	geom_smooth(formula = y ~ poly(x, 2), method = 'lm')+
	ylab(TeX('Exploitation time ($s^{-1}$)'))+
	xlab(latex2exp::TeX('Proportion of LR $(\\rho)$'))

ggplot(data = static_tps[, .(pnest = mean(pnest/12, na.rm = TRUE)), by = 'rho'],
       aes(factor(rho), pnest)) + geom_point(size = 3)+
	ylab('Probability of descovery from nest') +
	xlab(TeX('Proportion of LR ($\\rho$)'))
ggplot(data = static_tps[, .(pfood = mean(pfood/12, na.rm = TRUE)), by = 'rho'],
       aes(factor(rho), pfood)) + geom_point(size = 3)+
	ylab('Probability of descovery from food') +
	xlab(TeX('Proportion of LR ($\\rho$)'))


ggplot(data = static_tps[, .(p = mean(pfood/pnest, na.rm = TRUE)), by = 'rho'],
       aes(factor(rho), p)) + geom_point(size = 3) +
	
ggplot(data = static_tps[, .(p =pfood/pnest), by = 'rho'],
       aes(factor(rho), p)) + geom_boxplot()
ggplot(data = static_tps[, .(p =pnest/pfood), by = 'rho'],
       aes(factor(rho), p)) + geom_boxplot()

# ggplot(data = static_tps[, .(pnest = mean(pnest, na.rm = TRUE)), by = 'rho'],
#        aes(factor(rho), pnest)) + geom_point()
# ggplot(data = static_tps,
#        aes(factor(rho), pnest)) + geom_boxplot()
# ggplot(data = static_tps[, .(pfood = mean(pfood, na.rm = TRUE)), by = 'rho'],
#        aes(factor(rho), pfood)) + geom_point()
ggplot(data = static_tps,
       aes(factor(rho), pfood)) + geom_boxplot()
ggplot(data = static_tps[, .(pnotfound = mean(pnotfound, na.rm = TRUE)), by = 'rho'],
       aes(factor(rho), pnotfound)) + geom_point()



ggplot(data = rbind(det_sync[, c('rho', 'tp2_sincr')], static_tps[, c('rho', 'tp2_sincr')]),
       aes(factor(rho), 1/tp2_sincr)) + geom_boxplot()+
	coord_cartesian(ylim = c(0, 0.0075))
ggplot(data = static_tps,
       aes(factor(rho), tp2_min)) + geom_boxplot()





ggplot(data = static_tps[is.finite(tp2_min), .(tp = 1/mean(tp2_min)), by = 'rho'],
       aes(rho, tp)) + geom_point() +
	geom_smooth(formula = y ~ poly(x, 2), method = 'lm')+
	ylab('Efficiency (patch discovery to last exploitation)')+
	xlab(TeX('Proportion of LR ($\\rho$)'))


ggplot(data = static_tps[is.finite(tp2_min), .(tp = mean(tp2_min)), by = 'rho'],
       aes(rho, tp)) + geom_point() + geom_path()
ggplot(data = static_tps[is.finite(tp2_sincr), .(tp = mean(tp2_sincr)), by = 'rho'],
       aes(rho, tp)) + geom_point() + geom_path()



ggplot(data = static_tps[is.finite(tp_diff), .(tp = 1/mean(tp_diff)), by = 'rho'],
       aes(rho, tp)) + geom_point() + geom_path()
ggplot(data = static_tps[is.finite(tp_diff), .(tp = median(tp_diff)), by = 'rho'],
       aes(rho, tp)) + geom_point() + geom_path()




path <- '~/research/gits/AutomatAnts/results/2024/LR_det_rho_shutdown/'
f <- list.files(path)

library(arrow)
library(latex2exp)
l <- length(f[grepl('food', f)])
shutdown_tps <- rbindlist(lapply(f[grepl('food', f)], function(x){
	food <- data.table(read_parquet(paste0(path, x)))[, .(node = revert_node(node),
							      t = t, origin = origin)]
	if(any(is.finite(food[['t']]))){
		food[['patch']] <- c(rep(1, 6), rep(2, 6))
		p <- food[['origin']]
		rho <- round(as.numeric(strsplit(x, '_')[[1]][2]), 2)
		idx <- which.min(food[['t']])
		p1 <- food[idx, patch]
		mint_p1 <- food[idx, t]
		mint_p2 <- food[patch != p1, min(t)]
		tp2_sincr <- mint_p2 - mint_p1
		tp2_min <- min(mint_p2 - food[patch == p1 & t < mint_p2, t], na.rm = TRUE)
		cat(paste0('Finished iteration ', formatC(which(x == f[grepl('food', f)]), flag = '0', digits = 3), ' from ', l, '\r'))
		data.table(tp1 = mint_p1, 
			   tp2 = max(food[['t']], na.rm = TRUE) - mint_p1,
			   tp2_sincr = tp2_sincr, tp2_min = tp2_min, 
			   tp_diff = mean(diff(food[, sort(t)]), na.rm = TRUE),
			   pnest = sum(p == 'nest'), pfood = sum(p == 'food'), pnotfound = sum(!is.finite(food[['t']])),
			   rho = rho)
	}
}))


ggplot(data = shutdown_tps[, .(tp1 = 1/mean(tp1, na.rm = TRUE)), by = 'rho'],
       aes(factor(rho), tp1)) + geom_point()
ggplot(data = shutdown_tps[, .(tp2 = 1/mean(tp2, na.rm = TRUE)), by = 'rho'],
       aes(factor(rho), tp2)) + geom_point()

# ggplot(data = shutdown_tps,
#        aes(factor(rho), tp1)) + geom_boxplot()

## 

ggplot(data = rbind(det_sync[, c('rho', 'tp1')], shutdown_tps[rho %in% c(0, 1), c('rho', 'tp1')]),
       aes(factor(rho, levels = c('exp', '0', '1'),
       	   labels = c('Experiments', '0', '1')), 1/tp1, 
           fill = factor(rho, levels = c('exp', '0', '1')))) + 
	geom_boxplot(show.legend = FALSE, alpha = 0.6, outlier.shape = NA)+
	geom_jitter(size = 3, alpha = 0.25, show.legend = FALSE)+
	ylab(TeX('Exploration time ($s^{-1}$)'))+
	xlab(TeX('Proportion of LR $(\\rho)$'))+
	scale_fill_manual('',values = c('mediumpurple', 'gold3', 'brown4'))

# ggplot(data = shutdown_tps[rho %in% c(0, 1)],
#        aes(factor(rho), 1/tp1, fill = factor(rho))) + 
# 	geom_boxplot(show.legend = FALSE, alpha = 0.6, outlier.shape = NA)+
# 	geom_jitter(size = 3, alpha = 0.25, show.legend = FALSE)+
# 	ylab(TeX('Exploration time ($s^{-1}$)'))+
# 	xlab(TeX('Proportion of LR $(\\rho)$'))+
# 	scale_fill_manual('', values = c('mediumpurple', 'gold3'))

ggplot(data = shutdown_tps,
       aes(factor(rho), tp2)) + geom_boxplot()

ggplot(data = rbind(det_sync[, c('rho', 'tp2')], shutdown_tps[rho %in% c(0, 1), c('rho', 'tp2')]),
       aes(factor(rho, levels = c('exp', '0', '1'),
       	   labels = c('Experiments', '0', '1')), 1/tp2, 
           fill = factor(rho, levels = c('exp', '0', '1')))) + 
	geom_boxplot(show.legend = FALSE, alpha = 0.6, outlier.shape = NA)+
	geom_jitter(size = 3, alpha = 0.25, show.legend = FALSE)+
	ylab(TeX('Exploitation time ($s^{-1}$)'))+
	xlab(TeX('Proportion of LR $(\\rho)$'))+
	scale_fill_manual('',values = c('mediumpurple', 'gold3', 'brown4'))

# ggplot(data = shutdown_tps[rho %in% c(0, 1)],
#        aes(factor(rho), 1/tp2, fill = factor(rho))) + 
# 	geom_boxplot(show.legend = FALSE, alpha = 0.6, outlier.shape = NA)+
# 	geom_jitter(size = 3, alpha = 0.25, show.legend = FALSE)+
# 	ylab(TeX('Exploitation time ($s^{-1}$)'))+
# 	xlab(latex2exp::TeX('Proportion of LR $(\\rho)$'))+
# 	scale_fill_manual('', values = c('mediumpurple', 'gold3'))

ggplot(data = shutdown_tps[is.finite(tp1), .(tp = 1/mean(tp1)), by = 'rho'],
       aes(rho, tp)) + geom_point() + 
	# geom_hline(yintercept = 1/median(det_sync[['tp1']]), linetype = 'dashed')+
	geom_smooth(formula = y ~ x, method = 'lm')+
	ylab(TeX('Patch discovery ($s^{-1}$)'))+
	xlab(latex2exp::TeX('Proportion of LR $(\\rho)$'))

ggplot(data = shutdown_tps[is.finite(tp2_sincr), .(tp = 1/mean(tp2_sincr)), by = 'rho'],
       aes(rho, tp)) + geom_point() + 
	geom_hline(yintercept = 1/mean(det_sync[['tp2_sincr']]), linetype = 'dashed')+
	geom_smooth(formula = y ~ x, method = 'lm')+
	ylab(TeX('Patch syncronization ($s^{-1}$)'))+
	xlab(latex2exp::TeX('Proportion of LR $(\\rho)$'))

ggplot(data = shutdown_tps[is.finite(tp2), .(tp = 1/mean(tp2)), by = 'rho'],
       aes(rho, tp)) + geom_point() + 
	geom_hline(yintercept = 1/mean(det_sync[['tp2']]), linetype = 'dashed')+
	geom_smooth(formula = y ~ poly(x, 2), method = 'lm')+
	ylab(TeX('Exploitation time ($s^{-1}$)'))+
	xlab(latex2exp::TeX('Proportion of LR $(\\rho)$'))

ggplot(data = shutdown_tps[, .(pnest = mean(pnest/12, na.rm = TRUE)), by = 'rho'],
       aes(factor(rho), pnest)) + geom_point()
ggplot(data = shutdown_tps[, .(pfood = mean(pfood/12, na.rm = TRUE)), by = 'rho'],
       aes(factor(rho), pfood)) + geom_point()

ggplot(data = shutdown_tps[, .(p = mean(pfood/pnest, na.rm = TRUE)), by = 'rho'],
       aes(factor(rho), p)) + geom_point()
ggplot(data = shutdown_tps[, .(p =pfood/pnest), by = 'rho'],
       aes(factor(rho), p)) + geom_boxplot()
ggplot(data = shutdown_tps[, .(p =pnest/pfood), by = 'rho'],
       aes(factor(rho), p)) + geom_boxplot()

# ggplot(data = shutdown_tps[, .(pnest = mean(pnest, na.rm = TRUE)), by = 'rho'],
#        aes(factor(rho), pnest)) + geom_point()
# ggplot(data = shutdown_tps,
#        aes(factor(rho), pnest)) + geom_boxplot()
# ggplot(data = shutdown_tps[, .(pfood = mean(pfood, na.rm = TRUE)), by = 'rho'],
#        aes(factor(rho), pfood)) + geom_point()
ggplot(data = shutdown_tps,
       aes(factor(rho), pfood)) + geom_boxplot()
ggplot(data = shutdown_tps[, .(pnotfound = mean(pnotfound, na.rm = TRUE)), by = 'rho'],
       aes(factor(rho), pnotfound)) + geom_point()



ggplot(data = rbind(det_sync[, c('rho', 'tp2_sincr')], shutdown_tps[, c('rho', 'tp2_sincr')]),
       aes(factor(rho), 1/tp2_sincr)) + geom_boxplot()+
	coord_cartesian(ylim = c(0, 0.0075))
ggplot(data = shutdown_tps,
       aes(factor(rho), tp2_min)) + geom_boxplot()





ggplot(data = shutdown_tps[is.finite(tp2_min), .(tp = 1/mean(tp2_min)), by = 'rho'],
       aes(rho, tp)) + geom_point() +
	geom_smooth(formula = y ~ poly(x, 2), method = 'lm')


ggplot(data = shutdown_tps[is.finite(tp2_min), .(tp = mean(tp2_min)), by = 'rho'],
       aes(rho, tp)) + geom_point() + geom_path()
ggplot(data = shutdown_tps[is.finite(tp2_sincr), .(tp = mean(tp2_sincr)), by = 'rho'],
       aes(rho, tp)) + geom_point() + geom_path()



ggplot(data = shutdown_tps[is.finite(tp_diff), .(tp = 1/mean(tp_diff)), by = 'rho'],
       aes(rho, tp)) + geom_point() + geom_path()
ggplot(data = shutdown_tps[is.finite(tp_diff), .(tp = median(tp_diff)), by = 'rho'],
       aes(rho, tp)) + geom_point() + geom_path()




library(dtw)
library(dtwclust)
path <- '~/research/gits/AutomatAnts/results/2024/LR_det_rho/'

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
		       	!grepl('data', files) & !grepl('keys', files) & grepl('rho_0.001', files)]
	# files <- files[!grepl('food', files) & !grepl('position', files) &
	# 	       	!grepl('data', files) & !grepl('keys', files) & grepl('rho_1.001', files)]
	
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


load("/home/polfer/research/gits/AutomatAnts/data/zSeries.RData")
f <- list.files(path)
files <- f[!grepl('food', f) & !grepl('position', f) &
	       	!grepl('data', f) & !grepl('keys', f) & grepl('rho_0.001', f)]
data <- process_data(path)
knames <- c('norm', 'late', 'flat', 'low')

d <- rbindlist(lapply(data, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
}))
colnames(d) <- knames
k <- knames[apply(d, 1, which.min)]
plot_clusters(data, k)

full_data <- rbindlist(lapply(files, function(i){
	data.table(read_parquet(paste0(path, i)))
}), idcol = TRUE)[order(Frame)]
full_data_avg <- full_data[, .(N = mean(N)), by = 'Frame']

# foods <- rbindlist(lapply(files, function(i){
# 	data.table(read_parquet(paste0(path, i,'_food.parquet')))[, c('t')]
# }), idcol = TRUE)[, .(mint = min(t), maxt = max(t)), by = '.id'][is.finite(mint),
# 								 lapply(.SD, mean), .SDcols = c('mint', 'maxt')]
foods <- rbindlist(lapply(files, function(i){
	data.table(read_parquet(paste0(path, gsub('.parquet', '',i),'_food.parquet')))[, c('t')]
}), idcol = TRUE)[, .(mint = min(t), maxt = max(t)), by = '.id'][is.finite(mint),
								 lapply(.SD, median), .SDcols = c('mint', 'maxt')]


data_peak <- rbindlist(data, idcol = TRUE)

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


ggplot(data = data_peak, aes(Frame, N)) +
	geom_path(aes(group = .id), color = 'grey80', alpha = 0.2, linewidth = 0.8)+
	geom_path(data = det_N, aes(group = .id), color = 'grey30', alpha = 0.5, linewidth = 0.8)+
	
	geom_path(data = full_data_avg, color = 'black', linewidth = 3, alpha = 0.5)+
	geom_path(data = full_data_avg, aes(color = 'Simulations'), linewidth = 2, alpha = 0.5)+
	geom_path(data = det_avg, color = 'black', linewidth = 3, alpha = 0.5)+
	geom_path(data = det_avg, aes(color = 'Experiments'), linewidth = 2, alpha = 0.5)+
	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50)*120, labels = seq(0, 150, 50))+
	scale_y_continuous('Occupancy')+
	scale_color_manual('',values = c('grey30', 'grey80'))+
	geom_vline(xintercept = unlist(foods)*2, linetype = 2, linewidth = 1)+
	geom_vline(xintercept = unlist(foods_det)*2, linetype = 3, linewidth = 1)+
	theme(legend.position = c(0.8, 0.85),
	      legend.background = element_rect(fill = 'white', color = 'black'),
	      legend.title = element_blank())


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
