source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')

load('~/research/gits/AnTracks/data/det.RData')
path <- '/home/polfer/research/gits/AutomatAnts/results/2024/prob_rec/'


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
ggplot(data = static_tps[is.finite(tp1), .(tp2 = 1/mean(tp2, na.rm = TRUE)), by = 'rho'],
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
	geom_smooth(formula = y ~ poly(x,2), method = 'lm')+
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
