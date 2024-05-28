source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')

load("/home/polfer/research/gits/AnTracks/results/data_sims_rho_eps.RData")
result_LR <- data_sims_rho_eps[['LR']]
result_SR <- data_sims_rho_eps[['SR']]
result_both <- data_sims_rho_eps[['both']]

median_tp1 <- rbindlist(list(LR = result_LR[is.finite(tp1), .(tp1 = 1/median(tp1)), by = c('rho', 'eps')],
			     SR = result_SR[is.finite(tp1), .(tp1 = 1/median(tp1)), by = c('rho', 'eps')],
			     both = result_both[is.finite(tp1), .(tp1 = 1/median(tp1)), by = c('rho', 'eps')]), idcol = TRUE)

median_tp1 <- rbindlist(list(LR = enhance_land_res(result_LR[is.finite(tp1), .(tp1 = 1/median(tp1/2)),
							     by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
						   zvar = 'tp1', spar = 0.5, n = 200),
			     SR = enhance_land_res(result_SR[is.finite(tp1), .(tp1 = 1/median(tp1/2)),
			     				by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
			     				zvar = 'tp1', spar = 0.5, n = 200),
			     both = enhance_land_res(result_both[is.finite(tp1), .(tp1 = 1/median(tp1/2)),
			     				  by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
			     			zvar = 'tp1', spar = 0.5, n = 200)), idcol = TRUE)


# a <- ggplot(data = median_tp1,
#        aes(rho, eps, fill = tp1)) +
# 	geom_raster() + geom_contour(linetype = 4, aes(z = tp1), color = 'grey5') +
# 	scale_fill_viridis(TeX('Exploration efficiency ($\\cdot 10^{-3} s^{-1}$)'), option = 'C', breaks = seq(0.00125, 0.00225, 0.00025),
# 			   labels = vapply(seq(0.00125, 0.00225, 0.00025), function(i){
# 			   	as.character(i*1000)
# 			   }, character(1))) +
# 	xlab(TeX('Proportion of LR'))+ylab('')+ 
# 	theme(aspect.ratio = 1, strip.placement = 'outside', strip.background = element_blank(),
# 	      legend.position = 'top', legend.key.width = unit(1.5, 'cm'))+
# 	facet_wrap(~ factor(.id, levels = c('LR', 'SR', 'both')),
# 					   strip.position = 'left', 
# 					   labeller = as_labeller(c(both = 'Sociality | LR+SR', 
# 					   			 SR = 'Sociality | SR', 
# 					   			 LR = 'Sociality | LR')))+
# 	guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))

a <- ggplot(data = median_tp1,
	    aes(rho, eps, fill = tp1)) +
	geom_raster() + geom_contour(linetype = 4, aes(z = tp1), color = 'grey5') +
	scale_fill_viridis(TeX('Exploration\nefficiency ($s^{-1}$)'), 
			   option = 'C', breaks = seq(0.0025, 0.004, length.out = 4),
			   labels = sapply(seq(0.0025, 0.004, length.out = 4), function(i){
			   	# as.character(i*1000)
			   	TeX(paste0(as.character(i*1000), '$\\cdot 10^{-3}$'))
			   })) +
	xlab(TeX('Proportion of LR'))+ylab('')+ 
	theme(aspect.ratio = 1, strip.placement = 'outside', strip.background = element_blank(),
	      legend.position = 'right', legend.key.height = unit(1.5, 'cm'))+
	facet_wrap(~ factor(.id, levels = c('LR', 'SR', 'both')),
		   strip.position = 'left', 
		   labeller = as_labeller(c(both = 'Sociality | LR+SR', 
		   			 SR = 'Sociality | SR', 
		   			 LR = 'Sociality | LR')))+
	guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))



median_tp2 <- rbindlist(list(LR = enhance_land_res(result_LR[is.finite(tp2), .(tp2 = 1/median(tp2/2)),
							     by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
						   zvar = 'tp2', spar = 0.5, n = 200),
			     SR = enhance_land_res(result_SR[is.finite(tp2), .(tp2 = 1/median(tp2/2)),
			     				by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
			     				zvar = 'tp2', spar = 0.5, n = 200),
			     both = enhance_land_res(result_both[is.finite(tp2), .(tp2 = 1/median(tp2/2)),
			     				    by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
			     			zvar = 'tp2', spar = 0.5, n = 200)), idcol = TRUE)


# b <- ggplot(data = median_tp2,
#        aes(rho, eps, fill = tp2)) +
# 	geom_raster() + geom_contour(linetype = 4, aes(z = tp2), color = 'grey5', bins = 10) +
# 	scale_fill_viridis(TeX('Exploitation efficiency ($\\cdot 10^{-4} s^{-1}$)'), option = 'C', 
# 			   breaks = seq(0.0004, 0.0009, 0.0001),
# 			   labels = vapply(seq(0.0004, 0.0009, 0.0001), function(i){
# 			   	as.character(i*10000)
# 			   }, character(1))) +
# 	xlab(TeX('Proportion of LR'))+ylab('')+ 
# 	theme(aspect.ratio = 1, strip.placement = 'outside', strip.background = element_blank(),
# 	      legend.position = 'top', legend.key.width = unit(1.5, 'cm'))+
# 	facet_wrap(~ factor(.id, levels = c('LR', 'SR', 'both')),
# 		   strip.position = 'left', 
# 		   labeller = as_labeller(c(both = 'Sociality | LR+SR', 
# 		   			 SR = 'Sociality | SR', 
# 		   			 LR = 'Sociality | LR')))+
# 	guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))

b <- ggplot(data = median_tp2,
	    aes(rho, eps, fill = tp2)) +
	geom_raster() + geom_contour(linetype = 4, aes(z = tp2), color = 'grey5') +
	scale_fill_viridis(TeX('Exploitation\nefficiency ($s^{-1}$)'), 
			   option = 'C', breaks = seq(0.0008, 0.0016, length.out = 5),
			   labels = sapply(seq(0.0008, 0.0016, length.out = 5), function(i){
			   	# as.character(i*1000)
			   	TeX(paste0(as.character(i*1000), '$\\cdot 10^{-3}$'))
			   })) +
	xlab(TeX('Proportion of LR'))+ylab('')+ 
	theme(aspect.ratio = 1, strip.placement = 'outside', strip.background = element_blank(),
	      legend.position = 'right', legend.key.height = unit(1.5, 'cm'))+
	facet_wrap(~ factor(.id, levels = c('LR', 'SR', 'both')),
		   strip.position = 'left', 
		   labeller = as_labeller(c(both = 'Sociality | LR+SR', 
		   			 SR = 'Sociality | SR', 
		   			 LR = 'Sociality | LR')))+
	guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))




median_tradeoff <- rbindlist(list(LR = enhance_land_res(result_LR[, .(tradeoff = median(tp1/(tp1 +tp2))),
							     by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
						   zvar = 'tradeoff', spar = 0.5, n = 200),
			     SR = enhance_land_res(result_SR[, .(tradeoff = median(tp1/(tp1 +tp2))),
			     				by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
			     				zvar = 'tradeoff', spar = 0.5, n = 200),
			     both = enhance_land_res(result_both[, .(tradeoff = median(tp1/(tp1 +tp2))),
			     				    by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
			     			zvar = 'tradeoff', spar = 0.5, n = 200)), idcol = TRUE)


c <- ggplot(data = median_tradeoff,
       aes(rho, eps, fill = 100*tradeoff)) +
	geom_raster() + geom_contour(linetype = 4, aes(z = tradeoff), color = 'grey5', bins = 10) +
	scale_fill_viridis(TeX('Exploration\nduration (%)'), option = 'C', 
			   breaks = seq(20, 35, length.out = 6),
			   labels = sapply(seq(20, 35, length.out = 6), function(i){
			   	TeX(paste0(as.character(i)))
			   })) +
	xlab(TeX('Proportion of LR'))+ylab('')+ 
	theme(aspect.ratio = 1, strip.placement = 'outside', strip.background = element_blank(),
	      legend.position = 'right', legend.key.height = unit(1.5, 'cm'),
	      )+
	facet_wrap(~ factor(.id, levels = c('LR', 'SR', 'both')),
		   strip.position = 'left', 
		   labeller = as_labeller(c(both = 'Sociality | LR+SR', 
		   			 SR = 'Sociality | SR', 
		   			 LR = 'Sociality | LR')))+
	guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))



png('/home/polfer/research/gits/AnTracks/plots/figures_LiquidBrains/Fig_4.png', 6000,6000, res = 450)
grid.arrange(a+ggtitle('A'), b+ggtitle('B'), c+ggtitle('C'),
	     nrow = 3, ncol = 1)
# ggarrange(a, b, c, nrow = 3, ncol = 1, labels = c('A', 'B', 'C'), font.label = list(face = 'plain'))
dev.off()

median_total <- rbindlist(list(LR = enhance_land_res(result_LR[, .(total = 1/median((tp1 +tp2)/2)),
								  by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
							zvar = 'total', spar = 0.5, n = 200),
				  SR = enhance_land_res(result_SR[, .(total = 1/median((tp1 +tp2)/2)),
				  				by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
				  				zvar = 'total', spar = 0.5, n = 200),
				  both = enhance_land_res(result_both[, .(total = 1/median((tp1 +tp2)/2)),
				  				    by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
				  			zvar = 'total', spar = 0.5, n = 200)), idcol = TRUE)

d <- ggplot(data = median_total,
	    aes(rho, eps, fill = total)) +
	geom_raster() + geom_contour(linetype = 4, aes(z = total), color = 'grey5', bins = 10) +
	scale_fill_viridis(TeX('Foraging\nefficiency ($s^{-1}$)'), option = 'C',
			   breaks = seq(0.0006, 0.001, length.out = 5),
			   labels = sapply(seq(0.0006, 0.001, length.out = 5), function(i){
			   	if(i > 0.00099){
			   		TeX(paste0(as.character(i*1000), '$\\cdot 10^{-3}$'))
			   	} else {
			   		TeX(paste0(as.character(i*10000), '$\\cdot 10^{-4}$'))
			   	}
			   	
			   })) +
	xlab(TeX('Proportion of LR'))+ylab('')+ 
	theme(aspect.ratio = 1, strip.placement = 'outside', strip.background = element_blank(),
	      legend.position = 'right', legend.key.height = unit(1.5, 'cm'),
	)+
	facet_wrap(~ factor(.id, levels = c('LR', 'SR', 'both')),
		   strip.position = 'left', 
		   labeller = as_labeller(c(both = 'Sociality | LR+SR', 
		   			 SR = 'Sociality | SR', 
		   			 LR = 'Sociality | LR')))+
	guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))


png('/home/polfer/research/gits/AnTracks/plots/figures_LiquidBrains/Fig_4.png', 6000*2,8000*2, res = 460*2)
grid.arrange(a+ggtitle('A'), b+ggtitle('B'), c+ggtitle('C'), d +ggtitle('D'),
	     nrow = 4, ncol = 1)
# ggarrange(a, b, c, nrow = 3, ncol = 1, labels = c('A', 'B', 'C'), font.label = list(face = 'plain'))
dev.off()

# 
# 
# 
# ggarrange(
# ggplot(data = median_tp1[.id == 'LR'],
#        aes(rho, eps, fill = tp1)) +
# 	geom_raster(interpolate = TRUE) +
# 	scale_fill_viridis('Exploration efficiency', option = 'C') +
# 	xlab(TeX('Proportion of LR'))+ylab(TeX('Sociality | LR'))+ 
# 	theme(aspect.ratio = 1),
# ggplot(data = median_tp1[.id == 'SR'],
#        aes(rho, eps, fill = tp1)) +
# 	geom_raster(interpolate = TRUE) + scale_fill_viridis('Exploration efficiency',
# 							     option = 'C') +
# 	xlab(TeX('Proportion of LR'))+ylab(TeX('Sociality | SR'))+ 
# 	theme(aspect.ratio = 1),
# ggplot(data = median_tp1[.id == 'both'],
#        aes(rho, eps, fill = tp1)) +
# 	geom_raster(interpolate = TRUE) + scale_fill_viridis('Exploration efficiency',
# 							     option = 'C') +
# 	xlab(TeX('Proportion of LR'))+ylab(TeX('Sociality | LR + SR'))+ 
# 	theme(aspect.ratio = 1)
# , nrow = 1, ncol = 3, common.legend = TRUE, legend = 'top')
# 
# ggarrange(
# 	ggplot(data = median_tp2[.id == 'LR'],
# 	       aes(rho, eps, fill = tp2)) +
# 		geom_raster(interpolate = TRUE) +
# 		scale_fill_viridis('Exploitation efficiency', option = 'C') +
# 		xlab(TeX('Proportion of LR'))+ylab(TeX('Sociality | LR'))+ 
# 		theme(aspect.ratio = 1),
# 	ggplot(data = median_tp2[.id == 'SR'],
# 	       aes(rho, eps, fill = tp2)) +
# 		geom_raster(interpolate = TRUE) + scale_fill_viridis('Exploitation efficiency',
# 								     option = 'C') +
# 		xlab(TeX('Proportion of LR'))+ylab(TeX('Sociality | SR'))+ 
# 		theme(aspect.ratio = 1),
# 	ggplot(data = median_tp2[.id == 'both'],
# 	       aes(rho, eps, fill = tp2)) +
# 		geom_raster(interpolate = TRUE) + scale_fill_viridis('Exploitation efficiency',
# 								     option = 'C') +
# 		xlab(TeX('Proportion of LR'))+ylab(TeX('Sociality | LR + SR'))+ 
# 		theme(aspect.ratio = 1)
# 	, nrow = 1, ncol = 3, common.legend = TRUE, legend = 'top')
# 
# ggarrange(
# 	ggplot(data = median_tradeoff[.id == 'LR'],
# 	       aes(rho, eps, fill = tp2)) +
# 		geom_raster(interpolate = TRUE) +
# 		scale_fill_viridis('Exploration proportion', option = 'C') +
# 		xlab(TeX('Proportion of LR'))+ylab(TeX('Sociality | LR'))+ 
# 		theme(aspect.ratio = 1),
# 	ggplot(data = median_tradeoff[.id == 'SR'],
# 	       aes(rho, eps, fill = tp2)) +
# 		geom_raster(interpolate = TRUE) + scale_fill_viridis('Exploration proportion',
# 								     option = 'C') +
# 		xlab(TeX('Proportion of LR'))+ylab(TeX('Sociality | SR'))+ 
# 		theme(aspect.ratio = 1),
# 	ggplot(data = median_tradeoff[.id == 'both'],
# 	       aes(rho, eps, fill = tp2)) +
# 		geom_raster(interpolate = TRUE) + scale_fill_viridis('Exploration proportion',
# 								     option = 'C') +
# 		xlab(TeX('Proportion of LR'))+ylab(TeX('Sociality | LR + SR'))+ 
# 		theme(aspect.ratio = 1)
# 	, nrow = 1, ncol = 3, common.legend = TRUE, legend = 'top')
# 
# ggarrange(
# 	ggplot(data = median_total[.id == 'LR'],
# 	       aes(rho, eps, fill = tp2)) +
# 		geom_raster(interpolate = TRUE) +
# 		scale_fill_viridis('Total time', option = 'C') +
# 		xlab(TeX('Proportion of LR'))+ylab(TeX('Sociality | LR'))+ 
# 		theme(aspect.ratio = 1),
# 	ggplot(data = median_total[.id == 'SR'],
# 	       aes(rho, eps, fill = tp2)) +
# 		geom_raster(interpolate = TRUE) + scale_fill_viridis('Total time',
# 								     option = 'C') +
# 		xlab(TeX('Proportion of LR'))+ylab(TeX('Sociality | SR'))+ 
# 		theme(aspect.ratio = 1),
# 	ggplot(data = median_total[.id == 'both'],
# 	       aes(rho, eps, fill = tp2)) +
# 		geom_raster(interpolate = TRUE) + scale_fill_viridis('Total time',
# 								     option = 'C') +
# 		xlab(TeX('Proportion of LR'))+ylab(TeX('Sociality | LR + SR'))+ 
# 		theme(aspect.ratio = 1)
# 	, nrow = 1, ncol = 3, common.legend = TRUE, legend = 'top')
# 
# 
# median_tp2 <- rbindlist(list(LR = result_LR[is.finite(tp2), .(tp2 = 1/median(tp2)), by = c('rho', 'eps')],
# 			     SR = result_SR[is.finite(tp2), .(tp2 = 1/median(tp2)), by = c('rho', 'eps')],
# 			     both = result_both[is.finite(tp2), .(tp2 = 1/median(tp2)), by = c('rho', 'eps')]), idcol = TRUE)
# ggplot(data = median_tp2,
#        aes(rho, eps, fill = tp2)) +
# 	geom_raster(interpolate = TRUE) + scale_fill_viridis('Exploitation efficiency', breaks = seq(4e-4,8e-4, 2e-4),
# 							     option = 'C') +
# 	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+ 
# 	facet_wrap(~.id)+ theme(aspect.ratio = 1)
# 
# median_tradeoff <- rbindlist(list(LR = result_LR[is.finite(tp2), .(tp2 = median(tp1/(tp1 +tp2))), by = c('rho', 'eps')],
# 				  SR = result_SR[is.finite(tp2), .(tp2 = median(tp1/(tp1 +tp2))), by = c('rho', 'eps')],
# 				  both = result_both[is.finite(tp2), .(tp2 = median(tp1/(tp1 +tp2))), by = c('rho', 'eps')]), idcol = TRUE)
# ggplot(data = median_tradeoff,
#        aes(rho, eps, fill = tp2)) +
# 	geom_raster(interpolate = TRUE) + 
# 	scale_fill_viridis('Exploration-exploitation trade-off',
# 			   option = 'C') +
# 	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+ 
# 	facet_wrap(~.id)+ theme(aspect.ratio = 1)
# 
# median_total <- rbindlist(list(LR = result_LR[is.finite(tp2), .(tp2 = median((tp1 +tp2))/120), by = c('rho', 'eps')],
# 			       SR = result_SR[is.finite(tp2), .(tp2 = median((tp1 +tp2))/120), by = c('rho', 'eps')],
# 			       both = result_both[is.finite(tp2), .(tp2 = median((tp1 +tp2))/120), by = c('rho', 'eps')]), idcol = TRUE)
# ggplot(data = median_total,
#        aes(rho, eps, fill = tp2)) +
# 	geom_raster(interpolate = TRUE) + 
# 	scale_fill_viridis('Total time',
# 			   option = 'C') +
# 	xlab(TeX('Proportion of LR ($\\rho$)'))+ylab(TeX('Proportion of listeners ($\\epsilon$)'))+ 
# 	facet_wrap(~.id) + theme(aspect.ratio = 1)