source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')

load("/home/polfer/research/gits/AnTracks/results/results_sims_rho_eps.RData")
result_LR <- data_sims_rho_eps[['LR']]
result_SR <- data_sims_rho_eps[['SR']]

library(patchwork)


median_tp2 <- rbindlist(list(LR = enhance_land_res(result_LR[is.finite(tp2), .(tp2 = 1/median(tp2)),
							     by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
						   zvar = 'tp2', spar = 0.5, n = 200),
			     SR = enhance_land_res(result_SR[is.finite(tp2), .(tp2 = 1/median(tp2)),
			     				by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
			     				zvar = 'tp2', spar = 0.5, n = 200)), idcol = TRUE)


a <- ggplot(data = median_tp2,
	    aes(rho, eps, fill = tp2)) +
	geom_raster() + geom_contour(linetype = 4, aes(z = tp2), color = 'grey5') +
	scale_fill_viridis(TeX('Exploitation\nefficiency ($s^{-1}$)'), 
			   option = 'C', breaks = seq(0.0004, 0.0008, length.out = 5),
			   labels = sapply(seq(0.0004, 0.0008, length.out = 5), function(i){
			   	TeX(paste0(as.character(i*1000), '$\\cdot 10^{-3}$'))
			   })) +
	xlab(TeX('Proportion of scouts'))+ylab('')+ 
	theme(aspect.ratio = 1, strip.placement = 'outside', strip.background = element_blank(),
	      legend.position = 'right', legend.key.height = unit(1.5, 'cm'),
	      plot.title = element_text(size = 22))+
	facet_wrap(~ factor(.id, levels = c('LR', 'SR')),
		   strip.position = 'left', 
		   labeller = as_labeller(c(SR = 'Social copying | Recruits', 
		   			 LR = 'Social copying | Scouts')))+
	guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))+
	ggtitle('A')




median_total <- rbindlist(list(LR = enhance_land_res(result_LR[, .(total = 1/median((tp1 +tp2)/2)),
							       by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
						     zvar = 'total', spar = 0.5, n = 200),
			       SR = enhance_land_res(result_SR[, .(total = 1/median((tp1 +tp2)/2)),
			       				by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
			       				zvar = 'total', spar = 0.5, n = 200)), idcol = TRUE)

b <- ggplot(data = median_total,
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
	xlab(TeX('Proportion of scouts'))+ylab('')+ 
	theme(aspect.ratio = 1, strip.placement = 'outside', strip.background = element_blank(),
	      legend.position = 'right', legend.key.height = unit(1.5, 'cm'), 
	      plot.title = element_text(size = 22)
	)+
	facet_wrap(~ factor(.id, levels = c('LR', 'SR')),
		   strip.position = 'left', 
		   labeller = as_labeller(c(SR = 'Social copying | Recruits', 
		   			 LR = 'Social copying | Scouts')))+
	guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))+
	ggtitle('B')


a / b

png('/home/polfer/research/gits/AnTracks/plots/figures_LiquidBrains/Fig_suppl_heatmaps.png', 4400,4000, res = 460)
grid.arrange(a, b, nrow = 2, ncol = 1)
dev.off()

