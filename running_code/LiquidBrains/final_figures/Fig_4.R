source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
load('~/research/gits/AnTracks/data/det.RData')

load("/home/polfer/research/gits/AnTracks/results/results_sims_rho_eps.RData")
result_both <- data_sims_rho_eps[['both']]

## Enhance landscape resolution
median_total <- enhance_land_res(result_both[, .(z = 1/median((tp1 +tp2))),
					     by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
				 zvar = 'z', spar = 0.5, n = 200)
median_tp1 <- enhance_land_res(result_both[is.finite(tp1), .(z = 1/median(tp1)),
					   by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
					   zvar = 'z', spar = 0.5, n = 200)
median_tp2 <- enhance_land_res(result_both[is.finite(tp2), .(z = 1/median(tp2)),
					   by = c('rho', 'eps')], xvar = 'rho', yvar = 'eps', 
					   zvar = 'z', spar = 0.5, n = 200)

## Empirical data
det_tps <- rbindlist(lapply(det, function(i){
	t <- range(rbindlist(i@food)[['t']])
	data.frame(tp1 = t[1], tp2 = t[2] - t[1]) / 2
}), idcol = 'exp')


## Map empirical Z contours
cs_tp1 <- contourLines(z = matrix(median_tp1$z, ncol = 200), 
		   levels = 1/median(det_tps$tp1))
poly_tp1 <- rbindlist(lapply(cs_tp1, function(i){
	data.table(x = i$y, y = i$x, z = i$level
	)}))

cs_tp2 <- contourLines(z = matrix(median_tp2$z, ncol = 200), 
		       levels = 1/(median(det_tps$tp2)))
poly_tp2 <- rbindlist(lapply(cs_tp2, function(i){
	data.table(x = i$y, y = i$x, z = i$level
	)}))

cs_total <- contourLines(z = matrix(median_total$z, ncol = 200), 
		       levels = 1/(median(det_tps$tp2 + det_tps$tp1)))
poly_total <- rbindlist(lapply(cs_total, function(i){
	data.table(x = i$y, y = i$x, z = i$level
	)}))

multiplot <- ggarrange(
ggplot(data = median_tp1,
       aes(rho, eps, fill = z)) +
	geom_raster() + geom_contour(linetype = 4, aes(z = z), color = 'grey5') +
	scale_fill_viridis(TeX('Exploration efficiency ($10^{3} s^{-1}$)'), 
			   option = 'C',
			   breaks = seq(0.001, 0.0015, length.out = 5),
			   labels = sapply(seq(0.001, 0.0015, length.out = 5), 
			   		function(i){
			   			TeX(paste0(as.character(i*1000)))
			   		}))+
	xlab('Proportion of LR') + ylab ('Social copying')+
	ggtitle('A')+ theme(aspect.ratio = 1, legend.position = 'top', 
			    legend.key.width = unit(1.4, 'cm'))+
	guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))+
	geom_point(data = poly_tp1, 
		   aes(x, y), color = 'white', fill = NA, size = 0.3),

ggplot(data = median_tp2,
       aes(rho, eps, fill = z)) +
	geom_raster() + geom_contour(linetype = 4, aes(z = z), color = 'grey5') +
	scale_fill_viridis(TeX('Exploitation efficiency ($10^{3} s^{-1}$)'), 
			   option = 'C',
			   breaks = seq(0.00045, 0.00065, length.out = 5),
			   labels = sapply(seq(0.00045, 0.00065, length.out = 5), 
			   		function(i){
			   			TeX(paste0(as.character(i*1000)))
			   		}))+
	xlab('Proportion of LR') + ylab ('Social copying')+
	ggtitle('B') + theme(aspect.ratio = 1, legend.position = 'top', 
			     legend.key.width = unit(1.4, 'cm'))+
	guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))+ 
	geom_point(data = poly_tp2, 
		   aes(x, y), color = 'white', fill = NA, size = 0.3),
ggplot(data = median_total,
       aes(rho, eps, fill = z)) +
	geom_raster() + geom_contour(linetype = 4, aes(z = z), color = 'grey5') +
	scale_fill_viridis(TeX('Foraging efficiency ($10^{3} s^{-1}$)'), 
			   option = 'C',
			   breaks = seq(0.0003, 0.00045, length.out = 4),
			   labels = sapply(seq(0.0003, 0.00045, length.out = 4), 
			   		function(i){
			   			TeX(paste0(as.character(i*1000)))
			   		}))+
	xlab('Proportion of LR') + ylab ('Social copying')+
	ggtitle('C')+ theme(aspect.ratio = 1, legend.position = 'top', 
			    legend.key.width = unit(1.4, 'cm'))+
	guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))+ 
	geom_point(data = poly_total, 
		   aes(x, y), color = 'white', fill = NA, size = 0.3),
ncol = 3, nrow = 1)


png('/home/polfer/research/gits/AnTracks/plots/figures_LiquidBrains/Fig_4.png', 6000,2750, res = 500)
multiplot
dev.off()

