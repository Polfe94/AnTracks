source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
load('~/research/gits/AnTracks/results/mov_rho_R.RData')
source('~/research/gits/AnTracks/src/fit_functions.R')

library(arrow)

load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/nf.RData')

mov_rho_R <- rbindlist(mov_rho_R)

## ANALYSIS 1 : FPT ####

foodnodes <- unlist(lapply(det[1], function(i){
	p <- vapply(i@food[1], function(x){
		l <- pdist(as.matrix(x[, 1:2]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
		r <- range(l)
		avg <- mean(l)
		c(r, avg)
	}, numeric(3))
	data.table(mind = p[1], maxd = p[2], avgd = p[3])
}))

det_times <- unlist(lapply(det, function(i){
	f <- min(rbindlist(i@food)[['t']])
	f - setDT(i@data)[1, Frame]
}))
det_times <- c(range(det_times), mean(det_times))



R_exp <- seq(2.5, 20, 1.25) * 50


det_result <- rbindlist(lapply(det, function(p){
	data <- setDT(p@data)
	mdata <- merge(data[, c('node', 'Frame')], hex[, c('x', 'y', 'node')])[order(Frame)]
	dists <- pdist(as.matrix(hex[hex$node == 634, c('x', 'y')]), as.matrix(mdata[, c('x', 'y')]))
	indices <- vapply(R_exp, function(i){
		which.max(dists > i)
	}, numeric(1))
	times <- mdata[indices, Frame]
	data.table(R = R_exp, t = times - data[, min(Frame)], rho = 'det')
}))

nf_result <- rbindlist(lapply(nf[-c(1, 2)], function(p){
	data <- setDT(p@data)
	mdata <- merge(data[, c('node', 'Frame')], hex[, c('x', 'y', 'node')])[order(Frame)]
	dists <- pdist(as.matrix(hex[hex$node == 634, c('x', 'y')]), as.matrix(mdata[, c('x', 'y')]))
	indices <- vapply(R_exp, function(i){
		which.max(dists > i)
	}, numeric(1))
	times <- mdata[indices, Frame]
	data.table(R = R_exp, t = times - data[, min(Frame)], rho = 'nf')
}))



panel_a <- ggplot(data = rbindlist(list(det = det_result, nf = nf_result))[R %in% seq(2.5*50, 20*50, 50*2.5)],
		  aes(factor(R/10), t / 120, fill = factor(rho)))+
	geom_boxplot(outlier.shape = NA, show.legend = FALSE, alpha = 0.6)+
	geom_jitter(width = 0.15, alpha = 0.4, size = 3, show.legend = FALSE)+
	facet_wrap(~ factor(rho, labels = c('DET', 'NFD')))+
	scale_y_continuous('First passage time (min)', breaks = seq(0, 80, 15))+
	scale_x_discrete('Radius (cm)', breaks = seq(0, 20, 2.5)*5)+
	scale_fill_manual('', values = c('mediumpurple', 'gold3'))+ ggtitle('A')+
	theme(plot.title = element_text(size = 22))



# fit_DET <- nls.multstart::nls_multstart(t ~ SSlogis(R, Asym, xmid, scal),
# 					data = det_result[, .(t = median(t)), by = 'R'],
# 					supp_errors = 'Y', iter = 3e3,
# 					start_lower = c(Asym = 0, xmid = 400, scal = 10**-5),
# 					start_upper = c(Asym = 2500, xmid = 800, scal = 10**5))
# fit_NFD <- lm(t ~ poly(R,2), data = nf_result[, .(t = median(t)), by = 'R'])
fit_DET <- nls.multstart::nls_multstart(t ~ SSlogis(R, Asym, xmid, scal),
					data = det_result[, .(t = mean(t)), by = 'R'],
					supp_errors = 'Y', iter = 3e3,
					start_lower = c(Asym = 0, xmid = 400, scal = 10**-5),
					start_upper = c(Asym = 2500, xmid = 800, scal = 10**5))
fit_NFD <- lm(t ~ poly(R,2), data = nf_result[, .(t = mean(t)), by = 'R'])


R_exp_hd <- seq(2.5, 20, length.out = 150)*50

data_fits <- rbind(data.table(R = R_exp_hd, t = SSlogis(R_exp_hd, coef(fit_DET)[['Asym']],
							coef(fit_DET)[['xmid']], coef(fit_DET)[['scal']]),
							rho = 'det'),
		   data.table(R = R_exp_hd, t = predict(fit_NFD,
		   				     newdata = data.frame(R = R_exp_hd)),
		   	   rho = 'nf'))


panel_b <- ggplot(data = rbindlist(list(det = det_result, nf = nf_result))[, .(t = mean(t)), # .(t = median(t))
									   by = c('R', 'rho')],
		  aes(R/10, t / 120, fill = factor(rho), color = factor(rho)))+
	geom_polygon(fill = 'grey90', color = NA,
		     data = data.frame(x = foodnodes[c(1, 1, 2, 2)]/10,
		     		  y = c(-Inf, Inf, Inf, -Inf)),
		     aes(x = x, y = y)) +
	geom_polygon(fill = 'grey90', color = NA,
		     data = data.frame(x = c(-Inf, Inf, Inf, -Inf),
		     		  y = det_times[c(1, 1, 2, 2)]/120),
		     aes(x = x, y = y)) +
	geom_vline(xintercept = foodnodes[3]/10, linewidth = 1, linetype = 2, color = 'grey35')+
	geom_hline(yintercept = det_times[3]/120, linewidth = 1, linetype = 2, color = 'grey35')+
	scale_y_continuous('Median first passage time (min)', breaks = seq(0, 80, 10))+
	scale_x_continuous('Radius (cm)', breaks = seq(0, 20, 2.5)*5)+
	scale_fill_manual('', values = c('mediumpurple', 'gold3'),
			  labels = c('DET', 'NFD'))+
	scale_color_manual('', values = c('mediumpurple', 'gold3'),
			   labels = c('DET', 'NFD'))+
	geom_path(data = data_fits,
		  size = 1.75, alpha = 0.6, show.legend = FALSE)+
	geom_point(size = 5, shape = 21, color = 'black', alpha = 0.6)+
	theme(legend.position = c(0.2, 2/3),
	      legend.background = element_rect(fill = NA, color = 'black'),
	      legend.title = element_blank(), plot.title = element_text(size = 22)) +
	ggtitle('B')


grob_inset <- ggplot(data = mov_rho_R[,
				      .(N = mean(N, na.rm = TRUE)), by = c('rho', 'R')],
		     aes(R*5, N, color = rho, group = rho)) +
	geom_point(size = 3) + geom_path(show.legend = FALSE)+
	ylab('Mean first\npassage (iters)') + 
	# scale_x_continuous('Radius (cm)', breaks = seq(0, 20, 5)*5, limits = c(0, 20)*5)+
	scale_x_continuous('Radius (cm)', breaks = seq(0, 20, 2.5)*5)+
	scale_color_viridis_c('Proportion of LR')+
	theme(legend.position = c(0.33, 0.75), 
	      legend.background = element_rect(fill =NA, color = 'black', linewidth = 0.75),
	      plot.background = element_rect(fill =NA, color = 'black'),
	      legend.margin = margin(t = 5, l = 15, b = 5, r = 15),
	      aspect.ratio = 0.5, legend.direction = 'horizontal',
	      legend.title = element_text(size = 13), legend.key.width = unit(0.75, 'cm'),
	      legend.key.height = unit(0.2, 'cm'), legend.text = element_text(size = 11),
	      axis.text = element_text(size = 11),
	      axis.title = element_text(size = 13))+
	guides(color = guide_colorbar(title.position = 'top', title.hjust = 0.5))


panel_D <- ggplot(data = mov_rho_R[, .(N = mean(N), sd = sd(N)), by = c('scout', 'R')],
		  aes(R*5, N, color = factor(scout, levels = c('LR', 'SR')),
		      shape = factor(scout, levels = c('LR', 'SR')), 
		      fill = factor(scout, levels = c('LR', 'SR')))) +
	geom_segment(aes(y = N - sd, yend = N + sd, xend = R*5), 
		     linewidth = 3, alpha = 0.75, show.legend = FALSE) +
	# geom_pointrange(aes(ymin = N - sd, ymax = N + sd), size = 2, alpha = 0.75,
	# 		linewidth = 3) +
	geom_path(linewidth = 1, show.legend = FALSE)+
	geom_point(size = 6, color = 'black') +
	ylab('Mean first passage (Iterations)') + 
	# scale_x_continuous('Radius (cm)', breaks = seq(0, 20, 5)*5, limits = c(0, 20)*5)+
	scale_x_continuous('Radius (cm)', breaks = seq(0, 20, 2.5)*5)+
	scale_color_viridis_d('Scout behaviour', direction = -1, end = 0.85)+
	scale_fill_viridis_d('Scout behaviour', direction = -1, end = 0.85)+
	scale_shape_manual('Scout behaviour', values = c(21, 24))+
	theme(legend.position = c(0.7, 4/5),
	      legend.background = element_rect(fill =NA, color = 'black'),
	      legend.key.size = unit(1, 'cm'), aspect.ratio = 0.5, 
	      plot.title = element_text(size = 22),
	      legend.margin = margin(b = 0, t = 5, l = 5, r = 5))+
	# theme(legend.position = c(0.2, 2/3), 
	#       legend.background = element_rect(fill =NA, color = 'black'),
	#       legend.key.size = unit(1, 'cm'), aspect.ratio = 0.5)+
	# guides(color = guide_legend(override.aes = list(alpha = 1, linetype = 0)))+
	annotation_custom(grob = ggplotGrob(grob_inset),
			  xmin = 8, xmax = 65, ymin = 200, ymax = 525)+
	ggtitle('A')



estimated <- mov_rho_R[rho == 0 | rho == 1, .(t = mean(N)), by = c('R', 'rho')]
estimated[['rho']] <- ifelse(estimated[['rho']] == 0, 'SR', 'LR')
estimated[['t']] <- estimated[['t']] * (50/10) * 2
estimated[['R']] <- estimated[['R']] *50



fit_LR <- nls.multstart::nls_multstart(t ~ fexpC(R, a, k, C),
					data = estimated[rho == 'LR'],
					supp_errors = 'Y', iter = 3e3,
					start_lower = c(a = 50, k = 0, C = -10),
					start_upper = c(a = 100, k = 1, C = 5))
fit_SR <- nls.multstart::nls_multstart(t ~ fexpC(R, a, k, C),
				       data = estimated[rho == 'SR'],
				       supp_errors = 'Y', iter = 3e3,
				       start_lower = c(a = 50, k = 0, C = -10),
				       start_upper = c(a = 100, k = 1, C = 5))

data_fits_combined <- rbind(data_fits, data.table(R = sort(unique(estimated[['R']])),
			      t = fexpC(sort(unique(estimated[['R']])), 
			      	 coef(fit_LR)[['a']], coef(fit_LR)[['k']],
						  coef(fit_LR)[['C']]),
			      rho = 'LR'),
		   data.table(R = sort(unique(estimated[['R']])),
		   	   t = fexpC(sort(unique(estimated[['R']])), 
		   	   	 coef(fit_SR)[['a']], coef(fit_SR)[['k']],
		   	   	 coef(fit_SR)[['C']]), rho = 'SR'))

panel_b <- ggplot(data = rbind(rbindlist(list(det = det_result, nf = nf_result))[, .(t = mean(t)),
									   by = c('R', 'rho')],
			       estimated)[order(rho)],
		  aes(R/10, t / 120, fill = factor(rho, levels = c('det', 'nf', 'LR', 'SR')),
		      color = factor(rho, levels = c('det', 'nf', 'LR', 'SR')),
		      shape = factor(rho, levels = c('det', 'nf', 'LR', 'SR'))))+
	scale_y_continuous('Mean first passage time (min)', breaks = seq(0, 30, 5))+
	scale_x_continuous('Radius (cm)', breaks = seq(0, 20, 2.5)*5)+
	scale_fill_manual('', values = c('mediumpurple', 'gold3', 'brown4', 'brown4'),
			  breaks = c('det', 'nf', 'LR', 'SR'),
			  labels = c('DET', 'NFD', 'Model (LR)', 'Model (SR)'))+
	scale_color_manual('', values = c('mediumpurple', 'gold3', 'brown4', 'brown4'),
			   labels = c('DET', 'NFD', 'Model (LR)', 'Model (SR)'))+
	scale_linetype_manual('', values = c(1, 1, 2, 2),
			      breaks = c('det', 'nf', 'LR', 'SR'),
			      labels = c('DET', 'NFD', 'Model (LR)', 'Model (SR)'))+
	scale_shape_manual('', values = c(22, 22, 21, 24),
			   breaks = c('det', 'nf', 'LR', 'SR'),
			   labels = c('DET', 'NFD', 'Model (LR)', 'Model (SR)'))+
	geom_path(data = data_fits_combined, 
		  aes(linetype = factor(rho, levels = c('det', 'nf', 'LR', 'SR'))),
		  linewidth = 1.75, alpha = 0.6, show.legend = FALSE)+
	geom_point(size = 5, color = 'black', alpha = 0.8)+
	annotation_custom(grid::linesGrob(), xmin = -Inf, xmax = 8.5, 
			  ymax = det_times[1]/120, ymin = det_times[1]/120)+
	annotation_custom(grid::linesGrob(), xmin = -Inf, xmax = 8.5, 
			  ymax = det_times[2]/120, ymin = det_times[2]/120)+
	annotation_custom(grid::linesGrob(), xmin = -Inf, xmax = 8.5, 
			  ymax = det_times[3]/120, ymin = det_times[3]/120)+
	
	annotation_custom(grid::linesGrob(), ymin = -Inf, ymax = -1.2, 
			  xmax = foodnodes[['avgd']]/10, foodnodes[['avgd']]/10)+
	# annotation_custom(grid::textGrob(label = 'Average\nfood\ndistance'), xmin = 66.5,
	# 		  xmax = 68, ymin = -0.75, ymax = 1.2)+
	# annotation_custom(grid::textGrob(label = c(round(det_times[1]/120, 1))), xmin = 8.75,
	# 		  xmax = 9.25, ymin = det_times[1]/120, ymax = det_times[1]/120)+
	# annotation_custom(grid::textGrob(label = c(round(det_times[2]/120, 1))), xmin = 8.75,
	# 		  xmax = 9.25, ymin = det_times[2]/120, ymax = det_times[2]/120)+
	# annotation_custom(grid::textGrob(label = c(round(det_times[3]/120, 1))), xmin = 8.75,
	# 		  xmax = 9.25, ymin = det_times[3]/120, ymax = det_times[3]/120)+
	# geom_segment(data = data.frame(y = det_times/120, x = 11, rho = 'LR'),
	# 	     aes(x, y, xend = 12.5, yend = y), color = 'black', linetype = 1)+
	theme(legend.position = c(0.2, 2/3),
	      legend.background = element_rect(fill = NA, color = 'black'),
	      legend.title = element_blank(), plot.title = element_text(size = 22),
	      aspect.ratio = 0.5) +
	ggtitle('B')


# png('~/research/gits/AnTracks/plots/figures_LiquidBrains/Fig_3.png', 8000, 8000, res = 825)
# grid.arrange(panel_D, panel_b)
# dev.off()
png('~/research/gits/AnTracks/plots/figures_LiquidBrains/Fig_3.png', 7000, 8000, res = 800)
grid.arrange(panel_D, panel_b)
dev.off()

########### WITH MEDIANS ###########
# 
# fit_DET <- nls.multstart::nls_multstart(t ~ SSlogis(R, Asym, xmid, scal),
# 					data = det_result[, .(t = median(t)), by = 'R'],
# 					supp_errors = 'Y', iter = 3e3,
# 					start_lower = c(Asym = 0, xmid = 400, scal = 10**-5),
# 					start_upper = c(Asym = 2500, xmid = 800, scal = 10**5))
# fit_NFD <- lm(t ~ poly(R,2), data = nf_result[, .(t = median(t)), by = 'R'])
# 
# R_exp_hd <- seq(2.5, 20, length.out = 150)*50
# 
# data_fits <- rbind(data.table(R = R_exp_hd, t = SSlogis(R_exp_hd, coef(fit_DET)[['Asym']],
# 							coef(fit_DET)[['xmid']], coef(fit_DET)[['scal']]),
# 							rho = 'det'),
# 		   data.table(R = R_exp_hd, t = predict(fit_NFD,
# 		   				     newdata = data.frame(R = R_exp_hd)),
# 		   	   rho = 'nf'))
# 
# 
# estimated <- mov_rho_R[rho == 0 | rho == 1, .(t = median(N)), by = c('R', 'rho')]
# estimated[['rho']] <- ifelse(estimated[['rho']] == 0, 'SR', 'LR')
# estimated[['t']] <- estimated[['t']] * (50/12.5) * 2
# estimated[['R']] <- estimated[['R']] *50
# 
# panel_b <- ggplot(data = rbind(rbindlist(list(det = det_result, nf = nf_result))[, .(t = median(t)),
# 										 by = c('R', 'rho')],
# 			       estimated)[order(rho)],
# 		  aes(R/10, t / 120, fill = factor(rho, levels = c('det', 'nf', 'LR', 'SR')),
# 		      color = factor(rho, levels = c('det', 'nf', 'LR', 'SR')),
# 		      shape = factor(rho, levels = c('det', 'nf', 'LR', 'SR'))))+
# 	scale_y_continuous('Median first passage time (min)', breaks = seq(0, 80, 10))+
# 	scale_x_continuous('Radius (cm)', breaks = seq(0, 20, 2.5)*5)+
# 	scale_fill_manual('', values = c('mediumpurple', 'gold3', 'brown4', 'brown4'),
# 			  breaks = c('det', 'nf', 'LR', 'SR'),
# 			  labels = c('DET', 'NFD', 'Model (LR)', 'Model (SR)'))+
# 	scale_color_manual('', values = c('mediumpurple', 'gold3', 'brown4', 'brown4'),
# 			   labels = c('DET', 'NFD', 'Model (LR)', 'Model (SR)'))+
# 	
# 	scale_shape_manual('', values = c(22, 22, 21, 24),
# 			   breaks = c('det', 'nf', 'LR', 'SR'),
# 			   labels = c('DET', 'NFD', 'Model (LR)', 'Model (SR)'))+
# 	geom_path(data = data_fits,
# 		  size = 1.75, alpha = 0.6, show.legend = FALSE)+
# 	geom_point(size = 5, color = 'black', alpha = 0.6)+
# 	theme(legend.position = c(0.2, 2/3),
# 	      legend.background = element_rect(fill = NA, color = 'black'),
# 	      legend.title = element_blank(), plot.title = element_text(size = 22)) +
# 	ggtitle('B')
