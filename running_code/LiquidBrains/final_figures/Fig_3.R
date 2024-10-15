source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
source('~/research/gits/AnTracks/src/fit_functions.R')
load('~/research/gits/AnTracks/results/mov_rho_R.RData')

library(arrow)

load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/nf.RData')


path <- '/home/polfer/research/gits/AnTracks/results/sims/rho_075/'
files <- list.files(path)
files <- unique(unlist(regmatches(files,
				    gregexpr('rho_075_\\d{1,2}', files))))

R_sim <- seq(2.5, 20, 1.25)
R_exp <- R_sim * 50

mov_rho_R <- rbindlist(mov_rho_R)

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


## couple minutes to process
result_sim <- lapply(files, function(f){
	print(which(files == f))
	do.call('gc', args = list(verbose = FALSE))
	data <- setDT(read_parquet(paste0(path, f, '_data.parquet')))
	pos <- data[, .(node = unique(parse_nodes(pos))), by = 'Frame']
	
	
	mdata <- merge(pos, hex_sim[, c('x', 'y', 'node')])[order(Frame)]
	dists <- pdist(as.matrix(hex_sim[hex_sim$node == '(0, 22)', c('x', 'y')]), 
		       as.matrix(mdata[, c('x', 'y')]))
	indices <- vapply(R_sim, function(i){
		which.max(dists > i)
	}, numeric(1))
	times <- mdata[indices, Frame]
	data.table(R = R_exp, t = times - data[, min(Frame)], rho = 'sim')
})

estimated <- rbindlist(result_sim)

ggplot(data = rbind(rbindlist(list(det = det_result, nf = nf_result))[, .(t = median(t)),
								      by = c('R', 'rho')],
		    estimated[, .(t = median(t)), by = c('R', 'rho')]),
       aes(R/10, t / 120, fill = factor(rho, levels = c('det', 'nf', 'sim')),
           color = factor(rho, levels = c('det', 'nf', 'sim')),
           shape = factor(rho, levels = c('det', 'nf', 'sim'))))+
	
	
	
	scale_y_continuous('median first passage time (min)', breaks = seq(0, 30, 5))+
	scale_x_continuous('Radius (cm)', breaks = seq(0, 20, 2.5)*5)+
	scale_fill_manual('', values = c('mediumpurple', 'gold3', 'brown4'),
			  breaks = c('det', 'nf', 'sim'),
			  labels = c('DET', 'NFD', 'Model'))+
	scale_color_manual('', values = c('mediumpurple', 'gold3', 'brown4'),
			   labels = c('DET', 'NFD', 'Model'))+
	scale_linetype_manual('', values = c(1, 1, 2),
			      breaks = c('det', 'nf', 'sim'),
			      labels = c('DET', 'NFD', 'Model'))+
	scale_shape_manual('', values = c(22, 22, 24),
			   breaks = c('det', 'nf', 'sim'),
			   labels = c('DET', 'NFD', 'Model'))+
	# geom_path(data = data_fits_combined, 
	# 	  aes(linetype = factor(rho, levels = c('det', 'nf', 'LR', 'SR'))),
	# 	  linewidth = 1.75, alpha = 0.6, show.legend = FALSE)+
	geom_point(size = 5, color = 'black', alpha = 0.8)+
	
	theme(legend.position = c(0.2, 2/3),
	      legend.background = element_rect(fill = NA, color = 'black'),
	      legend.title = element_blank(), plot.title = element_text(size = 22),
	      aspect.ratio = 0.5) 


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




fit <- nls.multstart::nls_multstart(t ~ SSlogis(R, Asym, xmid, scal),
				      data = estimated[, .(t = mean(t)), by = 'R'],
				      supp_errors = 'Y', iter = 3e3,
				      start_lower = c(Asym = 0, xmid = 400, scal = 10**-5),
				      start_upper = c(Asym = 2500, xmid = 800, scal = 10**5))


data_fits <- rbind(data.table(R = R_exp_hd, t = SSlogis(R_exp_hd, coef(fit_DET)[['Asym']],
							coef(fit_DET)[['xmid']], coef(fit_DET)[['scal']]),
							rho = 'det'),
		   data.table(R = R_exp_hd, t = predict(fit_NFD,
		   				     newdata = data.frame(R = R_exp_hd)),
		   	   rho = 'nf'))

data_fits_combined <- rbind(data_fits, 
			    data.table(R = sort(unique(estimated[['R']])),
			    	   t = SSlogis(sort(unique(estimated[['R']])), 
			    	   	    coef(fit)[['Asym']],
			    	   	    coef(fit)[['xmid']],
			    	   	    coef(fit)[['scal']]),
			    	   rho = 'sim'))


panel_b <- ggplot(data = rbind(rbindlist(list(det = det_result, nf = nf_result))[, .(t = mean(t)),
								      by = c('R', 'rho')],
		    estimated[, .(t = mean(t)), by = c('R', 'rho')]),
       aes(R/10, t / 120, fill = factor(rho, levels = c('det', 'nf', 'sim')),
           color = factor(rho, levels = c('det', 'nf', 'sim')),
           shape = factor(rho, levels = c('det', 'nf', 'sim'))))+
	
	
	
	scale_y_continuous('Mean first passage time (min)', breaks = seq(0, 30, 5))+
	scale_x_continuous('Radius (cm)', breaks = seq(0, 20, 2.5)*5)+
	scale_fill_manual('', values = c('mediumpurple', 'gold3', 'brown4'),
			  breaks = c('det', 'nf', 'sim'),
			  labels = c('DET', 'NFD', 'Model'))+
	scale_color_manual('', values = c('mediumpurple', 'gold3','brown4'),
			   labels = c('DET', 'NFD', 'Model'))+
	scale_linetype_manual('', values = c(1, 1, 2, 2),
			      breaks = c('det', 'nf', 'sim'),
			      labels = c('DET', 'NFD', 'Model'))+
	scale_shape_manual('', values = c(22, 22, 24),
			   breaks = c('det', 'nf', 'sim'),
			   labels = c('DET', 'NFD', 'Model'))+
	
	geom_path(data = data_fits_combined, 
		  aes(linetype = factor(rho, levels = c('det', 'nf', 'sim'))),
		  linewidth = 1.75, alpha = 0.6, show.legend = FALSE)+
	geom_point(size = 5, color = 'black', alpha = 0.8)+
	
	theme(legend.position = c(0.2, 2/3),
	      legend.background = element_rect(fill = NA, color = 'black'),
	      legend.title = element_blank(), plot.title = element_text(size = 22),
	      aspect.ratio = 0.5) +
	ggtitle('B')


grob_inset <- ggplot(data = mov_rho_R[order(R, -scout),
				      .(N = mean(N, na.rm = TRUE)), by = c('rho', 'R')],
		     aes(R*5, N, color = rho, group = rho)) +
	geom_path(show.legend = FALSE)+
	geom_point(size = 3) + 
	ylab('Mean first\npassage (iters)') + 
	# scale_x_continuous('Radius (cm)', breaks = seq(0, 20, 5)*5, limits = c(0, 20)*5)+
	scale_x_continuous('Radius (cm)', breaks = seq(0, 20, 2.5)*5)+
	scale_color_viridis_c('Proportion of scouts')+
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


panel_D <- ggplot(data = mov_rho_R[order(R, -scout), .(N = mean(N), sd = sd(N)), 
				   by = c('scout', 'R')],
		  aes(R*5, N, color = factor(scout, levels = c('LR', 'SR'),
		  			   labels = c('Scouts', 'Recruits')),
		      shape = factor(scout, levels = c('LR', 'SR'),
		      	       labels = c('Scouts', 'Recruits')), 
		      fill = factor(scout, levels = c('LR', 'SR'),
		      	      labels = c('Scouts', 'Recruits')))) +
	geom_segment(aes(y = N - sd, yend = N + sd, xend = R*5), 
		     linewidth = 3, alpha = 0.75, show.legend = FALSE) +
	geom_path(linewidth = 1, show.legend = FALSE)+
	geom_point(size = 6, color = 'black') +
	ylab('Mean first passage (Iterations)') + 
	scale_x_continuous('Radius (cm)', breaks = seq(0, 20, 2.5)*5)+
	scale_color_viridis_d('Scout behaviour', direction = -1, end = 0.85)+
	scale_fill_viridis_d('Scout behaviour', direction = -1, end = 0.85)+
	scale_shape_manual('Scout behaviour', values = c(21, 24))+
	theme(legend.position = c(0.7, 4/5),
	      legend.background = element_rect(fill =NA, color = 'black'),
	      legend.key.size = unit(1, 'cm'), aspect.ratio = 0.5, 
	      plot.title = element_text(size = 22),
	      legend.margin = margin(b = 0, t = 5, l = 5, r = 5),
	      legend.title = element_blank())+
	annotation_custom(grob = ggplotGrob(grob_inset),
			  xmin = 8, xmax = 65, ymin = 200, ymax = 525)+
	ggtitle('A')

png('~/research/gits/AnTracks/plots/figures_LiquidBrains/Fig_3.png',
    2600*3, 3000*3, res = 300*3)
ggarrange(panel_D, panel_b, ncol = 1, nrow = 2)
dev.off()

# library(svglite)
# svglite::svglite('~/research/gits/AnTracks/plots/figures_LiquidBrains/Fig_3.svg',
# 		 8, 9)
# ggarrange(panel_D, panel_b, ncol = 1, nrow = 2)
# dev.off()
