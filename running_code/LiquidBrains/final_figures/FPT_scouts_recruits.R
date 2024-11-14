source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
load('~/research/gits/AnTracks/results/mov_rho_R.RData')

library(arrow)

# path <- '/home/polfer/research/gits/AnTracks/results/sims/movement_data/'
# files <- list.files(path)
# data_sims <- rbindlist(lapply(files, function(i){
# 	do.call('gc', args = list(verbose = FALSE))
# 	rho <- as.numeric(strsplit(i, split = '_')[[1]][2])
# 	print(rho)
# 	data <- data.table(read_parquet(paste0(path, i)))
# 	data[['t']] <- data[, .(t = seq_len(.N)), by = 'id'][['t']]
# 	pos <- data[, pos := sapply(pos, function(x) paste0(x, collapse = ","))]
# 	t <- pos[, .SD[1], by = c('pos', 'id')][['t']]
# 	data.table(t = t, rho = rho)
# }), idcol = 'exp')
# 
# save(data_sims, file = '~/research/gits/AnTracks/results/sims_visited_sites.RData')

load(file = '~/research/gits/AnTracks/results/sims_visited_sites.RData')

mov_rho_R <- rbindlist(mov_rho_R)

ggplot(data = data_sims[rho %in% c(0, 1) & t <= 1500], 
       aes(t, fill = factor(rho), color = factor(rho)))+
	geom_density(alpha = 0.5) 
	

ggplot(data = data_sims, aes(t, fill = factor(rho, levels = seq(1, 0, -0.1)),
			     color = factor(rho, levels = seq(1, 0, -0.1))))+
	geom_density(alpha = 0, linewidth = 1) + 
	coord_cartesian(xlim = c(0, 750))+
	scale_x_continuous('Time (iters)')+
	ylab('Visited sites (density)') +
	scale_fill_viridis_d('Proportion of scouts', direction = -1,
			     breaks = seq(0, 1, 0.1))+
	scale_color_viridis_d('Proportion of scouts', direction = -1,
			      breaks = seq(0, 1, 0.1))+
	theme(legend.position = c(0.66, 0.75), 
	      legend.background = element_rect(fill =NA, color = 'black', linewidth = 0.75),
	      legend.margin = margin(t = 5, l = 15, b = 5, r = 15),
	      aspect.ratio = 0.5, legend.direction = 'horizontal',
	      legend.title = element_text(size = 13), legend.key.width = unit(1, 'cm'),
	      legend.key.height = unit(0.2, 'cm'), legend.text = element_text(size = 11),
	      legend.spacing.x = unit(0, 'cm'),
	      axis.text = element_text(size = 11),
	      axis.title = element_text(size = 13))+
	guides(fill = guide_legend(override.aes = list(alpha = 1), nrow = 1,
				   title.position = 'top', title.hjust = 0.5,
				   label.position = 'bottom'))





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
	theme(legend.position = c(0.7, 0.85),
	      legend.key.size = unit(1, 'cm'), aspect.ratio = 0.5, 
	      plot.title = element_text(size = 22),
	      legend.margin = margin(b = 0, t = 5, l = 5, r = 5),
	      legend.title = element_blank())+
	annotation_custom(grob = ggplotGrob(grob_inset),
			  xmin = 8, xmax = 65, ymin = 200, ymax = 525)


















ldet <- lapply(det, function(i){
	x <- setDT(i@data)[, .SD[1], by = 'node'][order(Frame), Frame]
	y <- seq_along(x)
	data.table(t = x - i@data[, min(Frame)], visited_sites = y)
})




lnf <- lapply(nf[-c(1, 2)], function(i){
	x <- setDT(i@data)[, .SD[1], by = 'node'][order(Frame), Frame]
	y <- seq_along(x)
	data.table(t = x - i@data[, min(Frame)], visited_sites = y)
})


ggplot(data = rbindlist(list(det = rbindlist(ldet)[visited_sites < 414,
						   .(t = sum(t) / .N), by = 'visited_sites'],
			     
			     nf = rbindlist(lnf)[visited_sites < 414, .(t = sum(t) / .N), by = 'visited_sites']),
			idcol = 'exp'), aes(t / 120, visited_sites, color = exp))+
	geom_point() + geom_path()


ggplot(data = rbindlist(list(det = rbindlist(ldet),
			     
			     nf = rbindlist(lnf)),
			idcol = 'exp'), aes(t / 120, color = exp, fill = exp))+
	geom_density(linewidth = 1, alpha = 0.5, bw = 5)+
	scale_fill_manual('', values = c('mediumpurple', 'gold3'))+
	scale_color_manual('', values = c('mediumpurple', 'gold3'))+
	scale_x_continuous('Time (min)') + ylab('Density')+
	geom_vline(linetype = 2, linewidth = 1, 
		   xintercept = median(sapply(det, function(i){
		   	min(rbindlist(i@food)[['t']]) / 120
		   })))
