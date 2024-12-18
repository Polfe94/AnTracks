source('~/research/gits/AnTracks/src/Experiment.R')
load('~/research/gits/AnTracks/data/nf.RData')
load('~/research/gits/AnTracks/data/det.RData')

library(patchwork)

det <- lapply(det, class_ids)
nf <- lapply(nf[-c(1, 2)], class_ids)

det_stats <- rbindlist(lapply(det, function(i){
	x <- i
	inds <- x@data[type != 'UNKNOWN', .N, by = 'N_ind'][N > 10, N_ind]
	x@data <- x@data[N_ind %in% inds]
	move_stats(x)
}))

nf_stats <- rbindlist(lapply(nf, function(i){
	x <- i
	inds <- x@data[type != 'UNKNOWN', .N, by = 'N_ind'][N > 10, N_ind]
	x@data <- x@data[N_ind %in% inds]
	move_stats(x)
}))

stats <- rbindlist(list(det = det_stats, nf = nf_stats), idcol = 'exp')
stats <- na.omit(stats)


stats_mlt <- data.table(melt(stats[, c('exp', 'v', 'acc', 'd', 't', 'type')],
		  id.vars = c('exp', 'type')))
stats_mlt[['interaction']] <- apply(stats_mlt[, c('exp', 'type')], 1,
				    paste0, collapse = '', sep='')

## transformation of units 
stats_mlt[variable == 'v', value := value / 10] # to cm/s
stats_mlt[variable == 'd', value := value / 10] # to cm
stats_mlt[variable == 't', value := value / 60] # to min

pairwise.wilcox.test(stats_mlt[variable == 'v', value], 
		     stats_mlt[variable == 'v', interaction], 
		     p.adjust.method = 'bonferroni')

subpanel_1 <- ggplot(data = stats_mlt[variable == 'v'],
		     aes(factor(type, levels = c('Scout', 'Recruit')),
		         value, fill = exp)) + 
	
	geom_boxplot(alpha = 0.6, outlier.shape = 1)+
	annotate('text', x = 0.81 - 0.15, y = 2.25, label = 'a', size = 10) +
	annotate('text', x = 1.19 - 0.15, y = 2.25, label = 'a', size = 10) +
	annotate('text', x = 1.81 - 0.15, y = 2.25, label = 'b', size = 10) +
	annotate('text', x = 2.19 - 0.15, y = 2.25, label = 'c', size = 10) +
	scale_fill_manual('', labels = c('Food', 'No-Food'), 
			  values = c('mediumpurple','gold3'))+
	xlab('') +
	scale_y_continuous(TeX("Velocity ($cm\\cdot s^{-1}$)"),
			   limits = c(0, NA))+
	theme(legend.position = c(0.4, 0.8), # legend.position = c(0.15, 0.8), 
	      legend.background = element_rect(color = 'white', fill = 'white'), 
	      legend.justification = 'center', legend.title = element_blank(),
	      legend.text = element_text(size = 22),
	      legend.key.size = unit(1, 'cm'),
	      plot.title = element_text(size = 22))

pairwise.wilcox.test(stats_mlt[variable == 'd', value], 
		     stats_mlt[variable == 'd', interaction], 
		     p.adjust.method = 'bonferroni')

subpanel_3 <- ggplot(data = stats_mlt[variable == 'd'],
		     aes(factor(type, levels = c('Scout', 'Recruit')),
		         value, fill = exp)) + 
	
	geom_boxplot(alpha = 0.6, outlier.shape = 1, show.legend = FALSE)+
	annotate('text', x = 0.81 - 0.15, y = 140, label = 'a', size = 10) +
	annotate('text', x = 1.19 - 0.15, y = 140, label = 'b', size = 10) +
	annotate('text', x = 1.81 - 0.15, y = 140, label = 'c', size = 10) +
	annotate('text', x = 2.19 - 0.15, y = 140, label = 'c', size = 10) +
	scale_fill_manual('', labels = c('Experimental', 'Control'), 
			  values = c('mediumpurple','gold3'))+
	xlab('') +
	scale_y_continuous(TeX('Maximum distance (cm)'))+
	theme(legend.position = c(0.15, 0.8), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.justification = 'center', legend.title = element_blank(),
	      plot.title = element_text(size = 22))+
	coord_cartesian(ylim = c(0, 150))


pairwise.wilcox.test(stats_mlt[variable == 't', value], 
		     stats_mlt[variable == 't', interaction], 
		     p.adjust.method = 'bonferroni')

subpanel_4 <- ggplot(data = stats_mlt[variable == 't'],
		     aes(factor(type, levels = c('Scout', 'Recruit')),
		         value, fill = exp)) + 
	
	geom_boxplot(alpha = 0.6, outlier.shape = 1, show.legend = FALSE)+
	annotate('text', x = 0.81 - 0.15, y = 18, label = 'a', size = 10) +
	annotate('text', x = 1.19 - 0.15, y = 18, label = 'b', size = 10) +
	annotate('text', x = 1.81 - 0.15, y = 18, label = 'c', size = 10) +
	annotate('text', x = 2.19 - 0.15, y = 18, label = 'c', size = 10) +
	scale_fill_manual('', labels = c('Experimental', 'Control'), 
			  values = c('mediumpurple','gold3'))+
	xlab('') +
	scale_y_continuous(TeX('Time (min)'))+
	theme(legend.position = c(0.15, 0.8), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.justification = 'center', legend.title = element_blank(),
	      plot.title = element_text(size = 22))+
	coord_cartesian(ylim = c(0, 20))

pls_det <- rbindlist(lapply(det, function(i){
	data <- i@data[type != 'UNKNOWN']
	rbindlist(lapply(unique(data[['N_ind']]), function(i){
		a <- data[N_ind == i][c(1, diff(node)) != 0, node]
		m <- as.matrix(hex[a, c('x', 'y')])
		if(nrow(m) > 3){
			dirs <- get_directions_cpp(m)
			data.frame(pl = mean(find_pl(dirs)), type = data[N_ind == i, unique(type)])
		} else {
			data.frame()
		}
	}), idcol = 'N_ind')
}), idcol = 'exp')

pls_det[!is.finite(pl), 'pl'] <- 0

pls <- rbindlist(lapply(nf, function(i){
	data <- i@data[type != 'UNKNOWN']
	rbindlist(lapply(unique(data[['N_ind']]), function(i){
		a <- data[N_ind == i][c(1, diff(node)) != 0, node]
		m <- as.matrix(hex[a, c('x', 'y')])
		if(nrow(m) > 3){
			dirs <- get_directions_cpp(m)
			data.frame(pl = mean(find_pl(dirs)), type = data[N_ind == i, unique(type)])
		} else {
			data.frame()
		}
	}), idcol = 'N_ind')
}), idcol = 'exp')

pls[!is.finite(pl), 'pl'] <- 0

data_pl <- rbindlist(list(Experimental = pls_det, Control = pls), idcol = 'condition')

subpanel_5 <- ggplot(data = data_pl,
       aes(factor(type, levels = c('Scout', 'Recruit')),
           pl, fill = condition)) + 
	
	geom_boxplot(alpha = 0.6, outlier.shape = 1, show.legend = FALSE)+
	annotate('text', x = 0.81 - 0.15, y = 6, label = 'a', size = 10) +
	annotate('text', x = 1.19 - 0.15, y = 6, label = 'b', size = 10) +
	annotate('text', x = 1.81 - 0.15, y = 6, label = 'c', size = 10) +
	annotate('text', x = 2.19 - 0.15, y = 6, label = 'd', size = 10) +
	scale_fill_manual('', labels = c('Experimental', 'Control'), 
			  values = c('mediumpurple','gold3'))+
	xlab('') +
	scale_y_continuous(TeX('Directional persistence (steps)'))+
	coord_cartesian(ylim = c(0, 7))



panel_a <- (subpanel_1 + theme(axis.line.x = element_blank(),
			      axis.text.x = element_blank(),
			      axis.ticks.x = element_blank()) +
subpanel_5+ theme(axis.line.x = element_blank(),
		  axis.text.x = element_blank(),
		  axis.ticks.x = element_blank())) /
(subpanel_3+ subpanel_4)


draw_traffic_flow <- function(edges = NULL, add = NULL, ...){
	
	if(is.null(edges)){
		edges <- compute_edges(hex[hex$y > 1000, ])
	}
	
	if(!is.null(add)){
		pl <- add
	} else {
		pl <- ggplot()
	}
	
	pl <- pl + geom_segment(data = edges, 
				aes(x = x, xend = xend, y = y,
				    yend = yend), color = 'black',size = 2.9,
				show.legend = FALSE, ...) +
		geom_segment(data = edges, 
			     aes(x = x, xend = xend, y = y,
			         yend = yend, color = z), size = 2.5, ...)
	
	pl + xlab('') + ylab('') + theme(axis.text = element_blank(), axis.ticks = element_blank(),
					 axis.line = element_blank())
}


det_posits <- rbindlist(lapply(det, function(i){
	data <- i@data[type != 'UNKNOWN']
	data[['d']] <- get_segment(data[, c('Xmm', 'Ymm')])[, 2]
	colnames(data)[colnames(data) == 'node'] <- 'o'
	data
}), idcol = 'exp')
det_posits[['id']] <- apply(det_posits[, c('exp', 'N_ind')], 1, paste0, sep = '', collapse = '_')
det_density <- det_posits[, .(N = .N), by = c('type','o','d')][, z := rank(N)]
det_density[['condition']] <- 'DET'

nf_posits <- rbindlist(lapply(nf, function(i){
	do.call('gc', args = list(verbose = FALSE))
	data <- i@data[type != 'UNKNOWN']
	data[['d']] <- get_segment(data[, c('Xmm', 'Ymm')])[, 2]
	colnames(data)[colnames(data) == 'node'] <- 'o'
	data
}), idcol = 'exp')
nf_posits[['id']] <- apply(nf_posits[, c('exp', 'N_ind')], 1, paste0, sep = '', collapse = '_')
nf_density <- nf_posits[, .(N = .N), by = c('type','o','d')][, z := rank(N)]
nf_density[['condition']] <- 'NFD'

traffic_flow <- merge(rbind(det_density, nf_density), edges)


panel_c_1 <- draw_traffic_flow(traffic_flow[type == 'Scout' & condition == 'DET'], lineend = 'round',
			       add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round',
			       		    add = geom_foodpatches(fill = 'mediumpurple')))+
	scale_size_continuous(range = c(1, 5)) + 
	scico::scale_color_scico() + theme(legend.position = 'none',
					   strip.text.y = element_text(size = 15),
					   aspect.ratio = 0.5, plot.title = element_text(size =22))+
	geom_circle(hex[hex$node == 634, c('x', 'y')], 
		    r = 200, npoints = 500, linetype = 2, linewidth = 0.95) +
	scale_y_continuous(limits = c(1000, 1950)) + 
	geom_point(shape = 23, data = hex[hex$node == 634, c('x', 'y')], aes(x, y-10),
		   fill = 'purple3', size = 4)+
	ggtitle('Scouts (Food)') + theme(plot.title = element_text(hjust = 0.5, vjust = -4))

panel_c_2 <- draw_traffic_flow(traffic_flow[type == 'Recruit' & condition == 'DET'], lineend = 'round',
			       add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round',
			       		    add = geom_foodpatches(fill = 'mediumpurple')))+
	scale_size_continuous(range = c(1, 5)) + 
	scico::scale_color_scico() + theme(legend.position = 'none',
					   strip.text.y = element_text(size = 15),
					   aspect.ratio = 0.5, plot.title = element_text(size =22))+
	geom_circle(hex[hex$node == 634, c('x', 'y')], 
		    r = 200, npoints = 500, linetype = 2, linewidth = 0.95) +
	scale_y_continuous(limits = c(1000, 1950)) + 
	geom_point(shape = 23, data = hex[hex$node == 634, c('x', 'y')], aes(x, y-10),
		   fill = 'purple3', size = 4)+
ggtitle('Recruits (Food)') + theme(plot.title = element_text(hjust = 0.5, vjust = -4))

panel_c_3 <- draw_traffic_flow(traffic_flow[type == 'Scout' & condition == 'NFD'], lineend = 'round',
			       add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round',
			       		    add = geom_foodpatches(fill = 'grey60')))+
	scale_size_continuous(range = c(1, 5)) + 
	scico::scale_color_scico() + theme(legend.position = 'none',
					   strip.text.y = element_text(size = 15),
					   aspect.ratio = 0.5, plot.title = element_text(size =22))+
	geom_circle(hex[hex$node == 634, c('x', 'y')], 
		    r = 350, npoints = 500, linetype = 2, linewidth = 0.95) +
	scale_y_continuous(limits = c(1000, 1950)) + 
	geom_point(shape = 23, data = hex[hex$node == 634, c('x', 'y')], aes(x, y-10),
		   fill = 'purple3', size = 4)+
	
	ggtitle('Scouts (No-Food)') + theme(plot.title = element_text(hjust = 0.5, vjust = -4))

panel_c_4 <- draw_traffic_flow(traffic_flow[type == 'Recruit' & condition == 'NFD' & N > 11], lineend = 'round',
			       add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round',
			       		    add = geom_foodpatches(fill = 'grey60')))+
	scale_size_continuous(range = c(1, 5)) + 
	scico::scale_color_scico() + theme(legend.position = 'none',
					   strip.text.y = element_text(size = 15),
					   aspect.ratio = 0.5, plot.title = element_text(size =22))+
	geom_circle(hex[hex$node == 634, c('x', 'y')], 
		    r = 350, npoints = 500, linetype = 2, linewidth = 0.95) +
	scale_y_continuous(limits = c(1000, 1950)) + 
	geom_point(shape = 23, data = hex[hex$node == 634, c('x', 'y')], aes(x, y-10),
		   fill = 'purple3', size = 4)+
	
	ggtitle('Recruits (No-Food)') + theme(plot.title = element_text(hjust = 0.5, vjust = -4))


panel_legend <- cowplot::get_legend(draw_traffic_flow(traffic_flow[type == 'Scout' & condition == 'NFD'], lineend = 'round')+
				    	scico::scale_color_scico('Density',
				    				 breaks = c(97, 2580), 
				    				 labels = c('Low', 'High')) +
				    	theme(legend.position = 'bottom',
				    	      legend.key.width = unit(1.5, 'cm'))+
				    	guides(color = guide_colorbar(title.position = 'top',
				    				      title.hjust = 0.5)))
grid::grid.newpage()


grid::grid.draw(panel_legend)
# 
# grid.arrange(ggarrange(subpanel_1+ theme(axis.line.x = element_blank(),
# 					 axis.text.x = element_blank(),
# 					 axis.ticks.x = element_blank()), 
# 					 subpanel_5+ theme(axis.line.x = element_blank(),
# 					 		  axis.text.x = element_blank(),
# 					 		  axis.ticks.x = element_blank()),
# 		       subpanel_3, subpanel_4), panel_c_1, panel_c_2, panel_c_3, panel_c_4,
# 	     layout_matrix = rbind(c(1, 1, 2, 2),
# 	     		      c(1, 1, 2, 2),
# 	     		      c(1, 1, 3, 3),
# 	     		      c(1, 1,  3, 3),
# 	     		      c(NA, NA, 4, 4),
# 	     		      c(NA, NA, 4, 4),
# 	     		      c(NA, NA, 5, 5),
# 	     		      c(NA, NA, 5, 5)))
# 
# panel_a <- (subpanel_1 + theme(axis.line.x = element_blank(),
# 			       axis.text.x = element_blank(),
# 			       axis.ticks.x = element_blank()) +
# 	    	subpanel_5+ theme(axis.line.x = element_blank(),
# 	    			  axis.text.x = element_blank(),
# 	    			  axis.ticks.x = element_blank())) /
# 	(subpanel_3+ subpanel_4)

