source('~/research/gits/AnTracks/src/Experiment.R')
load('~/research/gits/AnTracks/data/nf.RData')
load('~/research/gits/AnTracks/data/det.RData')



ggplot(data = general_det, aes(d)) +
	# geom_line(stat = 'density', size = 1, aes(fill = type))
	geom_density(size = 1, aes(fill = type, color = type), alpha = 0.5, bw = bin_res)+
	geom_vline(xintercept = 200, linetype = 2, linewidth= 1)+
	scale_x_continuous('Distance to the nest (cm)', breaks = seq(0, ceiling(general_det[, max(d)] / 100) * 100, bin_res*3),
			   labels = seq(0, ceiling(general_det[, max(d)] / 100) * 100, bin_res*3)/10)+
	ylab('Density')+
	scale_fill_viridis_d('', direction = -1, end = 0.85) +
	scale_color_viridis_d('', direction = -1, end = 0.85) 

ggplot(data = general_det, aes(d)) +
	# geom_line(stat = 'density', size = 1, aes(fill = type))
	geom_density(size = 1, aes(fill = type, color = type), alpha = 0.5, bw = bin_res)+
	geom_vline(xintercept = 200, linetype = 2, linewidth= 1)+
	scale_x_continuous('Distance to the nest (cm)', breaks = seq(0, ceiling(general_det[, max(d)] / 100) * 100, bin_res*3),
			   labels = seq(0, ceiling(general_det[, max(d)] / 100) * 100, bin_res*3)/10)+
	ylab('Density')+
	scale_fill_viridis_d('', direction = -1, end = 0.85) +
	scale_color_viridis_d('', direction = -1, end = 0.85) 



nf <- lapply(nf[-c(1, 2)], class_ids)
general_nf <- rbindlist(lapply(nf, function(i){
	i@data[type != 'UNKNOWN']
}), idcol = 'exp')
general_nf[['d']] <- as.numeric(pdist(as.matrix(general_nf[, c('Xmm', 'Ymm')]), 
				       as.matrix(hex[hex$node == 634, c('x', 'y')])))

general_nf[['bin']] <- cut(general_nf[['d']],
			    breaks = seq(0, ceiling(general_nf[, max(d)] / 100) * 100, bin_res), 
			    right = FALSE, include.lowest = TRUE)


library(patchwork)


ggplot(data = general_det, aes(d)) +
	# geom_line(stat = 'density', size = 1, aes(fill = type))
	geom_density(size = 1, aes(fill = type, color = type), alpha = 0.5, bw = bin_res,
		     show.legend = FALSE)+
	geom_vline(xintercept = 200, linetype = 2, linewidth= 1)+
	scale_x_continuous('Distance to the nest (cm)', breaks = seq(0, ceiling(general_det[, max(d)] / 100) * 100, bin_res*3),
			   labels = seq(0, ceiling(general_det[, max(d)] / 100) * 100, bin_res*3)/10)+
	ylab('Density')+
	scale_fill_viridis_d('', direction = -1, end = 0.85) +
	scale_color_viridis_d('', direction = -1, end = 0.85) +
	theme(aspect.ratio = 0.6)


library(svglite)
svg('~/research/gits/AnTracks/plots/figures_LiquidBrains/fig_SX_site_fidelity.svg', 
    width = 20, height = 12, pointsize = 50)
ggplot(data = general_nf, aes(d)) +
	# geom_line(stat = 'density', size = 1, aes(fill = type))
	geom_density(size = 1, aes(fill = type, color = type), alpha = 0.5, bw = bin_res)+
	geom_vline(xintercept = 350, linetype = 2, linewidth= 1)+
	scale_x_continuous('Distance to the nest (cm)', breaks = seq(0, ceiling(general_nf[, max(d)] / 100) * 100, bin_res*3),
			   labels = seq(0, ceiling(general_nf[, max(d)] / 100) * 100, bin_res*3)/10)+
	ylab('Density')+
	scale_fill_viridis_d('', direction = -1, end = 0.85) +
	scale_color_viridis_d('', direction = -1, end = 0.85) +
	theme(aspect.ratio = 0.6, 
	      legend.position = c(0.75, 0.75))+
	guides(fill = guide_legend(override.aes = (list(alpha = 1))))
dev.off()




ggplot(data = rbindlist(list(det = general_det,
			     nf = general_nf), idcol = 'condition'), aes(d)) +
	# geom_line(stat = 'density', size = 1, aes(fill = type))
	geom_density(size = 1, aes(fill = factor(type, levels = c('Scout', 'Recruit'),
						 labels = c('Scouts', 'Recruits')),
				   color = factor(type, levels = c('Scout', 'Recruit'),
				   	       labels = c('Scouts', 'Recruits'))), alpha = 0.5, bw = bin_res)+
	geom_vline(aes(xintercept = splitd),
		   data = data.frame(splitd = c(200, 350),
		   		  condition = c('det', 'nf')), linetype = 2, linewidth= 1)+
	scale_x_continuous('Distance to the nest (cm)', breaks = seq(0, ceiling(general_nf[, max(d)] / 100) * 100, bin_res*4),
			   labels = seq(0, ceiling(general_nf[, max(d)] / 100) * 100, bin_res*4)/10)+
	ylab('Density')+
	facet_wrap(~ factor(condition, 
			    labels = c('Experimental', 'Control')))+
	scale_fill_viridis_d('', direction = -1, end = 0.85) +
	scale_color_viridis_d('', direction = -1, end = 0.85)  +
	theme(aspect.ratio = 0.75, 
	      legend.position = c(0.825, 0.75))+
	guides(fill = guide_legend(override.aes = (list(alpha = 1))))

####### OLD CODE -------- #########


maxt <- max(sapply(det, function(i){
	min(rbindlist(i@food)[['t']]) - min(i@data[['Frame']])
})) # 3000
maxd <- 500 # 500

library(trajr)
library(arrow)
library(egg)
library(latex2exp)

R_seq <- seq(100, 600, 20)

DET_fidelity <- rbindlist(lapply(det, function(i){
	data <- class_ids(i)@data[type != 'UNKNOWN']
	d <- as.numeric(pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')])))
	data[['d']] <- d
	inds <- data[, unique(N_ind)]
	rbindlist(lapply(inds, function(x){
		p <- vapply(R_seq, function(r){
			nrow(data[N_ind == x & d < r])
		}, numeric(1)) / nrow(data[N_ind == x]) 
		data.frame(r = R_seq, p = p, type = data[N_ind == x, unique(type)])
	}), idcol = 'N_ind')
}), idcol = 'exp')
DET_fidelity[['tag']] <- apply(DET_fidelity[, c('exp','N_ind')], 
			       1, paste0, collapse = '_', sep = '')

ggplot(data = DET_fidelity[order(r, p)],
       aes(r, 100*p, group = tag, fill = type, shape = type, color = type)) + 
	geom_line(show.legend = FALSE)+
	geom_point(size = 4, color = 'black',
		   alpha = 0.75) + 
	scale_y_continuous('Site fidelity (%)')+
	scale_x_continuous('Radius (mm)', breaks = seq(0, 600, 150))+
	scale_shape_manual('', values = c(24, 21))+
	scale_fill_viridis_d('', direction = -1, end = 0.85)+
	scale_color_viridis_d('', direction = -1, end = 0.85)+
	theme(legend.title = element_blank(), 
	      legend.background = element_rect(fill = NA, color = 'black'))+
	geom_vline(xintercept = 200, linetype = 2)

NFD_fidelity <- rbindlist(lapply(nf[-c(1, 2)], function(i){
	data <- class_ids(i)@data[type != 'UNKNOWN']
	d <- as.numeric(pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')])))
	data[['d']] <- d
	inds <- data[, unique(N_ind)]
	rbindlist(lapply(inds, function(x){
		p <- vapply(R_seq, function(r){
			nrow(data[N_ind == x & d < r])
		}, numeric(1)) / nrow(data[N_ind == x]) 
		data.frame(r = R_seq, p = p, type = data[N_ind == x, unique(type)])
	}), idcol = 'N_ind')
}), idcol = 'exp')
NFD_fidelity[['tag']] <- apply(NFD_fidelity[, c('exp','N_ind')], 
			       1, paste0, collapse = '_', sep = '')

# ggplot(data = NFD_fidelity[order(r, p)],
#        aes(r, 100*p, group = tag, fill = type, shape = type, color = type)) + 
# 	geom_line(show.legend = FALSE)+
# 	geom_point(size = 4, alpha = 0.75) + 
# 	scale_y_continuous('Site fidelity (%)')+
# 	scale_x_continuous('Radius (mm)', breaks = seq(0, 600, 150))+
# 	scale_shape_manual('', values = c(24, 21))+
# 	scale_fill_viridis_d('', direction = -1, end = 0.85)+
# 	scale_color_viridis_d('', direction = -1, end = 0.85)+
# 	theme(legend.title = element_blank(), 
# 	      legend.background = element_rect(fill = NA, color = 'black'))+
# 	geom_vline(xintercept = 350, linetype = 2)

ggplot(data = NFD_fidelity[order(r, p)],
       aes(r, 100*p, group = tag, fill = type, shape = type, color = type)) + 
	geom_line(show.legend = FALSE, alpha = 0.4)+
	geom_point(size = 1, alpha = 0.4) + 
	scale_y_continuous('Site fidelity (%)')+
	scale_x_continuous('Radius (mm)', breaks = seq(0, 600, 150))+
	scale_shape_manual('', values = c(24, 21))+
	scale_fill_viridis_d('', direction = -1, end = 0.85)+
	scale_color_viridis_d('', direction = -1, end = 0.85)+
	theme(legend.title = element_blank(), 
	      legend.background = element_rect(fill = NA, color = 'black'))+
	guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1)))+
	geom_vline(xintercept = 350, linetype = 2, linewidth = 2)



png('/home/polfer/research/gits/AnTracks/plots/figures_LiquidBrains/Fig_S1.png', 6000, 3000, res = 350)
ggplot(data = rbind(DET_scouts, NFD_scouts)[order(R, p), c('R', 'p', 'tag', 'scout', 'condition')],
       aes(R, 100*p, group = tag, color = scout, shape = scout)) + 
	geom_point(size = 4) + geom_line(show.legend = FALSE)+
	facet_wrap(~ factor(condition, labels = c('DET', 'NFD')))+
	scale_y_continuous('Site fidelity (%)')+
	scale_x_continuous('Radius (mm)', breaks = seq(0, 600, 150))+
	theme(legend.title = element_blank(), legend.background = element_rect(fill = NA, color = 'black'))
dev.off()


bin_res <- 50
det <- lapply(det, class_ids)
general_det <- rbindlist(lapply(det, function(i){
	i@data[type != 'UNKNOWN']
}), idcol = 'exp')
general_det[['d']] <- as.numeric(pdist(as.matrix(general_det[, c('Xmm', 'Ymm')]), 
				       as.matrix(hex[hex$node == 634, c('x', 'y')])))

general_det[['bin']] <- cut(general_det[['d']],
			    breaks = seq(0, ceiling(general_det[, max(d)] / 100) * 100, bin_res), 
			    right = FALSE, include.lowest = TRUE)



ggplot(data = general_det[order(bin), .N, by = c('type', 'bin')],
       aes(bin, N, color = type)) + 
	geom_path(aes(group = type)) + geom_point()


type_counts <- general_det[order(bin), .N, by = c('type', 'bin')]
type_totals <- general_det[, .N, by = c('type')]
type_counts <- merge(type_counts, type_totals, by = 'type', suffixes = c("", "_total"))
type_counts[, proportion := N / N_total]

ggplot(data = type_counts,
       aes(bin, proportion, color = type)) + 
	geom_path(aes(group = type)) + geom_point()