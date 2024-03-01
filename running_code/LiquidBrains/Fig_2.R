source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
source('~/research/gits/AnTracks/src/fit_functions.R')
load('~/research/gits/AnTracks/data/nf.RData')
load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/results/mov_rho_R.RData')


maxt <- 3000
maxd <- 500

library(trajr)
library(arrow)
library(egg)
library(latex2exp)

DET_scouts <- rbindlist(lapply(det, function(p){
	data <- setDT(p@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	set(data, j = 'd', value = as.numeric(dmatrix))
	LR <- unique(data[d > maxd, N_ind])
	SR <- unique(data[!N_ind %in% LR, N_ind])
	data_LR <- data[N_ind %in% LR]
	data_SR <- data[N_ind %in% SR]
	result_LR <- rbindlist(lapply(seq_along(LR), function(i){
		sbst <- data[N_ind == LR[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
		maxL <- max(Ls)
		if(maxL >= 100){
			idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
		} else {
			idxs <- nrow(xy)
		}
		straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
		data.frame(l = Ls[idxs], s = straightness, 
			   tag = paste0('det_LR_',i), scout = 'LR', condition = 'det')
		
	}), idcol = TRUE)
	result_SR <- rbindlist(lapply(seq_along(SR), function(i){
		sbst <- data[N_ind == SR[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 2){
			Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
			maxL <- max(Ls)
			if(maxL >= 100){
				idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
			} else {
				idxs <- nrow(xy)
			}
			straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
			data.frame(l = Ls[idxs], s = straightness, 
				   tag = paste0('det_SR_',i), scout = 'SR', condition = 'det')
		}
		
	}), idcol = TRUE)
	rbind(result_LR, result_SR)
}), idcol = TRUE)

NFD_scouts <- rbindlist(lapply(nf[-c(1, 2)], function(p){
	data <- setDT(p@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	set(data, j = 'd', value = as.numeric(dmatrix))
	LR <- unique(data[d > maxd, N_ind])
	SR <- unique(data[!N_ind %in% LR, N_ind])
	data_LR <- data[N_ind %in% LR]
	data_SR <- data[N_ind %in% SR]
	result_LR <- rbindlist(lapply(seq_along(LR), function(i){
		sbst <- data[N_ind == LR[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
		maxL <- max(Ls)
		if(maxL >= 100){
			idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
		} else {
			idxs <- nrow(xy)
		}
		straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
		data.frame(l = Ls[idxs], s = straightness, 
			   tag = paste0('nf_LR_',i), scout = 'LR', condition = 'nf')
		
	}), idcol = TRUE)
	result_SR <- rbindlist(lapply(seq_along(SR), function(i){
		sbst <- data[N_ind == SR[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 2){
			Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
			maxL <- max(Ls)
			if(maxL >= 100){
				idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
			} else {
				idxs <- nrow(xy)
			}
			straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
			data.frame(l = Ls[idxs], s = straightness, 
				   tag = paste0('nf_SR_',i), scout = 'SR', condition = 'nf')
		}
		
	}), idcol = TRUE)
	rbind(result_LR, result_SR)
}), idcol = TRUE)

data_list <- list(LR_det = DET_scouts[condition == 'det' & scout == 'LR',
				  .(y = mean(s)), by = .(x = round(l/100)*100)],
		  SR_det = DET_scouts[condition == 'det' & scout == 'SR',
		  		.(y = mean(s)), by = .(x = round(l/100)*100)],
		  LR_nf = NFD_scouts[condition == 'nf' & scout == 'LR',
		  		.(y = mean(s)), by = .(x = round(l/100)*100)],
		  SR_nf = NFD_scouts[condition == 'nf' & scout == 'SR',
		  		.(y = mean(s)), by = .(x = round(l/100)*100)])

fits <- lapply(seq_along(data_list), function(i){
	fit_doub_nexp(data_list[[i]], a_0 = 0.5, k_0 = 10**-3, a2_0 = 0.1, k2_0 = 10**-4,
		      alims = c(-10, 10), klims = c(-10, 10), a2lims = c(-10, 10), k2lims = c(-10, 10))
})
names(fits) <- names(data_list)
pred_data <- lapply(seq_along(data_list), function(i){
	data.table(x = data_list[[i]][['x']], 
		   y = doub_nexp(data_list[[i]][['x']], a = fits[[i]][['a']], 
		   	      k = fits[[i]][['k']], a2 = fits[[i]][['a2']], k2 = fits[[i]][['k2']]))
})
names(pred_data) <- names(data_list)


panel_C <- ggplot(data = rbindlist(data_list, idcol = TRUE),
       aes(x/10, y, color = factor(.id)), alpha = 0.5) + 
	geom_point(size = 4, shape = 21) +
	geom_path(data = rbindlist(pred_data, idcol = TRUE), linewidth = 1.5)+
	scale_y_continuous('<Net-to-Gross ratio>', limits = c(0, 1), breaks = seq(0, 1, 0.2))+
	scale_x_continuous('Travelled distance (cm)')+
	scale_color_manual('', values = c('mediumpurple', 'gold3', 'brown4', 'grey'), 
			   labels = c('LR DET', 'LR NFD', 'SR DET', 'SR NFD'))+
	theme(legend.title = element_blank(), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.position = c(3/4, 3/4))


grob_inset <- ggplot(data = mov_rho_R[rho %in% seq(0, 1, 0.25), .(N = mean(N)), by = c('rho', 'R')],
       aes(R*5, N, color = factor(rho))) +
	geom_point(size = 5) + geom_path(show.legend = FALSE, linewidth = 1.25)+
	ylab('Mean first passage (Iterations)') + 
	scale_x_continuous('Radius (cm)', breaks = seq(0, 20, 5)*5, limits = c(0, 20)*5)+
	scale_color_viridis_d('Proportion of LR')+
	theme(legend.position = c(0.2, 2/3), 
	      legend.background = element_rect(fill =NA, color = 'black'), 
	      plot.background = element_rect(fill = NA, color = 'black'))


panel_D <- ggplot(data = mov_rho_R[, .(N = mean(N), sd = sd(N)), by = c('scout', 'R')],
       aes(R*5, N, color = factor(scout))) +
	geom_pointrange(aes(ymin = N - sd, ymax = N + sd), size = 2, alpha = 0.75,
			linewidth = 3) +
	geom_path(show.legend = FALSE)+
	ylab('Mean first passage (Iterations)') + 
	scale_x_continuous('Radius (cm)', breaks = seq(0, 20, 5)*5, limits = c(0, 20)*5)+
	scale_color_viridis_d('Scout behaviour', direction = -1, end = 0.85)+
	theme(legend.position = c(0.8, 4/5), 
	      legend.background = element_rect(fill =NA, color = 'black'),
	      legend.key.size = unit(1, 'cm'))+
	guides(color = guide_legend(override.aes = list(alpha = 1, linetype = 0)))+
	annotation_custom(grob = ggplotGrob(grob_inset), xmin = 0, xmax = 12*5, ymin = 200, ymax = 500)

det_scouts <- rbindlist(lapply(det, function(p){
	data <- setDT(p@data)[Frame <= maxt]
	filter <- data[, .N, by = 'N_ind']
	data <- data[N_ind %in% filter[N >= 10, N_ind]]
	data[['d']] <- get_segment(data[, c('Xmm', 'Ymm')])[, 2]
	colnames(data)[colnames(data) == 'node'] <- 'o'
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	data[['dist']] <- as.numeric(dmatrix)
	LR <- unique(data[dist > maxd, N_ind])
	SR <- unique(data[!N_ind %in% LR, N_ind])
	data_LR <- data[N_ind %in% LR][, scout := 'LR']
	data_SR <- data[N_ind %in% SR][, scout := 'SR']
	rbindlist(list(LR = data_LR, SR = data_SR))
}), idcol = TRUE)
det_scouts[['id']] <- apply(det_scouts[, c('.id', 'N_ind')], 1, paste, collapse = '_', sep = '_')

nf_scouts <- rbindlist(lapply(nf[-c(1, 2)], function(p){ ## filter first two exps
	data <- setDT(p@data)[Frame <= maxt]
	filter <- data[, .N, by = 'N_ind']
	data <- data[N_ind %in% filter[N >= 10, N_ind]]
	data[['d']] <- get_segment(data[, c('Xmm', 'Ymm')])[, 2]
	colnames(data)[colnames(data) == 'node'] <- 'o'
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	data[['dist']] <- as.numeric(dmatrix)
	LR <- unique(data[dist > maxd, N_ind])
	SR <- unique(data[!N_ind %in% LR, N_ind])
	data_LR <- data[N_ind %in% LR][, scout := 'LR']
	data_SR <- data[N_ind %in% SR][, scout := 'SR']
	rbindlist(list(LR = data_LR, SR = data_SR))
}), idcol = TRUE)
nf_scouts[['id']] <- apply(nf_scouts[, c('.id', 'N_ind')], 1, paste, collapse = '_', sep = '_')

det_density <- det_scouts[, .(N = .N), by = c('scout','o','d')][, condition := 'DET'][, z := rank(N)]
nf_density <- nf_scouts[, .(N = .N), by = c('scout','o','d')][, condition := 'NFD'][, z := rank(N)]
traffic_flow <- merge(rbind(det_density, nf_density), edges)

LR_traffic <- merge(traffic_flow[scout == 'LR'], edges)
SR_traffic <- merge(traffic_flow[scout == 'SR'], edges)

panel_B <- draw_traffic_flow(traffic_flow, lineend = 'round',
		  add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round',
		  		    add = geom_foodpatches(fill = 'mediumpurple')))+
	scale_size_continuous(range = c(1, 5)) + 
	scico::scale_color_scico() + theme(legend.position = 'none',
					   strip.text.y = element_text(size = 15,
					   			    margin = margin(r = 5, l = 5)),
					   aspect.ratio = 0.5)+
	facet_grid(scout ~ condition)


DET_metrics <- rbindlist(lapply(det, function(p){
	data <- setDT(p@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	set(data, j = 'd', value = as.numeric(dmatrix))
	LR <- unique(data[d > maxd, N_ind])
	SR <- unique(data[!N_ind %in% LR, N_ind])
	data_LR <- data[N_ind %in% LR]
	data_SR <- data[N_ind %in% SR]
	result_LR <- rbindlist(lapply(seq_along(LR), function(i){
		sbst <- data[N_ind == LR[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 10){
			data.frame(s = TrajStraightness(xy),
				   s2 = Mod(TrajMeanVectorOfTurningAngles(xy)),
				   v = mean(Mod(TrajVelocity(xy)), na.rm = TRUE),
				   v2 = Mod(TrajMeanVelocity(xy)),
				   acc = mean(Mod(TrajAcceleration(xy)), na.rm = TRUE),
				   l = TrajLength(xy),
				   d = TrajDistance(xy),
				   t = TrajDuration(xy),
				   sin = TrajSinuosity2(xy), 
				   tag = paste0('det_LR_',i), scout = 'LR', condition = 'det')
		}
	}))
	result_SR <- rbindlist(lapply(seq_along(SR), function(i){
		sbst <- data[N_ind == SR[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 10){
			data.frame(s = TrajStraightness(xy),
				   s2 = Mod(TrajMeanVectorOfTurningAngles(xy)),
				   v = mean(Mod(TrajVelocity(xy)), na.rm = TRUE),
				   v2 = Mod(TrajMeanVelocity(xy)),
				   acc = mean(Mod(TrajAcceleration(xy)), na.rm = TRUE),
				   l = TrajLength(xy),
				   d = TrajDistance(xy),
				   t = TrajDuration(xy),
				   sin = TrajSinuosity2(xy), 
				   tag = paste0('det_SR_',i), scout = 'SR', condition = 'det')
		}

		
	}))
	rbind(result_LR, result_SR)
}), idcol = TRUE)

NFD_metrics <- rbindlist(lapply(nf, function(p){
	data <- setDT(p@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	set(data, j = 'd', value = as.numeric(dmatrix))
	LR <- unique(data[d > maxd, N_ind])
	SR <- unique(data[!N_ind %in% LR, N_ind])
	data_LR <- data[N_ind %in% LR]
	data_SR <- data[N_ind %in% SR]
	result_LR <- rbindlist(lapply(seq_along(LR), function(i){
		sbst <- data[N_ind == LR[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 10){
			data.frame(s = TrajStraightness(xy),
				   s2 = Mod(TrajMeanVectorOfTurningAngles(xy)),
				   v = mean(Mod(TrajVelocity(xy)), na.rm = TRUE),
				   v2 = Mod(TrajMeanVelocity(xy)),
				   acc = mean(Mod(TrajAcceleration(xy)), na.rm = TRUE),
				   l = TrajLength(xy),
				   d = TrajDistance(xy),
				   t = TrajDuration(xy),
				   sin = TrajSinuosity2(xy), 
				   tag = paste0('nf_LR_',i), scout = 'LR', condition = 'nf')
		}
	}))
	result_SR <- rbindlist(lapply(seq_along(SR), function(i){
		sbst <- data[N_ind == SR[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 10){
			data.frame(s = TrajStraightness(xy),
				   s2 = Mod(TrajMeanVectorOfTurningAngles(xy)),
				   v = mean(Mod(TrajVelocity(xy)), na.rm = TRUE),
				   v2 = Mod(TrajMeanVelocity(xy)),
				   acc = mean(Mod(TrajAcceleration(xy)), na.rm = TRUE),
				   l = TrajLength(xy),
				   d = TrajDistance(xy),
				   t = TrajDuration(xy),
				   sin = TrajSinuosity2(xy), 
				   tag = paste0('nf_SR_',i), scout = 'SR', condition = 'nf')
		}
		
		
	}))
	rbind(result_LR, result_SR)
}), idcol = TRUE)

mlt_det_metrics <- data.table(melt(DET_metrics, id.vars = c('tag', 'scout', 'condition', '.id')))
mlt_nf_metrics <- data.table(melt(NFD_metrics, id.vars = c('tag', 'scout', 'condition', '.id')))

## selected cols 
sdcols <- c('v', 'acc', 'l', 'd', 't')
general_metrics <- rbind(mlt_det_metrics[variable %in% sdcols], mlt_nf_metrics[variable %in% sdcols])
general_metrics[variable %in% sdcols[-5], value := value / 10]
general_metrics[variable %in% sdcols[5], value := value / 60]
# variable = c(v = TeX("Velocity ($cm*s^{-1}$)"),
# 	     acc = TeX("Acceleration ($cm*s^{-2}$)"),
# 	     l ="Cumulative distance (cm)",
# 	     d ='Max distance (cm)',
# 	     t = 'Time in arena (min)'), .default = label_parsed)

fcking_labeller <- function(string){
	TeX(paste0(string))
}
panel_A <- ggplot(data = general_metrics, aes(scout, value, fill = condition)) + 
	facet_wrap(~ factor(variable, 
			    labels = c("Velocity ($cm\\cdot s^{-1}$)",
						 "Acceleration ($cm\\cdot s^{-2}$)",
						 "`Cumulative distance (cm)`",
						 '`Max distance (cm)`',
						 '`Time in arena (min)`')),
		   scales = 'free_y', 
		   labeller = as_labeller(fcking_labeller, default = label_parsed)) + 
	geom_boxplot(alpha = 0.6, outlier.shape = NA)+
	geom_jitter(size = 3, alpha = 0.1, 
		    position = position_jitterdodge(dodge.width = 0.75,jitter.width = 0.15))+
	scale_fill_manual('Experimental condition', labels = c('DET', 'NFD'), values = c('mediumpurple','gold3'))+
	ylab('')+xlab('') +
	theme(legend.position = c(0.8, 0.25), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.justification = 'center')

grid.arrange(panel_A, panel_B, panel_C, panel_D, layout_matrix = rbind(c(rep(1, 5), rep(2, 5)),
								       c(rep(1, 5), rep(2, 5)),
								       c(rep(3, 5), rep(4, 5)),
								       c(rep(3, 5), rep(4, 5))))
