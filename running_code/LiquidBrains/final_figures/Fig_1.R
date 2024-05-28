source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
source('~/research/gits/AnTracks/src/fit_functions.R')
load('~/research/gits/AnTracks/data/nf.RData')
load('~/research/gits/AnTracks/data/det.RData')

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

maxt <- max(sapply(det, function(i){
	     min(rbindlist(i@food)[['t']]) - min(i@data[['Frame']])
	})) # 3000
maxd <- 500 # 500

library(trajr)
library(arrow)
library(egg)
library(latex2exp)

DET_scouts <- rbindlist(lapply(det, function(p){
	data <- setDT(p@data)[Frame <= (maxt + min(Frame))]
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
	# data <- setDT(p@data)[Frame <= (maxt + min(Frame))]
	# dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	# set(data, j = 'd', value = as.numeric(dmatrix))
	# LR <- unique(data[d > maxd, N_ind])
	# SR <- unique(data[!N_ind %in% LR, N_ind])
	# data_LR <- data[N_ind %in% LR]
	# data_SR <- data[N_ind %in% SR]
	data <- setDT(p@data)[Frame < maxt + min(Frame)]
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

data_list <- rbindlist(data_list, idcol = TRUE)
labels_unique <- do.call('rbind', strsplit(data_list[['.id']], '_'))
data_list[['scout']] <- labels_unique[, 1]
data_list[['condition']] <- labels_unique[, 2]

pred_data <- rbindlist(pred_data)
pred_data[['scout']] <- labels_unique[, 1]
pred_data[['condition']] <- labels_unique[, 2]

panel_B <- ggplot(data = data_list,
		  aes(x/10, y, color = factor(condition, labels = c('DET', 'NFD')),
		      shape = factor(scout, labels = c('LR', 'SR'))), alpha = 0.5) + 
	geom_point(size = 4) +
	geom_path(data = pred_data, linewidth = 1.5)+
	scale_y_continuous('<Net-to-Gross ratio>', limits = c(0, 1), breaks = seq(0, 1, 0.2))+
	scale_x_continuous('Travelled distance (cm)')+
	scale_color_manual('', values = c('mediumpurple', 'gold3'), 
			   labels = c('DET', 'NFD'))+
	scale_shape_manual('', values = c(21, 24))+
	theme(legend.title = element_blank(), 
	      legend.position = c(3/4, 3/4),
	      legend.box.background = element_rect(),
	      legend.background = element_blank(), plot.title = element_text(size = 22))+
	guides(shape = guide_legend(override.aes = list(shape = c(19, 17))))


det_scouts <- rbindlist(lapply(det, function(p){
	data <- setDT(p@data)[Frame <= (maxt + min(Frame))]
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
	data <- setDT(p@data)[Frame <= (maxt + min(Frame))]
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

# panel_D <- draw_traffic_flow(traffic_flow, lineend = 'round',
# 			     add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round',
# 			     		    add = geom_foodpatches(fill = 'mediumpurple')))+
# 	scale_size_continuous(range = c(1, 5)) + 
# 	scico::scale_color_scico() + theme(legend.position = 'none',
# 					   strip.text.y = element_text(size = 15,
# 					   			    margin = margin(r = 5, l = 5)),
# 					   aspect.ratio = 0.5)+
# 	facet_grid(scout ~ condition)
panel_D <- draw_traffic_flow(traffic_flow, lineend = 'round',
			     add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round',
			     		    add = geom_foodpatches(fill = 'mediumpurple')))+
	scale_size_continuous(range = c(1, 5)) + 
	scico::scale_color_scico() + theme(legend.position = 'none',
					   strip.text.y = element_text(size = 15,
					   			    margin = margin(r = 5, l = 5)),
					   aspect.ratio = 0.5, plot.title = element_text(size =22))+
	facet_grid(scout ~ condition)+ 
	geom_circle(hex[hex$node == 634, c('x', 'y')], 
		    r = maxd, npoints = 500, linetype = 2, linewidth = 0.95) +
	scale_y_continuous(limits = c(1000, 1950)) + 
	geom_point(shape = 24, data = hex[hex$node == 634, c('x', 'y')], aes(x, y-10),
		   fill = 'purple3', size = 4)


DET_metrics <- rbindlist(lapply(det, function(p){
	data <- setDT(p@data)[Frame <= (maxt + min(Frame))]
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
				   d = TrajDistance(xy, endIndex = which.max(sbst[['d']])),
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
				   d = TrajDistance(xy, endIndex = which.max(sbst[['d']])),
				   t = TrajDuration(xy),
				   sin = TrajSinuosity2(xy), 
				   tag = paste0('det_SR_',i), scout = 'SR', condition = 'det')
		}
		
		
	}))
	rbind(result_LR, result_SR)
}), idcol = TRUE)

NFD_metrics <- rbindlist(lapply(nf, function(p){
	data <- setDT(p@data)[Frame <= (maxt + min(Frame))]
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
				   d = TrajDistance(xy, endIndex = which.max(sbst[['d']])),
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
				   d = TrajDistance(xy, endIndex = which.max(sbst[['d']])),
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
# sdcols <- c('v', 'acc', 'l', 'd', 't')
sdcols <- c('v', 'acc', 'd', 't')
general_metrics <- rbind(mlt_det_metrics[variable %in% sdcols], mlt_nf_metrics[variable %in% sdcols])
general_metrics[variable %in% sdcols[-4], value := value / 10]
general_metrics[variable %in% sdcols[4], value := value / 60]
# variable = c(v = TeX("Velocity ($cm*s^{-1}$)"),
# 	     acc = TeX("Acceleration ($cm*s^{-2}$)"),
# 	     l ="Cumulative distance (cm)",
# 	     d ='Max distance (cm)',
# 	     t = 'Time in arena (min)'), .default = label_parsed)

labeller <- function(string){
	TeX(paste0(string))
}

## CODE TO OBTAIN THE PAIRWISED WILCOXON TESTS
# tests <- general_metrics
# tests[['label']] <- apply(tests[, c('scout', 'condition')], 1, paste0, collapse = '', sep='')
# 
# for(v in sdcols){
# 	print(v)
# 	print(
# 	pairwise.wilcox.test(x = tests[variable == v, value], 
# 			     g= tests[variable == v, label], p.adjust.method = 'bonferroni'))
# }


pvalues <- data.table(scout = rep(c('LR', 'SR'), 8), 
	   condition = rep(c(rep('det', 2), rep('nf', 2)), 4),
	   label = c('a', 'a', 'b', 'b',
	   	  'a', 'bc', 'ab', 'c',
	   	  'a', 'a', 'b', 'c',
	   	  'a', 'a', 'b', 'b'),
	   variable = c(vapply(sdcols, rep,FUN.VALUE = character(4), 4), recursive = T),
	   value = NA)
pvalues[['y']] <- sapply(1:nrow(pvalues), function(i) general_metrics[scout == pvalues[i, scout] & 
					 	variable == pvalues[i, variable]&
					 	condition == pvalues[i, condition],
					 min(quantile(value, p = 0.75) + 1.5 * (quantile(value, p = 0.75)-quantile(value, p = 0.25)),
					     max(value))  + median(value)* 0.1])
pvalues[['y']] <- pvalues[['y']] + c(rep(0.25, 4), 0.5, 1, 1,1, 7.5,20, 10, 15, 0.5,1.5, 1.25, 1.25)
# pvalues[['y']] <- pvalues[['y']] + c(rep(0, 5), 0.5, 0.5,0.5, 0,10, 0, 10, 0,0.5, 0.5, 0.5)

# pvalues[variable == 't', 'y'] <- pvalues[variable == 't', y / 60]
# pvalues[variable == 'd', 'y'] <- pvalues[variable == 'd', y / 100]
pvalues[['variable']] <- factor(pvalues[['variable']], levels = sdcols)

panel_A <- ggplot(data = general_metrics, aes(scout, value, fill = condition)) + 
 
	geom_boxplot(alpha = 0.6, outlier.shape = NA)+
	geom_jitter(size = 3, alpha = 0.1, 
		    position = position_jitterdodge(dodge.width = 0.75,jitter.width = 0.15),
		    show.legend = FALSE)+
	## for some reason labels are sorted wrong...
	geom_text(data = pvalues, aes(label = label, y = pvalues[order(variable, scout), y]),
		  size = 7, position = position_nudge(x = c(-0.19, -1+0.19, 1-0.19, 0.19)))+
	scale_fill_manual('', labels = c('DET', 'NFD'), values = c('mediumpurple','gold3'))+
	ylab('')+xlab('') +
	facet_wrap(~ factor(variable, 
			    labels = c("Velocity ($cm\\cdot s^{-1}$)",
			    	   "Acceleration ($cm\\cdot s^{-2}$)",
			    	   # "`Cumulative distance (cm)`",
			    	   '`Maximum distance (cm)`',
			    	   '`Time in arena (min)`')),
		   scales = 'free_y', 
		   labeller = as_labeller(labeller, default = label_parsed)) +
	theme(legend.position = c(0.35, 0.35), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.justification = 'center', legend.title = element_blank(),
	      plot.title = element_text(size = 22))


ref_raster <- expand.grid(c('R', 'B', 'L'), c('R', 'B', 'L'))
colnames(ref_raster) <- c('x', 'y')

## DET LR
mm_DET_LR = matrix(data = c(0.3770950,0.1703911, 0.4525140, 
			    0.46616541, 0.06015038, 0.47368421,
			    0.4373259, 0.1838440, 0.3788301),
		   nrow = 3, ncol = 3, 
		   dimnames = list(c('R', 'B', 'L'), c('R', 'B', 'L')))
apply(mm_DET_LR, 1, sum) / sum(mm_DET_LR)

## DET SR 
mm_DET_SR = matrix(data = c(0.3571429, 0.2362637, 0.4065934,
			    0.3928571,  0.0625000,  0.5446429,
			    0.3594470, 0.2580645, 0.3824885),
		   nrow = 3, ncol = 3, 
		   dimnames = list(c('R', 'B', 'L'), c('R', 'B', 'L')))
apply(mm_DET_SR, 1, sum) / sum(mm_DET_SR)

## NFD LR
mm_NFD_LR = matrix(data = c(0.3622642, 0.1716981, 0.4660377,
			    0.4120879, 0.1098901, 0.4780220,
			    0.4904580, 0.1335878, 0.3759542),
		   nrow = 3, ncol = 3, 
		   dimnames = list(c('R', 'B', 'L'), c('R', 'B', 'L')))
apply(mm_NFD_LR, 1, sum) / sum(mm_NFD_LR)

## NFD SR
mm_NFD_SR = matrix(data = c(0.3678571,0.2821429, 0.3500000,
			    0.58273381, 0.04316547, 0.37410072,
			    0.4377880, 0.2488479, 0.3133641),
		   nrow = 3, ncol = 3, 
		   dimnames = list(c('R', 'B', 'L'), c('R', 'B', 'L')))
apply(mm_NFD_SR, 1, sum) / sum(mm_NFD_SR)


full_raster <- rbindlist(replicate(4, ref_raster, simplify = FALSE))
full_raster[['z']] <- c(mm_DET_LR, mm_DET_SR, mm_NFD_LR, mm_NFD_SR, recursive =TRUE)
full_raster[['exp']] <- c(rep('DET', 18),
			  rep('NFD', 18))
full_raster[['scout']] <- c(rep('LR', 9), rep('SR', 9), rep('LR', 9), rep('SR', 9))

panel_C <- ggplot(data = full_raster, aes(x, y, fill = z)) + 
	geom_tile(color = 'black', show.legend = FALSE) +
	facet_grid(scout ~ exp)+
	geom_text(aes(label = format(round(z, 3), nsmall = 3)), size = 6)+
	scale_fill_gradient(limits = c(0, 1), low = 'white', high = scales::muted('blue'))+
	theme(axis.title = element_blank(),
	      strip.text.y = element_text(size = 15, margin = margin(l = 5, r = 5)),
	      plot.margin = margin(l = 35, b = 35), plot.title = element_text(size = 22))

layout_mat <- rbind(do.call('rbind', replicate(5, list(c(rep(1, 4), rep(2, 6))))),
      c(rep(1, 4), rep(4, 6)),
      do.call('rbind', replicate(5, list(c(rep(3, 4), rep(4, 6))))))

# rbind(c(rep(1, 4), rep(2, 6)),
#       c(rep(1, 4), rep(2, 6)),
#       c(rep(1, 4), rep(4, 6)),
#       c(rep(3, 4), rep(4, 6)),
#       c(rep(3, 4), rep(4, 6)))

png('~/research/gits/AnTracks/plots/figures_LiquidBrains/Fig_1.png', 6000, 4000, res = 350)
grid.arrange(panel_A + ggtitle('A'),
	     panel_B + ggtitle('B'),
	     panel_C + ggtitle('C'),
	     panel_D + ggtitle('D'), layout_matrix = layout_mat)
								       
dev.off()

# grid.arrange(panel_A, panel_B, panel_C, panel_D, layout_matrix = rbind(c(rep(1, 5), rep(2, 5)),
# 								       c(rep(1, 5), rep(2, 5)),
# 								       c(rep(3, 5), rep(4, 5)),
# 								       c(rep(3, 5), rep(4, 5))))
