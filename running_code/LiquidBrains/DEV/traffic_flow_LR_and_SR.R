source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
load('~/research/gits/AnTracks/data/nf.RData')
load('~/research/gits/AnTracks/data/det.RData')

library(trajr)

maxt <- 3000
maxd <- 500
LRs <- lapply(det, function(p){
	data <- setDT(p@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	set(data, j = 'd', value = as.numeric(dmatrix))
	LR <- unique(data[d > maxd, N_ind])
	SR <- unique(data[!N_ind %in% LR, N_ind])
	data_LR <- data[N_ind %in% LR]
	data_SR <- data[N_ind %in% SR]
	rbindlist(lapply(seq_along(LR), function(i){
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
		data.frame(l = Ls[idxs], s = straightness)
		
	}), idcol = TRUE)})
names(LRs) <- paste0('Det_', seq_along(det))
LRs <- rbindlist(LRs, idcol = TRUE)
LRs[['tag']] <- apply(LRs[, 1:2], 1, paste, collapse = '_')



# maxt <- 2400
maxd <- 500
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

# nf_scouts <- rbindlist(lapply(nf, function(p){
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
traffic_flow <- rbind(det_density, nf_density)

inds_det <- det_scouts[, c('id', 'scout', 'o', 'd')][, N := .N, by = c('id', 'o', 'd', 'scout')][, z:=rank(N)]



###### LR -- DET AND NFD #######
density_LR <- data.table(apply(do.call('rbind', 
				 strsplit(rle(apply(det_scouts[scout == 'LR', c('o','d')], 1,
				 		   paste, collapse = '_', sep = '_'))$values, '_')), 2, as.numeric))
colnames(density_LR) <- c('o','d')

data_flow <- merge(edges, density_LR[, N := .N, by = c('o', 'd')][, z := rank(N)])
draw_traffic_flow(data_flow, lineend = 'round',
		  add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round',
		  		    add = geom_foodpatches(fill = 'mediumpurple')))+
	scale_size_continuous(range = c(1, 5)) + 
	scico::scale_color_scico() + theme(legend.position = 'none')+
	ggtitle('LR - DET')



density_LR <- data.table(apply(do.call('rbind', 
				       strsplit(rle(apply(nf_scouts[scout == 'LR', c('o','d')], 1,
				       		   paste, collapse = '_', sep = '_'))$values, '_')), 2, as.numeric))
colnames(density_LR) <- c('o','d')

data_flow <- merge(edges, density_LR[, N := .N, by = c('o', 'd')][, z := rank(N)])
draw_traffic_flow(data_flow, lineend = 'round',
		  add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round',
		  		    add = geom_foodpatches(fill = 'mediumpurple')))+
	scale_size_continuous(range = c(1, 5)) + 
	scico::scale_color_scico() + theme(legend.position = 'none')+
	ggtitle('LR - NFD')

###### SR -- DET AND NFD #######
density_SR <- data.table(apply(do.call('rbind', 
				       strsplit(rle(apply(det_scouts[scout == 'SR', c('o','d')], 1,
				       		   paste, collapse = '_', sep = '_'))$values, '_')), 2, as.numeric))
colnames(density_SR) <- c('o','d')

data_flow <- merge(edges, density_SR[, N := .N, by = c('o', 'd')][, z := rank(N)])
draw_traffic_flow(data_flow, lineend = 'round',
		  add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round',
		  		    add = geom_foodpatches(fill = 'mediumpurple')))+
	scale_size_continuous(range = c(1, 5)) + 
	scico::scale_color_scico() + theme(legend.position = 'none')+
	ggtitle('SR - DET')



density_SR <- data.table(apply(do.call('rbind', 
				       strsplit(rle(apply(nf_scouts[scout == 'SR', c('o','d')], 1,
				       		   paste, collapse = '_', sep = '_'))$values, '_')), 2, as.numeric))
colnames(density_SR) <- c('o','d')

data_flow <- merge(edges, density_SR[, N := .N, by = c('o', 'd')][, z := rank(N)])
draw_traffic_flow(data_flow, lineend = 'round',
		  add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round',
		  		    add = geom_foodpatches(fill = 'mediumpurple')))+
	scale_size_continuous(range = c(1, 5)) + 
	scico::scale_color_scico() + theme(legend.position = 'none')+
	ggtitle('SR - NFD')













LR_traffic <- merge(traffic_flow[scout == 'LR'], edges)
SR_traffic <- merge(traffic_flow[scout == 'SR'], edges)

grid.arrange(
draw_traffic_flow(LR_traffic[condition == 'DET'], lineend = 'round',
		  add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round',
		  		    add = geom_foodpatches(fill = 'mediumpurple')))+
	scale_size_continuous(range = c(1, 5)) + 
	scico::scale_color_scico() + theme(legend.position = 'none')+
	ggtitle('LR - DET'),
draw_traffic_flow(LR_traffic[condition == 'NFD'], lineend = 'round',
		  add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round'))+
	scale_size_continuous(range = c(1, 5)) + 
	scico::scale_color_scico() + theme(legend.position = 'none')+
	ggtitle('LR - NFD'),
draw_traffic_flow(SR_traffic[condition == 'DET'], lineend = 'round',
		  add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round',
		  		    add = geom_foodpatches(fill = 'mediumpurple')))+
	scale_size_continuous(range = c(1, 5)) + 
	scico::scale_color_scico() + theme(legend.position = 'none')+
	ggtitle('SR - DET'),
draw_traffic_flow(SR_traffic[condition == 'NFD'], lineend = 'round',
		  add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round'))+
	scale_size_continuous(range = c(1, 5)) + 
	scico::scale_color_scico() + theme(legend.position = 'none')+
	ggtitle('SR - NFD')
)



draw_hexagons(edges = LR_traffic, z = rank(LR_traffic[['N']]), linewidth = 3, lineend = 'round',
	      add = draw_hexagons(linewidth = 3, color = 'white', lineend = 'round', 
	      		    add = draw_hexagons(linewidth = 4, lineend = 'round')))+
	aes(size = N)

draw_hexagons(add = add)
add = draw_hexagons(LR_traffic, z = rank(LR_traffic[['N']]), lineend = 'round') + aes(size = rank(N))


draw_hexagons(add = draw_hexagons(LR_traffic, z = rank(LR_traffic[['N']]),
				  lineend = 'round') + aes(size = rank(N)))
draw_hexagons
draw_hexagons(SR_traffic, z = rank(SR_traffic[['N']]), linewidth = 3, lineend = 'round',
	      add = draw_hexagons(linewidth = 3, color = 'white', lineend = 'round', 
	      		    add = draw_hexagons(linewidth = 4, lineend = 'round')))



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
				    yend = yend, size = z**1.2), color = 'black',
				show.legend = FALSE, ...) +
		geom_segment(data = edges, 
			     aes(x = x, xend = xend, y = y,
			         yend = yend, color = z, size = z), ...)
		
	pl + xlab('') + ylab('') + theme(axis.text = element_blank(), axis.ticks = element_blank(),
					 axis.line = element_blank())
}

anchored_trajectories <- function(edges = NULL, origin = 634L, add = NULL, ...){
	
	if(is.null(edges)){
		edges <- compute_edges(hex[hex$y > 1000, ])
	}
	
	if(!is.null(add)){
		pl <- add
	} else {
		pl <- ggplot()
	}
	
	x0y0 <- hex[hex$node == origin, c('x', 'y')]
	edges[['x']] <- edges[['x']] - x0y0[['x']]
	edges[['xend']] <- edges[['xend']] - x0y0[['x']]
	edges[['y']] <- edges[['y']] - x0y0[['y']]
	edges[['yend']] <- edges[['yend']] - x0y0[['y']]
	
	pl <- pl + geom_segment(data = edges, 
				aes(x = x, xend = xend, y = y,
				    yend = yend, size = z**1.2), color = 'black',
				show.legend = FALSE, ...) +
		geom_segment(data = edges, 
			     aes(x = x, xend = xend, y = y,
			         yend = yend, color = z, size = z), ...)
	
	pl + xlab('') + ylab('') + theme(axis.text = element_blank(), axis.ticks = element_blank(),
					 axis.line = element_blank())
}

merge(edges, inds_det)
anchored_trajectories(merge(edges, inds_det))

traffic_flow <- nf_scouts[, .(z = .N), by = c('.id','o','d')]
LR_traffic <- merge(traffic_flow[.id == 'LR'], edges)
SR_traffic <- merge(traffic_flow[.id == 'SR'], edges)
draw_hexagons(LR_traffic, z = rank(LR_traffic[['N']]), linewidth = 3, lineend = 'round',
	      add = draw_hexagons(linewidth = 3, color = 'white', lineend = 'round', 
	      		    add = draw_hexagons(linewidth = 4, lineend = 'round')))
draw_hexagons(SR_traffic, z = rank(SR_traffic[['N']]), linewidth = 3, lineend = 'round',
	      add = draw_hexagons(linewidth = 3, color = 'white', lineend = 'round', 
	      		    add = draw_hexagons(linewidth = 4, lineend = 'round')))+
	scico::scale_color_scico()




############### ANCHORED TRAJECTORIES ###############
data_test <- setDT(det[[6]]@data)[Frame <= maxt]
ind_test <- data_test[N_ind == 3]
xy_traj <- TrajFromCoords(ind_test[, c('Xmm', 'Ymm', 'Frame')], fps = 2, spatialUnits = 'mm')
xy0 <- as.numeric(hex[hex$node == 634, c('x', 'y')]) - as.numeric(xy_traj[1, 1:2])
xy0_1 <- as.numeric(hex[hex$node == 634, c('x', 'y')]) - as.numeric(data_test[N_ind == 2][1, c('Xmm', 'Ymm')])
draw_hexagons() + geom_point(data = data_test[N_ind == 3], aes(Xmm, Ymm), size = 2)+ 
	geom_point(data = data_test[N_ind == 3][, x := Xmm + xy0[1]][, y := Ymm + xy0[2]], 
		   aes(x, y), size = 2, color = 'mediumpurple')+
	geom_point(data = data_test[N_ind == 2][, x := Xmm + xy0_1[1]][, y := Ymm + xy0_1[2]], 
		   aes(x, y), size = 2, color = 'red')+
	geom_point(data = hex[hex$node == 634, c('x', 'y')], aes(x, y), size = 5, shape = 13)



ods <- get_segment(ind_test[, c('Xmm', 'Ymm')])
ind_test[['o']] <- ods[, 1]
ind_test[['d']] <- ods[, 2]

m1 <- merge(ind_test, compute_edges(hex))
m2 <- m1
x0y0 <- as.numeric(hex[hex$node == 634, c('x', 'y')]) - as.numeric(m2[which.min(Frame), c('x','y')])
m2[['x']] <- m2[['x']] + x0y0[1]
m2[['xend']] <- m2[['xend']] + x0y0[1]
m2[['y']] <- m2[['y']] + x0y0[2]
m2[['yend']] <- m2[['yend']] + x0y0[2]
draw_hexagons(m2, linewidth = 2, color = 'mediumpurple', 
	      add = draw_hexagons(m1, color = 'gold3', linewidth = 2, add = draw_hexagons(compute_edges(hex))))+
	geom_point(data = hex[hex$node == 634, ], aes(x, y), size = 4, shape = 13)


anchored_destiny_det <- rbindlist(lapply(det, function(p){
	x00 <- as.numeric(hex[hex$node == 634, c('x', 'y')])
	e <- compute_edges(hex, unique = FALSE)
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

	result_LR <- c()
	result_SR <- c()
	ind_counter <- 0
	for(i in unique(data_LR[['N_ind']])){
		x0y0 <- x00 - as.numeric(hex[hex$node == data_LR[N_ind == i][order(Frame)][1, o], c('x', 'y')])
		dta_merged <- merge(e, data_LR[N_ind == i, c('o', 'd')])

		dta_merged[['x']] <- dta_merged[['x']] + x0y0[1]
		dta_merged[['xend']] <- dta_merged[['xend']] + x0y0[1]
		dta_merged[['y']] <- dta_merged[['y']] + x0y0[2]
		dta_merged[['yend']] <- dta_merged[['yend']] + x0y0[2]
		dta_merged[['N_ind']] <- ind_counter
		result_LR <- rbind(result_LR, dta_merged)
		ind_counter <- ind_counter + 1
	}
	
	for(i in unique(data_SR[['N_ind']])){
		x0y0 <- x00 - as.numeric(hex[hex$node == data_SR[N_ind == i][order(Frame)][1, o], c('x', 'y')])
		dta_merged <- merge(e, data_SR[N_ind == i, c('o', 'd')])
		dta_merged[['x']] <- dta_merged[['x']] + x0y0[1]
		dta_merged[['xend']] <- dta_merged[['xend']] + x0y0[1]
		dta_merged[['y']] <- dta_merged[['y']] + x0y0[2]
		dta_merged[['yend']] <- dta_merged[['yend']] + x0y0[2]
		dta_merged[['N_ind']] <- ind_counter
		result_SR <- rbind(result_SR, dta_merged)
		ind_counter <- ind_counter + 1
	}
	rbindlist(list(LR = result_LR, SR = result_SR), idcol = TRUE)
}), idcol = TRUE)


anchored_destiny_nf <- rbindlist(lapply(nf[-c(1, 2)], function(p){
	
	x00 <- as.numeric(hex[hex$node == 634, c('x', 'y')])
	data <- setDT(p@data)[Frame <= maxt]
	e <- compute_edges(refcoords = hex, unique = FALSE)
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
	
	result_LR <- c()
	result_SR <- c()
	ind_counter <- 0
	
	for(i in unique(data_LR[['N_ind']])){
		x0y0 <- x00 - as.numeric(hex[hex$node == data_LR[N_ind == i][order(Frame)][1, o], c('x', 'y')])
		dta_merged <- merge(e, data_LR[N_ind == i, c('o', 'd')])
		dta_merged[['x']] <- dta_merged[['x']] + x0y0[1]
		dta_merged[['xend']] <- dta_merged[['xend']] + x0y0[1]
		dta_merged[['y']] <- dta_merged[['y']] + x0y0[2]
		dta_merged[['yend']] <- dta_merged[['yend']] + x0y0[2]
		dta_merged[['N_ind']] <- ind_counter
		result_LR <- rbind(result_LR, dta_merged)
		ind_counter <- ind_counter + 1
	}
	for(i in unique(data_SR[['N_ind']])){
		x0y0 <- x00 - as.numeric(hex[hex$node == data_SR[N_ind == i][order(Frame)][1, o], c('x', 'y')])
		dta_merged <- merge(e, data_SR[N_ind == i, c('o', 'd')])
		dta_merged[['x']] <- dta_merged[['x']] + x0y0[1]
		dta_merged[['xend']] <- dta_merged[['xend']] + x0y0[1]
		dta_merged[['y']] <- dta_merged[['y']] + x0y0[2]
		dta_merged[['yend']] <- dta_merged[['yend']] + x0y0[2]
		dta_merged[['N_ind']] <- ind_counter
		result_SR <- rbind(result_SR, dta_merged)
		ind_counter <- ind_counter + 1
	}
	rbindlist(list(LR = result_LR, SR = result_SR), idcol = TRUE)
}), idcol = TRUE)



draw_hexagons(anchored_destiny_nf[, 2:ncol(anchored_destiny_nf)][.id == 'LR'])
draw_hexagons(anchored_destiny_nf[, 2:ncol(anchored_destiny_nf)][.id == 'SR'])


draw_hexagons(anchored_destiny_det[, 2:ncol(anchored_destiny_det)][.id == 'LR'])
draw_hexagons(anchored_destiny_det[, 2:ncol(anchored_destiny_det)][.id == 'SR'])
