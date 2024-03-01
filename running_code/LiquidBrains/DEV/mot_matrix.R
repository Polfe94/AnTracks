source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
load('~/research/gits/AnTracks/data/nf.RData')
load('~/research/gits/AnTracks/data/det.RData')
library(trajr)

get_mot_matrix <- function(data){
	inds <- unique(data[, N_ind])
	M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
	
	for(j in seq_along(inds)){
		df_ind <- data[N_ind == inds[j], node]
		nds_ind <- sapply(2:length(df_ind), function(i){
			if(df_ind[i] != df_ind[i-1]){
				return(df_ind[i])
			}
		})
		nds_ind <- append(df_ind[1], nds_ind)
		nds_ind <- nds_ind[lapply(nds_ind, length) > 0]
		if(length(nds_ind) > 2){
			dirs <- vapply(3:length(nds_ind), function(i){
				x <- list(hex[hex$node == nds_ind[[i-2]], c('x', 'y')],
					  hex[hex$node == nds_ind[[i-1]], c('x', 'y')],
					  hex[hex$node == nds_ind[[i]], c('x', 'y')])
				get_direction(x)
			}, numeric(1))
			result <- na.omit(dirs)
			if(length(result) > 1){
				M_ind <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
				for(ii in 2:length(result)){
					M[as.character(result[ii-1]), as.character(result[ii])] <-
						M[as.character(result[ii-1]), as.character(result[ii])]+1
				}
				
				
				M <- M + M_ind
			}
		}
	}
	M
}

maxd <- 500 # 650
maxt <- 3000
min_pos <- 3


LR_det <- lapply(seq_along(det), function(i){
	
	data <- setDT(det[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]

	data[, tracklet := rleid(N_ind, Crossings)]
	data[['id']] <- as.character(gsub(' ', '', apply(data[, c('N_ind', 'tracklet')], 1, paste0, sep = '', collapse = '.')))
	ids <- data[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, as.numeric(id)]
	data <- data[id %in% ids]

	get_mot_matrix(data)
	
})
apply(Reduce(`+`, LR_det) / do.call('sum', LR_det), 2, function(i) i / sum(i))

SR_det <- lapply(seq_along(det), function(i){
	
	data <- setDT(det[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[!N_ind %in% inds]
	
	data[, tracklet := rleid(N_ind, Crossings)]
	data[['id']] <- as.character(gsub(' ', '', apply(data[, c('N_ind', 'tracklet')], 1, paste0, sep = '', collapse = '.')))
	ids <- data[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, as.numeric(id)]
	data <- data[id %in% ids]
	
	get_mot_matrix(data)
	
})
apply(Reduce(`+`, SR_det) / do.call('sum', SR_det), 2, function(i) i / sum(i))

# data[, tracklet := rleid(N_ind, Crossings)]
# data[['id']] <- as.character(gsub(' ', '', apply(data[, c('N_ind', 'tracklet')], 1, paste0, sep = '', collapse = '.')))
# ids <- data[, .N, by = c('id', 'node')][, .N, by = 'id'][N > min_pos, as.numeric(id)]
# data[id %in% ids]
# draw_hexagons() + geom_point(data = data[id %in% ids], aes(Xmm, Ymm, color = id, size = id), alpha = 0.2)+
# 	scale_color_viridis_d()
# 
# draw_hexagons() + geom_point(data = data, aes(Xmm, Ymm, color = id, size = id), alpha = 0.2)+
# 	scale_color_viridis_d()
# 
# draw_hexagons() + geom_point(data = merge(hex, data), aes(x, y, color = Frame, size = id), alpha = 0.2)+
# 	scale_color_viridis_c()














LR_nf <- lapply(seq_along(nf[-c(1, 2)]), function(i){
	
	data <- setDT(nf[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	
	data[, tracklet := rleid(N_ind, Crossings)]
	data[['id']] <- as.character(gsub(' ', '', apply(data[, c('N_ind', 'tracklet')], 1, paste0, sep = '', collapse = '.')))
	ids <- data[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, as.numeric(id)]
	data <- data[id %in% ids]
	
	get_mot_matrix(data)
	
})
apply(Reduce(`+`, LR_nf) / do.call('sum', LR_nf), 2, function(i) i / sum(i))

SR_nf <- lapply(seq_along(nf[-c(1, 2)]), function(i){
	
	data <- setDT(nf[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[!N_ind %in% inds]
	
	data[, tracklet := rleid(N_ind, Crossings)]
	data[['id']] <- as.character(gsub(' ', '', apply(data[, c('N_ind', 'tracklet')], 1, paste0, sep = '', collapse = '.')))
	ids <- data[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, as.numeric(id)]
	data <- data[id %in% ids]
	
	get_mot_matrix(data)
	
})
apply(Reduce(`+`, SR_nf) / do.call('sum', SR_nf), 2, function(i) i / sum(i))


LR_sto <- lapply(seq_along(sto), function(i){
	
	data <- setDT(sto[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	
	data[, tracklet := rleid(N_ind, Crossings)]
	data[['id']] <- as.character(gsub(' ', '', apply(data[, c('N_ind', 'tracklet')], 1, paste0, sep = '', collapse = '.')))
	ids <- data[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, as.numeric(id)]
	data <- data[id %in% ids]
	
	get_mot_matrix(data)
	
})
apply(Reduce(`+`, LR_sto) / do.call('sum', LR_sto), 2, function(i) i / sum(i))

SR_sto <- lapply(seq_along(sto), function(i){
	
	data <- setDT(sto[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[!N_ind %in% inds]
	
	data[, tracklet := rleid(N_ind, Crossings)]
	data[['id']] <- as.character(gsub(' ', '', apply(data[, c('N_ind', 'tracklet')], 1, paste0, sep = '', collapse = '.')))
	ids <- data[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, as.numeric(id)]
	data <- data[id %in% ids]
	get_mot_matrix(data)
	
})
apply(Reduce(`+`, SR_sto) / do.call('sum', SR_sto), 2, function(i) i / sum(i))


####### straightness and other metrics ....####

LR_det <- rbindlist(lapply(seq_along(det), function(i){
	
	data <- setDT(det[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	
	data[, tracklet := rleid(N_ind, Crossings)]
	data[['id']] <- as.character(gsub(' ', '', apply(data[, c('N_ind', 'tracklet')], 1, paste0, sep = '', collapse = '.')))
	ids <- data[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, id]
	data <- data[id %in% ids]
	mdata <- data.table(merge(hex, data, by = 'node'))
	
	df <- data.frame(t(vapply(seq_along(ids), function(ii){
		xy <- TrajFromCoords(mdata[id == ids[ii], c('x', 'y', 'Frame')], fps = 2)

		c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
		  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
		  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
		
	}, numeric(6))))
	colnames(df)<- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')
	df[['expl']] <- 'LR'
	df
	
}))
SR_det <- rbindlist(lapply(seq_along(det), function(i){
	
	data <- setDT(det[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[!N_ind %in% inds]
	data[['tracklet']] <- data[, rleid(N_ind, Crossings)]
	data[['id']] <- as.character(gsub(' ', '', apply(data[, c('N_ind', 'tracklet')], 1, paste0, sep = '', collapse = '.')))
	ids <- data[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, id]
	data <- data[id %in% ids]
	mdata <- data.table(merge(hex, data, by = 'node'))
	
	df <- data.frame(t(vapply(seq_along(ids), function(ii){
		xy <- TrajFromCoords(mdata[id == ids[ii], c('x', 'y', 'Frame')], fps = 2)
		c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
		  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
		  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)

		
	}, numeric(6))))
	colnames(df)<- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')
	df[['expl']] <- 'SR'
	df
	
}))

LR_nf <- rbindlist(lapply(seq_along(nf[-c(1, 2)]), function(i){
	
	data <- setDT(nf[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	
	data[, tracklet := rleid(N_ind, Crossings)]
	data[['id']] <- as.character(gsub(' ', '', apply(data[, c('N_ind', 'tracklet')], 1, paste0, sep = '', collapse = '.')))
	ids <- data[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, id]
	data <- data[id %in% ids]
	mdata <- data.table(merge(hex, data, by = 'node'))
	
	df <- data.frame(t(vapply(seq_along(ids), function(ii){
		xy <- TrajFromCoords(mdata[id == ids[ii], c('x', 'y', 'Frame')], fps = 2)
		c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
		  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
		  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
		
	}, numeric(6))))
	colnames(df)<- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')

	if(nrow(df)){
		df[['expl']] <- 'LR'
		df
	} else {NULL}
	
}))
SR_nf <- rbindlist(lapply(seq_along(nf[-c(1, 2)]), function(i){
	
	data <- setDT(nf[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[!N_ind %in% inds]
	data[['tracklet']] <- data[, rleid(N_ind, Crossings)]
	data[['id']] <- as.character(gsub(' ', '', apply(data[, c('N_ind', 'tracklet')], 1, paste0, sep = '', collapse = '.')))
	ids <- data[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, id]
	data <- data[id %in% ids]
	mdata <- data.table(merge(hex, data, by = 'node'))
	
	df <- data.frame(t(vapply(seq_along(ids), function(ii){
		xy <- TrajFromCoords(mdata[id == ids[ii], c('x', 'y', 'Frame')], fps = 2)
		c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
		  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
		  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
		
		
	}, numeric(6))))
	colnames(df)<- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')
	df[['expl']] <- 'SR'
	df
	
}))

LR_sto <- rbindlist(lapply(seq_along(sto[-c(1, 2)]), function(i){
	
	data <- setDT(sto[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	
	data[, tracklet := rleid(N_ind, Crossings)]
	data[['id']] <- as.character(gsub(' ', '', apply(data[, c('N_ind', 'tracklet')], 1, paste0, sep = '', collapse = '.')))
	ids <- data[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, id]
	data <- data[id %in% ids]
	mdata <- data.table(merge(hex, data, by = 'node'))
	
	df <- data.frame(t(vapply(seq_along(ids), function(ii){
		xy <- TrajFromCoords(mdata[id == ids[ii], c('x', 'y', 'Frame')], fps = 2)
		c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
		  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
		  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
		
	}, numeric(6))))
	colnames(df)<- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')
	
	if(nrow(df)){
		df[['expl']] <- 'LR'
		df
	} else {NULL}
	
}))
SR_sto <- rbindlist(lapply(seq_along(sto), function(i){
	
	data <- setDT(sto[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[!N_ind %in% inds]
	data[['tracklet']] <- data[, rleid(N_ind, Crossings)]
	data[['id']] <- as.character(gsub(' ', '', apply(data[, c('N_ind', 'tracklet')], 1, paste0, sep = '', collapse = '.')))
	ids <- data[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, id]
	data <- data[id %in% ids]
	mdata <- data.table(merge(hex, data, by = 'node'))
	
	df <- data.frame(t(vapply(seq_along(ids), function(ii){
		xy <- TrajFromCoords(mdata[id == ids[ii], c('x', 'y', 'Frame')], fps = 2)
		c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
		  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
		  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
		
		
	}, numeric(6))))
	colnames(df)<- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')
	df[['expl']] <- 'SR'
	df
	
}))


ggplot(data = rbindlist(list(LR_det = LR_det[,.(straightness = mean(Straightness)), 
					     by = .(traj_length = round(Distance/100)*100)],
			     SR_det = SR_det[,.(straightness = mean(Straightness)), 
			     		by = .(traj_length = round(Distance/100)*100)],
			     LR_nf = LR_nf[,.(straightness = mean(Straightness)), 
			     		by = .(traj_length = round(Distance/100)*100)],
			     SR_nf = SR_nf[,.(straightness = mean(Straightness)), 
			     		by = .(traj_length = round(Distance/100)*100)]),
			idcol = TRUE)[order(traj_length)],
       aes(traj_length/ 10, straightness, color = factor(.id)), alpha = 0.5) + 
	# aes(log(traj_length/ 10), straightness, color = factor(.id)), alpha = 0.5) + 
	geom_point(size = 4, shape = 21) +
	geom_path()+
	scale_y_continuous('<Net-to-Gross ratio>', limits = c(0, 1), breaks = seq(0, 1, 0.2))+
	scale_x_continuous('Travelled distance (cm)')+
	scale_color_manual('', values = c('mediumpurple', 'gold3', 'brown4', 'grey'), 
			   labels = c('LR DET', 'LR NFD', 'SR DET', 'SR NFD'))+
	theme(legend.title = element_blank(), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.position = c(3/4, 3/4))


ggplot(data = rbindlist(list(DET = rbind(LR_det[Mean_v > 0], 
				     SR_det[Mean_v > 0]),
			     NFD = rbind(LR_nf[Mean_v > 0], 
			     	    SR_nf[Mean_v > 0]),
			     STO = rbind(LR_sto[Mean_v > 0], 
			     	    SR_sto[Mean_v > 0])), idcol = TRUE),
	       aes(expl, Straightness, fill = .id)) + 
	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
	geom_jitter(aes(color = .id), 
		    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
		    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition',
			  values = c('mediumpurple', 'gold3', 'brown4'),
			  labels = c('DET', 'NFD', 'STO'))+
	scale_color_manual(values = rep('black', 3))+
	theme(legend.title = element_blank())+
	ylab('Net-to-Gross ratio')+
	scale_x_discrete('', labels = c('Long-range scouts', 'Short-range scouts'))
ggplot(data = rbindlist(list(DET = rbind(LR_det[Mean_v > 0], 
					 SR_det[Mean_v > 0]),
					 NFD = rbind(LR_nf[Mean_v > 0], 
					 	    SR_nf[Mean_v > 0]),
					 STO = rbind(LR_sto[Mean_v > 0], 
					 	    SR_sto[Mean_v > 0])), idcol = TRUE),
       aes(expl, Mean_v, fill = .id)) + 
	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
	geom_jitter(aes(color = .id), 
		    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
		    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition',
			  values = c('mediumpurple', 'gold3', 'brown4'),
			  labels = c('DET', 'NFD', 'STO'))+
	scale_color_manual(values = rep('black', 3))+
	theme(legend.title = element_blank())+
	ylab('<Velocity> (mm/s)')+
	scale_x_discrete('', labels = c('Long-range scouts', 'Short-range scouts'))

ggplot(data = rbindlist(list(DET = rbind(LR_det[Mean_v > 0], 
					 SR_det[Mean_v > 0]),
					 NFD = rbind(LR_nf[Mean_v > 0], 
					 	    SR_nf[Mean_v > 0]),
					 STO = rbind(LR_sto[Mean_v > 0], 
					 	    SR_sto[Mean_v > 0])), idcol = TRUE),
       aes(expl, Mean_acc, fill = .id)) + 
	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
	geom_jitter(aes(color = .id), 
		    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
		    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition',
			  values = c('mediumpurple', 'gold3', 'brown4'),
			  labels = c('DET', 'NFD', 'STO'))+
	scale_color_manual(values = rep('black', 3))+
	theme(legend.title = element_blank())+
	ylab(latex2exp::TeX('<Acceleration> (mm / s$^{2}$)'))+
	scale_x_discrete('', labels = c('Long-range scouts', 'Short-range scouts'))
ggplot(data = rbindlist(list(DET = rbind(LR_det[Mean_v > 0], 
					 SR_det[Mean_v > 0]),
					 NFD = rbind(LR_nf[Mean_v > 0], 
					 	    SR_nf[Mean_v > 0]),
					 STO = rbind(LR_sto[Mean_v > 0], 
					 	    SR_sto[Mean_v > 0])), idcol = TRUE),
       aes(expl, Diffusion, fill = .id)) + 
	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
	geom_jitter(aes(color = .id), 
		    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
		    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition',
			  values = c('mediumpurple', 'gold3', 'brown4'),
			  labels = c('DET', 'NFD', 'STO'))+
	scale_color_manual(values = rep('black', 3))+
	theme(legend.title = element_blank())+
	ylab('Diffusion')+
	scale_x_discrete('', labels = c('Long-range scouts', 'Short-range scouts'))
ggplot(data = rbindlist(list(DET = rbind(LR_det[Mean_v > 0], 
					 SR_det[Mean_v > 0]),
					 NFD = rbind(LR_nf[Mean_v > 0], 
					 	    SR_nf[Mean_v > 0]),
					 STO = rbind(LR_sto[Mean_v > 0], 
					 	    SR_sto[Mean_v > 0])), idcol = TRUE),
       aes(expl, Distance, fill = .id)) + 
	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
	geom_jitter(aes(color = .id), 
		    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
		    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition',
			  values = c('mediumpurple', 'gold3', 'brown4'),
			  labels = c('DET', 'NFD', 'STO'))+
	scale_color_manual(values = rep('black', 3))+
	theme(legend.title = element_blank())+
	ylab('Travel distance (mm)')+
	scale_x_discrete('', labels = c('Long-range scouts', 'Short-range scouts'))
ggplot(data = rbindlist(list(DET = rbind(LR_det[Mean_v > 0], 
					 SR_det[Mean_v > 0]),
					 NFD = rbind(LR_nf[Mean_v > 0], 
					 	    SR_nf[Mean_v > 0]),
					 STO = rbind(LR_sto[Mean_v > 0], 
					 	    SR_sto[Mean_v > 0])), idcol = TRUE),
       aes(expl, Time, fill = .id)) + 
	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
	geom_jitter(aes(color = .id), 
		    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
		    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition',
			  values = c('mediumpurple', 'gold3', 'brown4'),
			  labels = c('DET', 'NFD', 'STO'))+
	scale_color_manual(values = rep('black', 3))+
	theme(legend.title = element_blank())+
	ylab('Time spent exploring (min)')+
	scale_x_discrete('', labels = c('Long-range scouts', 'Short-range scouts'))


# 
# ggplot(data = rbindlist(list(LR_det[Mean_v > 0], 
# 			     SR_det[Mean_v > 0])),
#        aes(expl, Mean_v)) + 
# 	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
# 	geom_jitter(aes(color = .id), 
# 		    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
# 		    size = 3, show.legend = FALSE, alpha = 0.25)+
# 	theme(legend.title = element_blank())+
# 	xlab('') + ylab('< Velocity >')
# 
# ggplot(data = rbindlist(list(LR_det[Mean_v > 0], 
# 			     SR_det[Mean_v > 0])),
#        aes(expl, Diffusion)) + 
# 	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
# 	geom_jitter(aes(color = .id), 
# 		    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
# 		    size = 3, show.legend = FALSE, alpha = 0.25)+
# 	theme(legend.title = element_blank())+
# 	xlab('') + ylab('Diffusion')
# 
# ggplot(data = rbindlist(list(LR_det[Mean_v > 0], 
# 			     SR_det[Mean_v > 0])),
#        aes(expl, Distance)) + 
# 	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
# 	geom_jitter(aes(color = .id), 
# 		    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
# 		    size = 3, show.legend = FALSE, alpha = 0.25)+
# 	theme(legend.title = element_blank())+
# 	xlab('') + ylab('Travel distance')
# 
# ggplot(data = rbindlist(list(LR_det[Mean_v > 0], 
# 			     SR_det[Mean_v > 0])),
#        aes(expl, Mean_acc)) + 
# 	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
# 	geom_jitter(aes(color = .id), 
# 		    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
# 		    size = 3, show.legend = FALSE, alpha = 0.25)+
# 	theme(legend.title = element_blank())+
# 	xlab('') + ylab('< Acceleration >')
# 
# ggplot(data = rbindlist(list(LR_det[Mean_v > 0], 
# 			     SR_det[Mean_v > 0])),
#        aes(expl, Time)) + 
# 	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
# 	geom_jitter(aes(color = .id), 
# 		    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
# 		    size = 3, show.legend = FALSE, alpha = 0.25)+
# 	theme(legend.title = element_blank())+
# 	xlab('') + ylab('Time in arena (min)')
# 
# 
# 


det_scouts <- lapply(seq_along(det), function(i){
	
	data <- setDT(det[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data[, tracklet := rleid(N_ind, Crossings)]
	data[['id']] <- as.character(gsub(' ', '', apply(data[, c('N_ind', 'tracklet')], 1, paste0, sep = '', collapse = '.')))
	data[['d']] <- get_segment(data[, c('Xmm', 'Ymm')])[, 2]
	data[['o']] <- data[['node']]
	data_LR <- data[N_ind %in% inds][, scout := 'LR']
	data_SR <- data[!N_ind %in% inds][, scout := 'SR']
	ids_LR <- data_LR[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, as.numeric(id)]
	ids_SR <- data_SR[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, as.numeric(id)]
	rbindlist(list(LR = data_LR[id %in% ids_LR], SR = data_SR[id %in% ids_SR]))
})
det_density <- rbindlist(det_scouts)[, .(N = .N), by = c('scout','o','d')][, condition := 'DET'][, z := rank(N)]
density_LR <- data.table(apply(do.call('rbind', 
				       strsplit(rle(apply(rbindlist(det_scouts)[scout == 'LR', c('o','d')], 1,
				       		   paste, collapse = '_', sep = '_'))$values, '_')), 2, as.numeric))
colnames(density_LR) <- c('o', 'd')
data_flow <- merge(edges, density_LR[, N := .N, by = c('o', 'd')][, z := rank(N)])
draw_traffic_flow(data_flow, lineend = 'round',
		  add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round',
		  		    add = geom_foodpatches(fill = 'mediumpurple')))+
	scale_size_continuous(range = c(1, 5)) + 
	scico::scale_color_scico() + theme(legend.position = 'none')+
	ggtitle('LR - DET')

density_SR <- data.table(apply(do.call('rbind', 
				       strsplit(rle(apply(rbindlist(det_scouts)[scout == 'SR', c('o','d')], 1,
				       		   paste, collapse = '_', sep = '_'))$values, '_')), 2, as.numeric))
colnames(density_SR) <- c('o', 'd')
data_flow <- merge(edges, density_SR[, N := .N, by = c('o', 'd')][, z := rank(N)])
draw_traffic_flow(data_flow, lineend = 'round',
		  add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round',
		  		    add = geom_foodpatches(fill = 'mediumpurple')))+
	scale_size_continuous(range = c(1, 5)) + 
	scico::scale_color_scico() + theme(legend.position = 'none')+
	ggtitle('SR - DET')




nf_scouts <- lapply(seq_along(nf[-c(1, 2)]), function(i){
	
	data <- setDT(nf[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data[, tracklet := rleid(N_ind, Crossings)]
	data[['id']] <- as.character(gsub(' ', '', apply(data[, c('N_ind', 'tracklet')], 1, paste0, sep = '', collapse = '.')))
	data[['d']] <- get_segment(data[, c('Xmm', 'Ymm')])[, 2]
	data[['o']] <- data[['node']]
	data_LR <- data[N_ind %in% inds][, scout := 'LR']
	data_SR <- data[!N_ind %in% inds][, scout := 'SR']
	ids_LR <- data_LR[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, as.numeric(id)]
	ids_SR <- data_SR[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, as.numeric(id)]
	rbindlist(list(LR = data_LR[id %in% ids_LR], SR = data_SR[id %in% ids_SR]))
})
nf_density <- rbindlist(nf_scouts)[, .(N = .N), by = c('scout','o','d')][, condition := 'nf'][, z := rank(N)]
density_LR <- data.table(apply(do.call('rbind', 
				       strsplit(rle(apply(rbindlist(nf_scouts)[scout == 'LR', c('o','d')], 1,
				       		   paste, collapse = '_', sep = '_'))$values, '_')), 2, as.numeric))
colnames(density_LR) <- c('o', 'd')
data_flow <- merge(edges, density_LR[, N := .N, by = c('o', 'd')][, z := rank(N)])
draw_traffic_flow(data_flow, lineend = 'round',
		  add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round',
		  		    add = geom_foodpatches(fill = 'mediumpurple')))+
	scale_size_continuous(range = c(1, 5)) + 
	scico::scale_color_scico() + theme(legend.position = 'none')+
	ggtitle('LR - nf')

density_SR <- data.table(apply(do.call('rbind', 
				       strsplit(rle(apply(rbindlist(nf_scouts)[scout == 'SR', c('o','d')], 1,
				       		   paste, collapse = '_', sep = '_'))$values, '_')), 2, as.numeric))
colnames(density_SR) <- c('o', 'd')
data_flow <- merge(edges, density_SR[, N := .N, by = c('o', 'd')][, z := rank(N)])
draw_traffic_flow(data_flow, lineend = 'round',
		  add = draw_hexagons(linewidth = 1, color = 'grey80', lineend = 'round',
		  		    add = geom_foodpatches(fill = 'mediumpurple')))+
	scale_size_continuous(range = c(1, 5)) + 
	scico::scale_color_scico() + theme(legend.position = 'none')+
	ggtitle('SR - nf')














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




#####################################################
data <- setDT(nf[[i]]@data)[Frame <= maxt]
dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
set(data, j = 'd', value = dmatrix)
data[, tracklet := rleid(N_ind, Crossings)]
data[['id']] <- as.character(gsub(' ', '', apply(data[, c('N_ind', 'tracklet')], 1, paste0, sep = '', collapse = '.')))
data[['d']] <- get_segment(data[, c('Xmm', 'Ymm')])[, 2]
data[['o']] <- data[['node']]
data_LR <- data[N_ind %in% inds][, scout := 'LR']
data_SR <- data[!N_ind %in% inds][, scout := 'SR']
ids_LR <- data_LR[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, as.numeric(id)]
ids_SR <- data_SR[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, as.numeric(id)]
LR = data_LR[id %in% ids_LR]
SR = data_SR[id %in% ids_SR]

maxt <- 3000
scouts__DET <- lapply(det, function(p){
	data <- setDT(p@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data[, tracklet := rleid(N_ind, Crossings)]
	data[['id']] <- as.character(gsub(' ', '', apply(data[, c('N_ind', 'tracklet')], 1, paste0, sep = '', collapse = '.')))
	data[['d']] <- get_segment(data[, c('Xmm', 'Ymm')])[, 2]
	data[['o']] <- data[['node']]
	data_LR <- data[N_ind %in% inds][, scout := 'LR']
	data_SR <- data[!N_ind %in% inds][, scout := 'SR']
	ids_LR <- data_LR[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, as.numeric(id)]
	ids_SR <- data_SR[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, as.numeric(id)]
	LR = data_LR[id %in% ids_LR]
	SR = data_SR[id %in% ids_SR]
	LR_list <- rbindlist(lapply(seq_along(ids_LR), function(i){
		xy <- TrajFromCoords(LR[id == ids_LR[i], c('Xmm', 'Ymm', 'Frame')], fps = 2)
		Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
		maxL <- max(Ls)
		if(!is.na(maxL)){
			if(maxL >= 100){
				idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
			} else {
				idxs <- nrow(xy)
			}
			straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
			data.frame(l = Ls[idxs], s = straightness)
		} else {NULL}

	}), idcol = TRUE)
	
	SR_list <- rbindlist(lapply(seq_along(ids_SR), function(i){
		xy <- TrajFromCoords(SR[id == ids_SR[i], c('Xmm', 'Ymm', 'Frame')], fps = 2)
		Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
		maxL <- max(Ls)
		if(!is.na(maxL)){
			if(maxL >= 100){
				idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
			} else {
				idxs <- nrow(xy)
			}
			straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
			data.frame(l = Ls[idxs], s = straightness)
		} else {NULL}
		
	}), idcol = TRUE)
	if(nrow(LR_list)){
		LR_list[['scout']] <- 'LR'
	}
	if(nrow(SR_list)){
		SR_list[['scout']] <- 'SR'
	}
	
	rbind(LR_list, SR_list)
})
names(scouts__DET) <- paste0('Det_', seq_along(det))
scouts__DET <- rbindlist(scouts__DET, idcol = TRUE)
scouts__DET[['tag']] <- apply(scouts__DET[, 1:2], 1, paste, collapse = '_')

maxt <- 3000
scouts__NFD <- lapply(nf[-c(1, 2)], function(p){
	data <- setDT(p@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data[, tracklet := rleid(N_ind, Crossings)]
	data[['id']] <- as.character(gsub(' ', '', apply(data[, c('N_ind', 'tracklet')], 1, paste0, sep = '', collapse = '.')))
	data[['d']] <- get_segment(data[, c('Xmm', 'Ymm')])[, 2]
	data[['o']] <- data[['node']]
	data_LR <- data[N_ind %in% inds][, scout := 'LR']
	data_SR <- data[!N_ind %in% inds][, scout := 'SR']
	ids_LR <- data_LR[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, as.numeric(id)]
	ids_SR <- data_SR[, .N, by = c('id', 'node')][, .N, by = 'id'][N >= min_pos, as.numeric(id)]
	LR = data_LR[id %in% ids_LR]
	SR = data_SR[id %in% ids_SR]
	LR_list <- rbindlist(lapply(seq_along(ids_LR), function(i){
		xy <- TrajFromCoords(LR[id == ids_LR[i], c('Xmm', 'Ymm', 'Frame')], fps = 2)
		Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
		maxL <- max(Ls)
		if(!is.na(maxL)){
			if(maxL >= 100){
				idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
			} else {
				idxs <- nrow(xy)
			}
			straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
			data.frame(l = Ls[idxs], s = straightness)
		} else {NULL}
		
	}), idcol = TRUE)
	
	SR_list <- rbindlist(lapply(seq_along(ids_SR), function(i){
		xy <- TrajFromCoords(SR[id == ids_SR[i], c('Xmm', 'Ymm', 'Frame')], fps = 2)
		Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
		maxL <- max(Ls)
		if(!is.na(maxL)){
			if(maxL >= 100){
				idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
			} else {
				idxs <- nrow(xy)
			}
			straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
			data.frame(l = Ls[idxs], s = straightness)
		} else {NULL}
		
	}), idcol = TRUE)
	if(nrow(LR_list)){
		LR_list[['scout']] <- 'LR'
	}
	if(nrow(SR_list)){
		SR_list[['scout']] <- 'SR'
	}

	rbind(LR_list, SR_list)
})
names(scouts__NFD) <- paste0('nf_', seq_along(nf[-c(1,2)]))
scouts__NFD <- rbindlist(scouts__NFD, idcol = TRUE)
scouts__NFD[['tag']] <- apply(scouts__NFD[, 1:2], 1, paste, collapse = '_')

ggplot(data = scouts__NFD[,.(straightness = mean(s)), by = .(traj_length = round(l/100)*100)],
       aes(traj_length, straightness)) + geom_point()




ggplot(data = rbindlist(DET= scouts__DET, NFD =scouts__NFD, idcol = TRUE)[,.(straightness = mean(s)), 
					      by = .(traj_length = round(l/100)*100,
					             scout = scout,
					             type = .id)]),
       aes(traj_length/ 10, straightness, color = factor(type)), alpha = 0.5) + 
	geom_point(size = 4, shape = 21) +
	geom_smooth(formula = y ~ log(x), se = FALSE)+
	scale_y_continuous('<Net-to-Gross ratio>', limits = c(0, 1), breaks = seq(0, 1, 0.2))+
	scale_x_continuous('Travelled distance (cm)')+
	scale_color_manual('', values = c('mediumpurple', 'gold3'), 
			   labels = c('LR scouts', 'SR scouts'))+
	theme(legend.title = element_blank(), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.position = c(3/4, 3/4))


















SRs <- lapply(det, function(p){
	data <- setDT(p@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	set(data, j = 'd', value = as.numeric(dmatrix))
	SR <- unique(data[d > maxd, N_ind])
	SR <- unique(data[!N_ind %in% SR, N_ind])
	data_SR <- data[N_ind %in% SR]
	data_SR <- data[N_ind %in% SR]
	rbindlist(lapply(seq_along(SR), function(i){
		sbst <- data[N_ind == SR[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 3){
			Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
			maxL <- max(Ls)
			if(maxL >= 100){
				idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
			} else {
				idxs <- nrow(xy)
			}
			straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
			data.frame(l = Ls[idxs], s = straightness)
		}
		
	}), idcol = TRUE)})
names(SRs) <- paste0('Det_', seq_along(det))
SRs <- rbindlist(SRs, idcol = TRUE)
SRs[['tag']] <- apply(SRs[, 1:2], 1, paste, collapse = '_')

ggplot(data = SRs[,.(straightness = mean(s)), by = .(traj_length = round(l/100)*100)],
       aes(traj_length, straightness)) + geom_point()







ggplot(data = rbindlist(list(LR = LRs[,.(straightness = mean(s)), by = .(traj_length = round(l/100)*100)],
			     SR = SRs[,.(straightness = mean(s)), by = .(traj_length = round(l/100)*100)]),
			idcol = TRUE),
       aes(traj_length/ 10, straightness, color = factor(.id)), alpha = 0.5) + 
	geom_point(size = 4, shape = 21) +
	geom_smooth(formula = y ~ log(x), se = FALSE)+
	scale_y_continuous('<Net-to-Gross ratio>', limits = c(0, 1), breaks = seq(0, 1, 0.2))+
	scale_x_continuous('Travelled distance (cm)')+
	scale_color_manual('', values = c('mediumpurple', 'gold3'), 
			   labels = c('LR scouts', 'SR scouts'))+
	theme(legend.title = element_blank(), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.position = c(3/4, 3/4))






maxt <- 2400
LRs_nf <- lapply(nf, function(p){
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
names(LRs_nf) <- paste0('nf_', seq_along(nf))
LRs_nf <- rbindlist(LRs_nf, idcol = TRUE)
LRs_nf[['tag']] <- apply(LRs_nf[, 1:2], 1, paste, collapse = '_')


SRs_nf <- lapply(nf[-c(1, 2)], function(p){
	data <- setDT(p@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	set(data, j = 'd', value = as.numeric(dmatrix))
	SR <- unique(data[d > maxd, N_ind])
	SR <- unique(data[!N_ind %in% SR, N_ind])
	data_SR <- data[N_ind %in% SR]
	data_SR <- data[N_ind %in% SR]
	rbindlist(lapply(seq_along(SR), function(i){
		sbst <- data[N_ind == SR[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 3){
			Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
			maxL <- max(Ls)
			if(maxL >= 100){
				idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
			} else {
				idxs <- nrow(xy)
			}
			straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
			data.frame(l = Ls[idxs], s = straightness)
		}
		
	}), idcol = TRUE)})
names(SRs_nf) <- paste0('nf_', seq_along(nf[-c(1, 2)]))
SRs_nf <- rbindlist(SRs_nf, idcol = TRUE)
SRs_nf[['tag']] <- apply(SRs_nf[, 1:2], 1, paste, collapse = '_')




ggplot(data = rbindlist(list(LR = LRs[,.(straightness = mean(s)), by = .(traj_length = round(l/100)*100)],
			     SR = SRs[,.(straightness = mean(s)), by = .(traj_length = round(l/100)*100)],
			     LR_nf = LRs_nf[,.(straightness = mean(s)), by = .(traj_length = round(l/100)*100)],
			     SR_nf = SRs_nf[,.(straightness = mean(s)), by = .(traj_length = round(l/100)*100)]),
			idcol = TRUE),
       aes(traj_length/ 10, straightness, color = factor(.id)), alpha = 0.5) + 
	# aes(log(traj_length/ 10), straightness, color = factor(.id)), alpha = 0.5) + 
	geom_point(size = 4, shape = 21) +
	geom_smooth(formula = y ~ log(x), se = FALSE)+
	scale_y_continuous('<Net-to-Gross ratio>', limits = c(0, 1), breaks = seq(0, 1, 0.2))+
	scale_x_continuous('Travelled distance (cm)')+
	scale_color_manual('', values = c('mediumpurple', 'gold3', 'brown4', 'grey'), 
			   labels = c('LR DET', 'LR NFD', 'SR DET', 'SR NFD'))+
	theme(legend.title = element_blank(), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.position = c(3/4, 3/4))