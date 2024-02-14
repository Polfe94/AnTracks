source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
load('~/research/gits/AnTracks/data/nf.RData')
load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/sto.RData')



#### FUNCTIONS

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


################################# MUTUAL INFORMATION #################################
load("/home/polfer/research/gits/AnTracks/results/MI_ballistic_norm.RData")
load("/home/polfer/research/gits/AnTracks/results/MI_det.RData")
load("/home/polfer/research/gits/AnTracks/results/MI_default.RData")

mi_bal <- rbindlist(MI_ballistic)[!is.na(t)]
mi_det <- rbindlist(lapply(det_mi, function(i) data.frame(t = as.numeric(names(i)), x = unlist(i))))
mi_exp <-  rbindlist(lapply(mutual_info_result[lapply(mutual_info_result, length)> 1],
			    function(i) data.frame(t = as.numeric(names(i)), x = as.numeric(unlist(i)))))

MI <- rbindlist(list(mi_bal[, .(x = mean(x), cond = 'ballistic'), by = 't'], 
		     mi_det[, .(x = mean(x), cond = 'det'), by = 't'], 
		     mi_exp[, .(x = mean(x), cond = 'empiric_move'), by = 't']))

ggplot(data = MI, aes(t, x, color = cond)) + 
	geom_point(shape = 21, fill = NA, size = 4)+
	geom_line(size = 1) + ylab('<MI>') +
	scale_x_continuous('Time (min)', breaks = 120*seq(0, 150, 50), labels = seq(0, 150, 50))+
	theme(legend.title = element_blank())+
	scale_color_manual(values = c('mediumpurple', 'gold3', 'brown4'),
			   labels = c('Ballistic exploration', 'Experiments', 'Simulations'))

################################# DET EXPERIMENTS #################################

foods <- rbindlist(lapply(det, function(i){
	d <- rbindlist(i@food)
	idx <- which.min(d[['t']])
	d[idx]
}))[, node := get_node(x, y), by = seq_along(det)]

library(trajr)
exp_trajs <- vapply(seq_along(det), function(i){
	xy <- TrajFromCoords(ids[[i]][['data']][Frame <= foods[i, t], c('Xmm', 'Ymm', 'Frame')], fps = 2)
	TrajStraightness(xy)
}, numeric(1))

exploration_trajs <- do.call('rbind.data.frame', lapply(seq_along(det), function(i){
	data <- setDT(det[[i]]@data)[Frame <= foods[i, t]]
	inds <- unique(data[['N_ind']])
	t(vapply(seq_along(inds), function(ii){
		xy <- TrajFromCoords(data[N_ind == inds[ii], c('Xmm', 'Ymm', 'Frame')], fps = 2)
		c(TrajStraightness(xy), nrow(xy), TrajDistance(xy))
	}, numeric(3)))
}))
colnames(exploration_trajs) <- c('straightness', 'length', 'distance')

ggplot(data = exploration_trajs,aes('Explorers', straightness)) + geom_boxplot() +
	geom_jitter(width = 0.2, size = 3, alpha = 0.2)

ggplot(data = as.data.frame(exploration_trajs), aes(length, straightness))+geom_point()
ggplot(data = as.data.frame(exploration_trajs), aes('', distance))+geom_boxplot(outlier.shape = NA)+
	geom_jitter()
ggplot(data = as.data.frame(exploration_trajs), aes(distance))+ 
	geom_histogram(bins = 50, fill = 'grey80', color = 'black')+
	scale_x_continuous('Travelled distance (mm)', breaks = seq(0, 600, 150))+
	geom_segment(x = 200, xend = 200, yend = Inf, y = 0, linetype = 2)+ylab('Frequency')


################################# STO EXPERIMENTS #################################
library(trajr)

foods <- rbindlist(lapply(sto, function(i){
	d <- rbindlist(i@food)
	idx <- which.min(d[['t']])
	d[idx]
}))[, node := get_node(x, y), by = seq_along(sto)]

ids <- lapply(seq_along(sto), function(i){
	id <- setDT(sto[[i]]@data)[Frame == foods[i, t] & node == foods[i, node]][['N_ind']]
	data <- sto[[i]]@data[N_ind == id]
	Dt <- diff(data[, Frame])
	stops <- which(Dt != 1)
	stop_duration <- Dt[Dt != 1]
	list(data = data, stops = stops, stop_duration = stop_duration, interactions = sum(data[['Crossings']]))
})


exp_trajs <- vapply(seq_along(sto), function(i){
	xy <- TrajFromCoords(ids[[i]][['data']][Frame <= foods[i, t], c('Xmm', 'Ymm', 'Frame')], fps = 2)
	TrajStraightness(xy)
}, numeric(1))

exploration_trajs <- do.call('rbind.data.frame', lapply(seq_along(sto), function(i){
	data <- setDT(sto[[i]]@data)[Frame <= foods[i, t]]
	inds <- unique(data[['N_ind']])
	t(vapply(seq_along(inds), function(ii){
		xy <- TrajFromCoords(data[N_ind == inds[ii], c('Xmm', 'Ymm', 'Frame')], fps = 2)
		c(TrajStraightness(xy), nrow(xy), TrajDistance(xy))
	}, numeric(3)))
}))
colnames(exploration_trajs) <- c('straightness', 'length', 'distance')

ggplot(data = exploration_trajs,aes('Explorers', straightness)) + geom_boxplot() +
	geom_jitter(width = 0.2, size = 3, alpha = 0.2)

ggplot(data = as.data.frame(exploration_trajs), aes(length, straightness))+geom_point()
ggplot(data = as.data.frame(exploration_trajs), aes('', distance))+geom_boxplot(outlier.shape = NA)+
	geom_jitter()
ggplot(data = as.data.frame(exploration_trajs), aes(distance))+ 
	geom_histogram(bins = 50, fill = 'grey80', color = 'black')+
	scale_x_continuous('Travelled distance (mm)', breaks = seq(0, round(max(exploration_trajs[['distance']])*150)/150, 150))+
	geom_segment(x = 200, xend = 200, yend = Inf, y = 0, linetype = 2)+ylab('Frequency')


################################# NF EXPERIMENTS #################################
library(trajr)
# foods <- rbindlist(lapply(det, function(i){
# 	d <- rbindlist(i@food)
# 	idx <- which.min(d[['t']])
# 	d[idx]
# }))[, node := get_node(x, y), by = seq_along(det)]
maxt <- 2400
# ids <- lapply(seq_along(nf), function(i){
# 	id <- setDT(nf[[i]]@data)[Frame == maxt & node == foods[i, node]][['N_ind']]
# 	data <- nf[[i]]@data[N_ind == id]
# 	Dt <- diff(data[, Frame])
# 	nfps <- which(Dt != 1)
# 	nfp_duration <- Dt[Dt != 1]
# 	list(data = data, nfps = nfps, nfp_duration = nfp_duration, interactions = sum(data[['Crossings']]))
# })
# 
# 
# exp_trajs <- vapply(seq_along(nf), function(i){
# 	xy <- TrajFromCoords(ids[[i]][['data']][Frame <= foods[i, t], c('Xmm', 'Ymm', 'Frame')], fps = 2)
# 	TrajStraightness(xy)
# }, numeric(1))

exploration_trajs <- do.call('rbind.data.frame', lapply(seq_along(nf), function(i){
	data <- setDT(nf[[i]]@data)[Frame <= maxt]
	inds <- unique(data[['N_ind']])
	t(vapply(seq_along(inds), function(ii){
		xy <- TrajFromCoords(data[N_ind == inds[ii], c('Xmm', 'Ymm', 'Frame')], fps = 2)
		c(TrajStraightness(xy), nrow(xy), TrajDistance(xy))
	}, numeric(3)))
}))
colnames(exploration_trajs) <- c('straightness', 'length', 'distance')

ggplot(data = exploration_trajs,aes('Explorers', straightness)) + geom_boxplot() +
	geom_jitter(width = 0.2, size = 3, alpha = 0.2)

ggplot(data = as.data.frame(exploration_trajs), aes(length, straightness))+geom_point()
ggplot(data = as.data.frame(exploration_trajs), aes('', distance))+geom_boxplot(outlier.shape = NA)+
	geom_jitter()
ggplot(data = as.data.frame(exploration_trajs), aes(distance))+ 
	geom_histogram(bins = 50, fill = 'grey80', color = 'black')+
	scale_x_continuous('Travelled distance (mm)', breaks = seq(0, 600, 150))+
	geom_segment(x = 200, xend = 200, yend = Inf, y = 0, linetype = 2)+ylab('Frequency')


#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################

##############   SIMULATIONS AND EXPERIMENTS COMPARISON EXPLORATION ##########

path <- '~/research/gits/AutomatAnts/results/2024/nest_avoid_explor/'
f <- list.files(path)
files <- f[grepl('data',f)]

maxd <- 10 # 650
maxt <- 2400
draw_hexagons(sim_edges) + 
	geom_circle(center = as.data.frame(hex_sim[node == '(0, 22)', c('x', 'y')]), r = 10, npoints = 100) +
	geom_point(data = hex_sim[node == '(0, 22)', c('x', 'y')], aes(x, y)) + 
	scale_y_continuous(limits = c(0, 20))

revert_node <- function(n){
	n <- do.call('rbind', n)
	apply(n, 1, function(i) paste0('(', paste0(i, collapse = ', '), ')'))
}
parse_nodes <- function(nodes){
	unlist(strsplit(nodes, ';'))
}

parse_ids <- function(ids){
	as.integer(unlist(strsplit(ids, ',')))
}

## calculate velocity, acceleration, straightness, distance travelled, time in arena...
exploration_trajs_sim <- do.call('rbind.data.frame', lapply(seq_along(files), function(i){
	library(arrow)
	do.call('gc', args = list(verbose = FALSE))
	data <- setDT(read_parquet(paste0(path, files[i])))[Frame <= maxt, c('pos', 'id_out', 'Frame')][, .(id = parse_ids(id_out), node = parse_nodes(pos)), by = 'Frame']
	mdata <- merge(hex_sim, data, by = 'node')
	dmatrix <- pdist(as.matrix(mdata[, c('x', 'y')]), as.matrix(hex_sim[node == '(0, 22)', c('x', 'y')]))
	inds <- unique(mdata[as.numeric(dmatrix) > maxd, id])
	# filtr <- unique(data[N_ind %in% inds & Crossings > 0, N_ind])
	set(mdata, j = 'd', value = dmatrix)
	mdata <- mdata[id %in% inds]
	idx <- mdata[, .(idx = which.max((d - maxd)> 0)), by = 'id'][['idx']]
	
	t(vapply(seq_along(inds), function(ii){
		sbst <- mdata[id == inds[ii]]
		xy <- TrajFromCoords(sbst[1:idx[ii], c('x', 'y', 'Frame')], fps = 2)
		if(nrow(xy) > 10){
			c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
			  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
			  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
		} else {rep(0, 6)}

	}, numeric(6)))
}))
colnames(exploration_trajs_sim) <- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')

high_explrs <- rbindlist(list(det = exploration_trajs_det, 
			      nf = exploration_trajs_nf, 
			      sto = exploration_trajs_sto, 
			      sim = exploration_trajs_sim), idcol = TRUE)

ggplot(data = high_explrs,
       aes(.id, Straightness, fill = .id)) + 
	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9),
		     show.legend = FALSE)+
	geom_jitter(aes(color = .id), 
		    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
		    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition', 
			  values = c('mediumpurple', 'gold3', 'brown4', 'grey'))+
	scale_color_manual(values = rep('black', 4))+
	theme(legend.title = element_blank())+
	scale_x_discrete('', labels = c('DET', 'NFD', 'Simulations', 'STO')) +
	ylab('Net-to-Gross ratio') + theme(legend.title = element_blank())

ggplot(data = high_explrs,
       aes(.id, Diffusion, fill = .id)) + 
	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9),
		     show.legend = FALSE)+
	geom_jitter(aes(color = .id), 
		    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
		    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition', 
			  values = c('mediumpurple', 'gold3', 'brown4', 'grey'))+
	scale_color_manual(values = rep('black', 4))+
	theme(legend.title = element_blank())+
	scale_x_discrete('', labels = c('DET', 'NFD', 'Simulations', 'STO')) +
	ylab('MSD ??') + theme(legend.title = element_blank())



#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################

#################################################################
################### EXPLORATION NF ##############################
#################################################################
draw_hexagons(add =geom_foodpatches()) + 
	geom_circle(hex[hex$node == 634, c('x', 'y')], r = 500, npoints = 200) +
	geom_point(data = hex[hex$node == 634, c('x', 'y')], aes(x, y)) + 
	scale_y_continuous(limits = c(1000, 2000))
maxd <- 500 # 650
maxt <- 2400

## calculate velocity, acceleration, straightness, distance travelled, time in arena...
exploration_trajs_nf <- do.call('rbind.data.frame', lapply(seq_along(nf), function(i){
	data <- setDT(nf[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	# filtr <- unique(data[N_ind %in% inds & Crossings > 0, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	idx <- data[, .(idx = which.max((d - maxd)> 0)), by = 'N_ind'][['idx']]
	# inds <- unique(data[['N_ind']])
	# t(vapply(seq_along(inds), function(ii){
	# 	xy <- TrajFromCoords(data[N_ind == inds[ii], c('Xmm', 'Ymm', 'Frame')], fps = 2)
	# 	c(TrajStraightness(xy), nrow(xy), TrajDistance(xy))
	# }, numeric(3)))
	
	t(vapply(seq_along(inds), function(ii){
		sbst <- data[N_ind == inds[ii]]
		xy <- TrajFromCoords(sbst[1:idx[ii], c('Xmm', 'Ymm', 'Frame')], fps = 2)
		c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
		  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
		  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
	}, numeric(6)))
}))
colnames(exploration_trajs_nf) <- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')
low_expl_trajs_nf <- do.call('rbind.data.frame', lapply(seq_along(nf), function(i){
	data <- setDT(nf[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	inds <- unique(data[!N_ind %in% inds, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	# idx <- data[, .(idx = which.max((d - maxd)> 0)), by = 'N_ind'][['idx']]
	t(vapply(seq_along(inds), function(ii){
		sbst <- data[N_ind == inds[ii]]
		if(nrow(sbst) > 10){
			xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
			c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
			  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
			  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
		} else {rep(0, 6)}
		
	}, numeric(6)))
}))
colnames(low_expl_trajs_nf) <- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')

######## !!!!! WITH MODIFICATIONS !!!!! ########
low_expl_trajs_nf <- do.call('rbind.data.frame', lapply(seq_along(nf), function(i){
	data <- setDT(nf[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	inds <- unique(data[!N_ind %in% inds, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	idx <- data[, .(idx = which.min(abs((d - maxd)))), by = 'N_ind'][['idx']]
	t(vapply(seq_along(inds), function(ii){
		sbst <- data[N_ind == inds[ii]]
		xy <- TrajFromCoords(sbst[1:idx[[ii]], c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 10){
			c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
			  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
			  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
		} else {rep(0, 6)}
		
	}, numeric(6)))
}))
colnames(low_expl_trajs_nf) <- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')

#################################################################
#################################################################

#################################################################
#################################################################


#################################################################
################### EXPLORATION DET #############################
#################################################################
# draw_hexagons(add =geom_foodpatches()) +
# 	geom_circle(hex[hex$node == 634, c('x', 'y')], r = 500, npoints = 200) +
# 	geom_point(data = hex[hex$node == 634, c('x', 'y')], aes(x, y)) +
# 	scale_y_continuous(limits = c(1000, 2000))
maxd <- 500 # 650
maxt <- 2400

## calculate velocity, acceleration, straightness, distance travelled, time in arena...
exploration_trajs_det <- do.call('rbind.data.frame', lapply(seq_along(det), function(i){
	data <- setDT(det[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	filtr <- unique(data[N_ind %in% inds & Crossings > 0, N_ind])
	inds <- inds[!inds %in% filtr]
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	idx <- data[, .(idx = which.max((d - maxd)> 0)), by = 'N_ind'][['idx']]
	# inds <- unique(data[['N_ind']])
	# t(vapply(seq_along(inds), function(ii){
	# 	xy <- TrajFromCoords(data[N_ind == inds[ii], c('Xmm', 'Ymm', 'Frame')], fps = 2)
	# 	c(TrajStraightness(xy), nrow(xy), TrajDistance(xy))
	# }, numeric(3)))
	
	t(vapply(seq_along(inds), function(ii){
		sbst <- data[N_ind == inds[ii]]
		xy <- TrajFromCoords(sbst[1:idx[ii], c('Xmm', 'Ymm', 'Frame')], fps = 2)
		c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
		  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
		  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
	}, numeric(6)))
}))
colnames(exploration_trajs_det) <- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')


######## WITHOUT MODIFICATIONS ########
low_expl_trajs_det <- do.call('rbind.data.frame', lapply(seq_along(det), function(i){
	data <- setDT(det[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	inds <- unique(data[!N_ind %in% inds, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	# idx <- data[, .(idx = which.max((d - maxd)> 0)), by = 'N_ind'][['idx']]
	t(vapply(seq_along(inds), function(ii){
		sbst <- data[N_ind == inds[ii]]
		if(nrow(sbst) > 20){
			xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
			c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
			  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
			  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
		} else {rep(0, 6)}

	}, numeric(6)))
}))
colnames(low_expl_trajs_det) <- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')



######## !!!!! WITH MODIFICATIONS !!!!! ########
low_expl_trajs_det <- do.call('rbind.data.frame', lapply(seq_along(det), function(i){
	data <- setDT(det[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	inds <- unique(data[!N_ind %in% inds, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	idx <- data[, .(idx = which.min(abs((d - maxd)))), by = 'N_ind'][['idx']]
	t(vapply(seq_along(inds), function(ii){
		sbst <- data[N_ind == inds[ii]]
		xy <- TrajFromCoords(sbst[1:idx[[ii]], c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 20){
			c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
			  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
			  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
		} else {rep(0, 6)}
		
	}, numeric(6)))
}))
colnames(low_expl_trajs_det) <- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')

#################################################################
################### EXPLORATION STO #############################
#################################################################

maxd <- 500 # 650
maxt <- 2400

## calculate velocity, acceleration, straightness, distance travelled, time in arena...
exploration_trajs_sto <- do.call('rbind.data.frame', lapply(seq_along(sto), function(i){
	data <- setDT(sto[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	filtr <- unique(data[N_ind %in% inds & Crossings > 0, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	idx <- data[, .(idx = which.max((d - maxd)> 0)), by = 'N_ind'][['idx']]
	# inds <- unique(data[['N_ind']])
	# t(vapply(seq_along(inds), function(ii){
	# 	xy <- TrajFromCoords(data[N_ind == inds[ii], c('Xmm', 'Ymm', 'Frame')], fps = 2)
	# 	c(TrajStraightness(xy), nrow(xy), TrajDistance(xy))
	# }, numeric(3)))
	
	t(vapply(seq_along(inds), function(ii){
		sbst <- data[N_ind == inds[ii]]
		xy <- TrajFromCoords(sbst[1:idx[ii], c('Xmm', 'Ymm', 'Frame')], fps = 2)
		c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
		  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
		  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
	}, numeric(6)))
}))
colnames(exploration_trajs_sto) <- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')
low_expl_trajs_sto <- do.call('rbind.data.frame', lapply(seq_along(sto), function(i){
	data <- setDT(sto[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	inds <- unique(data[!N_ind %in% inds, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	# idx <- data[, .(idx = which.max((d - maxd)> 0)), by = 'N_ind'][['idx']]
	t(vapply(seq_along(inds), function(ii){
		sbst <- data[N_ind == inds[ii]]
		if(nrow(sbst) > 20){
			xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
			c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
			  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
			  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
		} else {rep(0, 6)}
		
	}, numeric(6)))
}))
colnames(low_expl_trajs_sto) <- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')

######## !!!!! WITH MODIFICATIONS !!!!! ########
low_expl_trajs_sto <- do.call('rbind.data.frame', lapply(seq_along(sto), function(i){
	data <- setDT(sto[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	inds <- unique(data[!N_ind %in% inds, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	idx <- data[, .(idx = which.min(abs((d - maxd)))), by = 'N_ind'][['idx']]
	t(vapply(seq_along(inds), function(ii){
		sbst <- data[N_ind == inds[ii]]
		xy <- TrajFromCoords(sbst[1:idx[[ii]], c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 10){
			c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
			  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
			  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
		} else {rep(0, 6)}
		
	}, numeric(6)))
}))
colnames(low_expl_trajs_sto) <- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')

#################################################################
#################################################################

#################################################################
#################################################################
exploration_trajs <- rbindlist(list(det = exploration_trajs_det, 
				    sto = exploration_trajs_sto,
				    nf = exploration_trajs_nf), idcol = TRUE)
set(exploration_trajs, j = 'expl', value = 'Long-range scouts')

# ggplot(data = exploration_trajs, aes(.id, Straightness, fill = .id)) + geom_boxplot()

low_exploration_trajs <- rbindlist(list(det = low_expl_trajs_det, 
				    sto = low_expl_trajs_sto,
				    nf = low_expl_trajs_nf), idcol = TRUE)

set(low_exploration_trajs, j = 'expl', value = 'Short-range scouts')

ggplot(data = low_exploration_trajs, aes(.id, Straightness, fill = .id)) + geom_boxplot()

ggplot(data = rbindlist(list(exploration_trajs, 
			     low_exploration_trajs[Mean_v > 0])),
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
	xlab('') + ylab('Net-to-Gross ratio')


# ggplot(data = rbindlist(list(exploration_trajs[, expl := 'Long-range scouts'], 
# 			     low_exploration_trajs[Mean_v > 0, expl := 'Short-range scouts'])),
#        aes(expl, Straightness, fill = .id)) + 
# 	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
# 	geom_jitter(aes(color = .id), 
# 		    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
# 		    size = 3, show.legend = FALSE, alpha = 0.25)+
# 	scale_fill_manual('Experimental condition', 
# 			  values = c('mediumpurple', 'gold3', 'brown4'),
# 			  labels = c('DET', 'NFD', 'STO'))+
# 	scale_color_manual(values = rep('black', 3))+
# 	theme(legend.title = element_blank())+
# 	xlab('') + ylab('Net-to-Gross ratio')
	
	
ggplot(data = rbindlist(list(exploration_trajs, 
			     low_exploration_trajs[Mean_v > 0])),
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
	xlab('') + scale_y_continuous('Distance travelled (mm)', 
				      breaks = seq(0, 6000, 1500))

ggplot(data = rbindlist(list(exploration_trajs, 
			     low_exploration_trajs[Mean_v > 0])),
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
	xlab('') + scale_y_continuous('MSD ??')
	

ggplot(data = rbindlist(list(exploration_trajs, 
			     low_exploration_trajs[Mean_v > 0])),
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
	xlab('') + scale_y_continuous('< Velocity > (mm/s)')

ggplot(data = rbindlist(list(exploration_trajs, 
			     low_exploration_trajs[Mean_v > 0])),
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
	xlab('') + scale_y_continuous(latex2exp::TeX('< Acceleration > ($mm/s^{2}$)'))

ggplot(data = rbindlist(list(exploration_trajs, 
			     low_exploration_trajs[Mean_v > 0])),
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
	xlab('') + scale_y_continuous('Exploration time (min)', breaks = seq(0, 10, 1.5))


#################################################################
#################################################################

#################################################################
##################### NTG - DISTANCE FOR NF EXPS ################
#################################################################
maxd <- 500 # 650
maxt <- 2400

nf_NTG_d <- lapply(seq(50, 1000, 50), function(d){
	maxd <- d
	## calculate velocity, acceleration, straightness, distance travelled, time in arena...
	do.call('rbind.data.frame', lapply(seq_along(nf), function(i){
		data <- setDT(nf[[i]]@data)[Frame <= maxt]
		dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
		inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
		set(data, j = 'd', value = dmatrix)
		data <- data[N_ind %in% inds]
		idx <- data[, .(idx = which.max((d - maxd)> 0)), by = 'N_ind'][['idx']]
		t(vapply(seq_along(inds), function(ii){
			sbst <- data[N_ind == inds[ii]]
			data_sbst <- sbst[1:idx[ii], c('Xmm', 'Ymm', 'Frame')]
			if(nrow(data_sbst) > 10){
				xy <- TrajFromCoords(data_sbst, fps = 2)
				c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
				  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
				  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
			} else {rep(0, 6)}

		}, numeric(6)))
	}))
})

nf_dt <- rbindlist(nf_NTG_d, idcol = TRUE)
colnames(nf_dt) <- c('.id', 'Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')

ggplot(data = nf_dt,
       aes(.id * 50, Straightness, group=.id)) + 
	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
	# geom_jitter(
	# 	    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
	# 	    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition', 
			  values = c('mediumpurple', 'gold3', 'brown4'),
			  labels = c('DET', 'NFD', 'STO'))+
	scale_color_manual(values = rep('black', 3))+
	theme(legend.title = element_blank())+
	xlab('') + ylab('Net-to-Gross ratio')

ggplot(data = nf_dt[, .(Straightness = median(Straightness)), by = '.id'],
       aes(.id * 50, Straightness, group=.id)) + 
	geom_point()+
	# geom_jitter(
	# 	    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
	# 	    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition', 
			  values = c('mediumpurple', 'gold3', 'brown4'),
			  labels = c('DET', 'NFD', 'STO'))+
	scale_color_manual(values = rep('black', 3))+
	theme(legend.title = element_blank())+
	xlab('') + ylab('Net-to-Gross ratio')

nf_NTG_d_low <- lapply(seq(50, 1000, 50), function(d){
	maxd <- d
	## calculate velocity, acceleration, straightness, distance travelled, time in arena...
	do.call('rbind.data.frame', lapply(seq_along(nf), function(i){
		data <- setDT(nf[[i]]@data)[Frame <= maxt]
		dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
		inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
		inds <- unique(data[!N_ind %in% inds, N_ind])
		set(data, j = 'd', value = dmatrix)
		data <- data[N_ind %in% inds]
		idx <- data[, .(idx = which.min(abs((d - maxd)))), by = 'N_ind'][['idx']]
		t(vapply(seq_along(inds), function(ii){
			sbst <- data[N_ind == inds[ii]]
			xy <- TrajFromCoords(sbst[1:idx[[ii]], c('Xmm', 'Ymm', 'Frame')], fps = 2)
			if(nrow(xy) > 10){
				c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
				  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
				  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
			} else {rep(0, 6)}
			
		}, numeric(6)))
	}))
})

nf_dt_low <- rbindlist(nf_NTG_d_low, idcol = TRUE)
colnames(nf_dt_low) <- c('.id', 'Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')

ggplot(data = nf_dt_low,
       aes(.id * 50, Straightness, group=.id)) + 
	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
	# geom_jitter(
	# 	    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
	# 	    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition', 
			  values = c('mediumpurple', 'gold3', 'brown4'),
			  labels = c('DET', 'NFD', 'STO'))+
	scale_color_manual(values = rep('black', 3))+
	theme(legend.title = element_blank())+
	xlab('') + ylab('Net-to-Gross ratio')

ggplot(data = nf_dt_low[, .(Straightness = median(Straightness)), by = '.id'],
       aes(.id * 50, Straightness, group=.id)) + 
	geom_point()+
	# geom_jitter(
	# 	    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
	# 	    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition', 
			  values = c('mediumpurple', 'gold3', 'brown4'),
			  labels = c('DET', 'NFD', 'STO'))+
	scale_color_manual(values = rep('black', 3))+
	theme(legend.title = element_blank())+
	xlab('') + ylab('Net-to-Gross ratio')

## calculate velocity, acceleration, straightness, distance travelled, time in arena...
exploration_trajs_nf <- do.call('rbind.data.frame', lapply(seq_along(nf), function(i){
	data <- setDT(nf[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	inds <- unique(data[!N_ind %in% inds, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	idx <- data[, .(idx = which.min(abs((d - maxd)))), by = 'N_ind'][['idx']]
	t(vapply(seq_along(inds), function(ii){
		sbst <- data[N_ind == inds[ii]]
		xy <- TrajFromCoords(sbst[1:idx[[ii]], c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 10){
			c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
			  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
			  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
		} else {rep(0, 6)}
		
	}, numeric(6)))
}))
#################################################################
#################################################################

#################################################################
#################### NTG - DISTANCE FOR DET EXPS ################
#################################################################
maxd <- 500 # 650
maxt <- 2400

det_NTG_d <- lapply(seq(50, 1000, 50), function(d){
	maxd <- d
	## calculate velocity, acceleration, straightness, distance travelled, time in arena...
	do.call('rbind.data.frame', lapply(seq_along(det), function(i){
		data <- setDT(det[[i]]@data)[Frame <= maxt]
		dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
		inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
		set(data, j = 'd', value = dmatrix)
		data <- data[N_ind %in% inds]
		idx <- data[, .(idx = which.max((d - maxd)> 0)), by = 'N_ind'][['idx']]
		t(vapply(seq_along(inds), function(ii){
			sbst <- data[N_ind == inds[ii]]
			data_sbst <- sbst[1:idx[ii], c('Xmm', 'Ymm', 'Frame')]
			if(nrow(data_sbst) > 10){
				xy <- TrajFromCoords(data_sbst, fps = 2)
				c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
				  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
				  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
			} else {rep(0, 6)}
			
		}, numeric(6)))
	}))
})

det_dt <- rbindlist(det_NTG_d, idcol = TRUE)
colnames(det_dt) <- c('.id', 'Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')

ggplot(data = det_dt,
       aes(.id * 50, Straightness, group=.id)) + 
	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
	# geom_jitter(
	# 	    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
	# 	    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition', 
			  values = c('mediumpurple', 'gold3', 'brown4'),
			  labels = c('DET', 'detD', 'STO'))+
	scale_color_manual(values = rep('black', 3))+
	theme(legend.title = element_blank())+
	xlab('') + ylab('Net-to-Gross ratio')

ggplot(data = det_dt[, .(Straightness = median(Straightness)), by = '.id'],
       aes(.id * 50, Straightness, group=.id)) + 
	geom_point()+
	# geom_jitter(
	# 	    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
	# 	    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition', 
			  values = c('mediumpurple', 'gold3', 'brown4'),
			  labels = c('DET', 'detD', 'STO'))+
	scale_color_manual(values = rep('black', 3))+
	theme(legend.title = element_blank())+
	xlab('') + ylab('Net-to-Gross ratio')

det_NTG_d_low <- lapply(seq(50, 1000, 50), function(d){
	maxd <- d
	## calculate velocity, acceleration, straightness, distance travelled, time in arena...
	do.call('rbind.data.frame', lapply(seq_along(det), function(i){
		data <- setDT(det[[i]]@data)[Frame <= maxt]
		dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
		inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
		inds <- unique(data[!N_ind %in% inds, N_ind])
		set(data, j = 'd', value = dmatrix)
		data <- data[N_ind %in% inds]
		idx <- data[, .(idx = which.min(abs((d - maxd)))), by = 'N_ind'][['idx']]
		t(vapply(seq_along(inds), function(ii){
			sbst <- data[N_ind == inds[ii]]
			xy <- TrajFromCoords(sbst[1:idx[[ii]], c('Xmm', 'Ymm', 'Frame')], fps = 2)
			if(nrow(xy) > 10){
				c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
				  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
				  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
			} else {rep(0, 6)}
			
		}, numeric(6)))
	}))
})

det_dt_low <- rbindlist(det_NTG_d_low, idcol = TRUE)
colnames(det_dt_low) <- c('.id', 'Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')

ggplot(data = det_dt_low,
       aes(.id * 50, Straightness, group=.id)) + 
	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
	# geom_jitter(
	# 	    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
	# 	    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition', 
			  values = c('mediumpurple', 'gold3', 'brown4'),
			  labels = c('DET', 'detD', 'STO'))+
	scale_color_manual(values = rep('black', 3))+
	theme(legend.title = element_blank())+
	xlab('') + ylab('Net-to-Gross ratio')

ggplot(data = det_dt_low[, .(Straightness = median(Straightness)), by = '.id'],
       aes(.id * 50, Straightness, group=.id)) + 
	geom_point()+
	# geom_jitter(
	# 	    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
	# 	    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition', 
			  values = c('mediumpurple', 'gold3', 'brown4'),
			  labels = c('DET', 'detD', 'STO'))+
	scale_color_manual(values = rep('black', 3))+
	theme(legend.title = element_blank())+
	xlab('') + ylab('Net-to-Gross ratio')

#################################################################
#################################################################

#################################################################
#################### NTG - DISTANCE FOR STO EXPS ################
#################################################################
maxd <- 500 # 650
maxt <- 2400

sto_NTG_d <- lapply(seq(50, 1000, 50), function(d){
	maxd <- d
	## calculate velocity, acceleration, straightness, distance travelled, time in arena...
	do.call('rbind.data.frame', lapply(seq_along(sto), function(i){
		data <- setDT(sto[[i]]@data)[Frame <= maxt]
		dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
		inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
		set(data, j = 'd', value = dmatrix)
		data <- data[N_ind %in% inds]
		idx <- data[, .(idx = which.max((d - maxd)> 0)), by = 'N_ind'][['idx']]
		t(vapply(seq_along(inds), function(ii){
			sbst <- data[N_ind == inds[ii]]
			data_sbst <- sbst[1:idx[ii], c('Xmm', 'Ymm', 'Frame')]
			if(nrow(data_sbst) > 10){
				xy <- TrajFromCoords(data_sbst, fps = 2)
				c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
				  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
				  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
			} else {rep(0, 6)}
			
		}, numeric(6)))
	}))
})

sto_dt <- rbindlist(sto_NTG_d, idcol = TRUE)
colnames(sto_dt) <- c('.id', 'Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')

ggplot(data = sto_dt,
       aes(.id * 50, Straightness, group=.id)) + 
	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
	# geom_jitter(
	# 	    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
	# 	    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition', 
			  values = c('mediumpurple', 'gold3', 'brown4'),
			  labels = c('sto', 'stoD', 'STO'))+
	scale_color_manual(values = rep('black', 3))+
	theme(legend.title = element_blank())+
	xlab('') + ylab('Net-to-Gross ratio')

ggplot(data = sto_dt[, .(Straightness = median(Straightness)), by = '.id'],
       aes(.id * 50, Straightness, group=.id)) + 
	geom_point()+
	# geom_jitter(
	# 	    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
	# 	    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition', 
			  values = c('mediumpurple', 'gold3', 'brown4'),
			  labels = c('sto', 'stoD', 'STO'))+
	scale_color_manual(values = rep('black', 3))+
	theme(legend.title = element_blank())+
	xlab('') + ylab('Net-to-Gross ratio')

sto_NTG_d_low <- lapply(seq(50, 1000, 50), function(d){
	maxd <- d
	## calculate velocity, acceleration, straightness, distance travelled, time in arena...
	do.call('rbind.data.frame', lapply(seq_along(sto), function(i){
		data <- setDT(sto[[i]]@data)[Frame <= maxt]
		dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
		inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
		inds <- unique(data[!N_ind %in% inds, N_ind])
		set(data, j = 'd', value = dmatrix)
		data <- data[N_ind %in% inds]
		idx <- data[, .(idx = which.min(abs((d - maxd)))), by = 'N_ind'][['idx']]
		t(vapply(seq_along(inds), function(ii){
			sbst <- data[N_ind == inds[ii]]
			xy <- TrajFromCoords(sbst[1:idx[[ii]], c('Xmm', 'Ymm', 'Frame')], fps = 2)
			if(nrow(xy) > 10){
				c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
				  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
				  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
			} else {rep(0, 6)}
			
		}, numeric(6)))
	}))
})

sto_dt_low <- rbindlist(sto_NTG_d_low, idcol = TRUE)
colnames(sto_dt_low) <- c('.id', 'Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')

ggplot(data = sto_dt_low,
       aes(.id * 50, Straightness, group=.id)) + 
	geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.9))+
	# geom_jitter(
	# 	    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
	# 	    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition', 
			  values = c('mediumpurple', 'gold3', 'brown4'),
			  labels = c('sto', 'stoD', 'STO'))+
	scale_color_manual(values = rep('black', 3))+
	theme(legend.title = element_blank())+
	xlab('') + ylab('Net-to-Gross ratio')

ggplot(data = sto_dt_low[, .(Straightness = median(Straightness)), by = '.id'],
       aes(.id * 50, Straightness, group=.id)) + 
	geom_point()+
	# geom_jitter(
	# 	    position = position_jitterdodge(dodge.width = 0.9,jitter.width = 0.15),
	# 	    size = 3, show.legend = FALSE, alpha = 0.25)+
	scale_fill_manual('Experimental condition', 
			  values = c('mediumpurple', 'gold3', 'brown4'),
			  labels = c('sto', 'stoD', 'STO'))+
	scale_color_manual(values = rep('black', 3))+
	theme(legend.title = element_blank())+
	xlab('') + ylab('Net-to-Gross ratio')

maxd <- 500 # 650
maxt <- 2400

## calculate velocity, acceleration, straightness, distance travelled, time in arena...
exploration_trajs_nf <- do.call('rbind.data.frame', lapply(seq_along(nf), function(i){
	data <- setDT(nf[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	idx <- data[, .(idx = which.max((d - maxd)> 0)), by = 'N_ind'][['idx']]
	t(vapply(seq_along(inds), function(ii){
		sbst <- data[N_ind == inds[ii]]
		xy <- TrajFromCoords(sbst[1:idx[ii], c('Xmm', 'Ymm', 'Frame')], fps = 2)
		c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
		  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
		  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE), nrow(xy) / 120)
	}, numeric(6)))
}))


#################################################################
#################################################################

#################################################################
#################################################################

tmax <- 20*120
min_pos <- 15
## ALL NF EXPERIMENTS
M_nf <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
for(e in seq_along(nf)){
	exp <- nf[[e]]
	dt <- setDT(exp@data)[Frame <= (min(Frame) + tmax), c('N_ind', 'node')]
	Ns <- dt[, .N, N_ind]
	inds <- Ns[N >= min_pos, N_ind]
	result_list <- vector('list', length(inds))
	M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
	for(j in seq_along(inds)){
		df_ind <- exp@data[N_ind == inds[j], node]
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
					M_ind[as.character(result[ii-1]), as.character(result[ii])] <-
						M_ind[as.character(result[ii-1]), as.character(result[ii])]+1
					M_nf[as.character(result[ii-1]), as.character(result[ii])] <-
						M_nf[as.character(result[ii-1]), as.character(result[ii])]+1
				}
				
				
				M_ind / sum(M_ind)
				result_list[[j]] <- apply(M_ind / sum(M_ind), 2, function(i) i / sum(i))
			}
		}
	}
	M
	apply(M / sum(M), 2, function(i) i / sum(i))
}
apply(M_nf / sum(M_nf), 2, function(i) i / sum(i))


#################################################################
## ALL DETERMINIST EXPS


min_pos <- 15

M_det <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
for(e in seq_along(det)){
	exp <- det[[e]]
	tmax <- min(rbindlist(exp@food)[['t']])
	dt <- setDT(exp@data)[Frame <= tmax, c('N_ind', 'node')]
	Ns <- dt[, .N, N_ind]
	inds <- Ns[N >= min_pos, N_ind]
	result_list <- vector('list', length(inds))
	M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
	for(j in seq_along(inds)){
		df_ind <- exp@data[N_ind == inds[j], node]
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
					M_ind[as.character(result[ii-1]), as.character(result[ii])] <-
						M_ind[as.character(result[ii-1]), as.character(result[ii])]+1
					M_det[as.character(result[ii-1]), as.character(result[ii])] <-
						M_det[as.character(result[ii-1]), as.character(result[ii])]+1
				}
				
				
				M_ind / sum(M_ind)
				result_list[[j]] <- apply(M_ind / sum(M_ind), 2, function(i) i / sum(i))
			}
		}
	}
	M
	apply(M / sum(M), 2, function(i) i / sum(i))
}
apply(M_det / sum(M_det), 2, function(i) i / sum(i))


#################################################################
## ALL DETERMINIST EXPS

foods <- rbindlist(lapply(det, function(i){
	d <- rbindlist(i@food)
	idx <- which.min(d[['t']])
	d[idx]
}))[, node := get_node(x, y), by = seq_along(det)]


ids <- lapply(seq_along(det), function(i){
	id <- setDT(det[[i]]@data)[Frame == foods[i, t] & node == foods[i, node]][['N_ind']]
	data <- det[[i]]@data[N_ind == id]
	Dt <- diff(data[, Frame])
	stops <- which(Dt != 1)
	stop_duration <- Dt[Dt != 1]
	list(data = data, stops = stops, stop_duration = stop_duration, interactions = sum(data[['Crossings']]))
})



min_pos <- 15

M_det <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
full_det_result <- list()
for(e in seq_along(det)){
	exp <- det[[e]]
	tmax <- min(rbindlist(exp@food)[['t']])
	dt <- setDT(exp@data)[Frame <= tmax, c('N_ind', 'node')]
	Ns <- dt[, .N, N_ind]
	inds <- Ns[N >= min_pos, N_ind]
	result_list <- vector('list', length(inds))
	M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
	for(j in seq_along(inds)){
		df_ind <- exp@data[N_ind == inds[j], node]
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
					M_ind[as.character(result[ii-1]), as.character(result[ii])] <-
						M_ind[as.character(result[ii-1]), as.character(result[ii])]+1
					M_det[as.character(result[ii-1]), as.character(result[ii])] <-
						M_det[as.character(result[ii-1]), as.character(result[ii])]+1
				}
				
				
				# M_ind / sum(M_ind)
				result_list[[j]] <- apply(M_ind / sum(M_ind), 2, function(i) i / sum(i))
				full_det_result <- append(full_det_result, list(M_ind))
			}
		}
	}
	# M
	# apply(M / sum(M), 2, function(i) i / sum(i))
}

apply(M_det / sum(M_det), 2, function(i) i / sum(i))
na <- apply(full_det_result[[4]] / sum(full_det_result[[4]]), 2, function(i) i / sum(i))
na_data <- melt(na)
ggplot(data = na_data, aes(Var1, Var2)) + 
	geom_raster(aes(fill = value)) +
	geom_text(aes(label = round(value, 2)))+
	scale_y_continuous('Previous movement', breaks = seq(-1, 1))+
	scale_x_continuous('Posterior movement', breaks = seq(-1, 1), labels = c('Left', 'Back', 'Right'))

#################################################################
## ALL STOCHASTIC EXPS
min_pos <- 15

M_sto <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
for(e in seq_along(sto)){
	exp <- sto[[e]]
	tmax <- min(rbindlist(exp@food)[['t']])
	dt <- setDT(exp@data)[Frame <= tmax, c('N_ind', 'node')]
	Ns <- dt[, .N, N_ind]
	inds <- Ns[N >= min_pos, N_ind]
	result_list <- vector('list', length(inds))
	M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
	for(j in seq_along(inds)){
		df_ind <- exp@data[N_ind == inds[j], node]
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
					M_ind[as.character(result[ii-1]), as.character(result[ii])] <-
						M_ind[as.character(result[ii-1]), as.character(result[ii])]+1
					M_sto[as.character(result[ii-1]), as.character(result[ii])] <-
						M_sto[as.character(result[ii-1]), as.character(result[ii])]+1
				}
				
				
				M_ind / sum(M_ind)
				result_list[[j]] <- apply(M_ind / sum(M_ind), 2, function(i) i / sum(i))
			}
		}
	}
	M
	apply(M / sum(M), 2, function(i) i / sum(i))
}
apply(M_sto / sum(M_sto), 2, function(i) i / sum(i))



##################################################################

source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
library(arrow)
library(latex2exp)

load('~/research/gits/AnTracks/data/det.RData')

foods <- rbindlist(lapply(det, function(i){
	d <- rbindlist(i@food)
	idx <- which.min(d[['t']])
	d[idx]
}))[, node := get_node(x, y), by = seq_along(det)]


ids <- lapply(seq_along(det), function(i){
	id <- setDT(det[[i]]@data)[Frame == foods[i, t] & node == foods[i, node]][['N_ind']]
	data <- det[[i]]@data[N_ind == id]
	Dt <- diff(data[, Frame])
	stops <- which(Dt != 1)
	stop_duration <- Dt[Dt != 1]
	list(data = data, stops = stops, stop_duration = stop_duration, interactions = sum(data[['Crossings']]))
})

plots <- lapply(seq_along(det), function(i){
	draw_hexagons(add = geom_foodpatches()) +
		geom_point(data = ids[[i]][['data']],
			   aes(Xmm, Ymm, size = 1/Frame, color = Frame), show.legend = FALSE)+
		scale_color_viridis_c()+
		scale_size_continuous(range = c(3, 9))+
		geom_point(data = ids[[i]][['data']][ids[[i]][['stops']]+1][, d := ids[[i]][['stop_duration']]],
			   aes(Xmm, Ymm, size = d), shape = 13, 
			   show.legend = FALSE)+
		ggtitle(paste0('N of interactions: ', ids[[i]][['interactions']]))
	
})	

ggarrange(plotlist = plots)







# nm <- as.numeric(hex[hex$node == foods[1, node], c('x', 'y')] - hex[hex$node == 662, c('x', 'y')])
# nt <- as.numeric(tail(ids[[1]][['data']], 1)[, c('Xmm', 'Ymm')]) - as.numeric(hex[hex$node == 662, c('x', 'y')])

# v_norm <- function(x){
# 	sqrt(x[1]**2 + x[2] **2)
# }
# 
# cross_vec <- function(x, y){
# 	x[1]*y[2] - y[1]*x[2]
# }
# 
# v_norm(cross(nt, nm))



# dist2d <- function(x, origin = as.numeric(hex[hex$node == 662, c('x', 'y')]), destiny) {
# 	if(is.data.frame(x)){
# 		x <- as.numeric(x)
# 	}
# 	if(is.data.frame(origin)){
# 		origin <- as.numeric(origin)
# 	}
# 	if(is.data.frame(destiny)){
# 		destiny <- as.numeric(destiny)
# 	}
# 	v1 <- destiny - origin
# 	v2 <- x - origin
# 	m <- cbind(v1,v2)
# 	abs(det(m))/sqrt(sum(v1*v1))
# }
# 
# ds <- vapply(1:nrow(ids[[1]][['data']]), function(i){
# 	dist2d(x = ids[[1]][['data']][i, c('Xmm', 'Ymm')], destiny = foods[1, c('x', 'y')])
# }, numeric(1))
# 
# boxplot(ds)
# 
# dist2d(x = as.numeric(tail(ids[[1]][['data']], 1)[, c('Xmm', 'Ymm')]),
#        destiny = as.numeric(hex[hex$node == 910, c('x', 'y')]))


library(trajr)
# xy <- TrajFromCoords(ids[[1]][['data']][, c('Xmm', 'Ymm', 'Frame')], fps = 2)
# xy <- TrajFromCoords(ids[[2]][['data']][, c('Xmm', 'Ymm', 'Frame')], fps = 2)
# xy <- TrajFromCoords(ids[[2]][['data']][Frame <= foods[2, t], c('Xmm', 'Ymm', 'Frame')], fps = 2)
# TrajStraightness(xy)
# 
exp_trajs <- vapply(seq_along(det), function(i){
	xy <- TrajFromCoords(ids[[i]][['data']][Frame <= foods[i, t], c('Xmm', 'Ymm', 'Frame')], fps = 2)
	TrajStraightness(xy)
}, numeric(1))




path <- '/home/polfer/research/gits/AutomatAnts/results/2024/nest_avoid/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('ballistic_\\d{1,2}', files))))

revert_node <- function(n){
	n <- do.call('rbind', n)
	apply(n, 1, function(i) paste0('(', paste0(i, collapse = ', '), ')'))
}
parse_nodes <- function(nodes){
	unlist(strsplit(nodes, ';'))
}

parse_ids <- function(ids){
	as.integer(unlist(strsplit(ids, ',')))
}

# 
# food_example <- data.table(arrow::read_parquet(paste0(path, files[2], '_food.parquet')))[, .(node = revert_node(node)), by = 't']
# mint_example <- food_example[which.min(t), .(t, node = node)]
# 
# data_example <- data.table(arrow::read_parquet(paste0(path, files[2],
# 						      '_data.parquet')))[`T` <= mint_example[['t']],
# 						      		   c('T', 'id_out', 'pos')]
# data_example <- data_example[, .(node = parse_nodes(pos), id = parse_ids(id_out)), by = 'T'][, Frame := round(`T`*2)]
# xy <- TrajFromCoords(merge(hex_sim, data_example[id == data_example[node == mint_example[['node']],
# 								    id]] & Frame == mint_example[['Frame']], by = 'node')[, c('x', 'y', 'Frame')], fps = 2)
# TrajStraightness(xy)


sim_trajs <- vapply(seq_along(files), function(i){
	food <- data.table(arrow::read_parquet(paste0(path, files[i], '_food.parquet')))[, .(node = revert_node(node)), by = 't']
	mint <- food[which.min(t), .(Frame = round(t*2), node = node)]
	
	data <- data.table(arrow::read_parquet(paste0(path, files[i],
						      '_data.parquet')))[Frame <= mint[['Frame']],
						      		   c('Frame', 'id_out', 'pos')]
	
	data <- data[, .(node = parse_nodes(pos), id = parse_ids(id_out)), by = 'Frame']
	data <- data[id == unique(data[node == mint[['node']], id])]
	ref <- data.table(Frame = (min(data[['Frame']]):(max(data[['Frame']]))))
	# dt <- merge(hex_sim, data, by = 'node')[, c('x', 'y', 'Frame', 'node')]
	dt <- merge(hex_sim, data, by = 'node')[, c('x', 'y', 'Frame')]
	final_data <- merge(dt, ref, all = TRUE, by = 'Frame')
	set(final_data, j = 'x', value = fillVec(final_data[['x']]))
	set(final_data, j = 'y', value = fillVec(final_data[['y']]))
	xy <- TrajFromCoords(final_data[, c('x', 'y', 'Frame')], fps = 2)
	TrajStraightness(xy)
}, numeric(1)) 

path <- '/home/polfer/research/gits/AutomatAnts/results/2024/default/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('default_\\d{1,2}', files))))

model_trajs <- vapply(seq_along(files), function(i){
	food <- data.table(arrow::read_parquet(paste0(path, files[i], '_food.parquet')))[, .(node = revert_node(node)), by = 't']
	if(any(is.finite(food[['t']]))){
		mint <- food[which.min(t), .(Frame = round(t*2), node = node)]
		
		data <- data.table(arrow::read_parquet(paste0(path, files[i],
							      '_data.parquet')))[Frame <= mint[['Frame']],
							      		   c('Frame', 'id_out', 'pos')]
		
		data <- data[, .(node = parse_nodes(pos), id = parse_ids(id_out)), by = 'Frame']
		data <- data[id == unique(data[node == mint[['node']], id])]
		ref <- data.table(Frame = (min(data[['Frame']]):(max(data[['Frame']]))))
		dt <- merge(hex_sim, data, by = 'node')[, c('x', 'y', 'Frame')]
		final_data <- merge(dt, ref, all = TRUE, by = 'Frame')
		set(final_data, j = 'x', value = fillVec(final_data[['x']]))
		set(final_data, j = 'y', value = fillVec(final_data[['y']]))
		xy <- TrajFromCoords(final_data[, c('x', 'y', 'Frame')], fps = 2)
		TrajStraightness(xy)
	} else {0}
}, numeric(1)) 
# 
# model_trajs <- vapply(seq_along(files), function(i){
# 	
# 	food <- data.table(arrow::read_parquet(paste0(path, files[i], '_food.parquet')))[, .(node = revert_node(node)), by = 't']
# 	if(any(is.finite(food[['t']]))){
# 		mint <- food[which.min(t), .(Frame = round(t*2), node = node)]
# 		
# 		data <- data.table(arrow::read_parquet(paste0(path, files[i],
# 							      '_data.parquet')))[Frame <= mint[['Frame']],
# 							      		   c('Frame', 'id_out', 'pos')]
# 		
# 		data <- data[, .(node = parse_nodes(pos), id = parse_ids(id_out)), by = 'Frame']
# 		
# 		xy <- TrajFromCoords(merge(hex_sim, data[id == unique(data[node == mint[['node']], id])] , by = 'node')[, c('x', 'y', 'Frame')], fps = 2)
# 		TrajStraightness(xy)
# 	} else {0}
# 	
# }, numeric(1)) 
model_trajs[model_trajs == 0] <- NA

ggplot(data = data.frame(type = c(rep('Experiments', length(exp_trajs)),
				  rep('Simulations', length(sim_trajs))), y = c(exp_trajs, sim_trajs)),
       aes(type, y))+
	
	geom_boxplot(aes(fill = type), alpha = 0.6, show.legend = FALSE) +
	geom_jitter(width = 0.15, alpha = 0.25, size = 3)+
	scale_fill_manual('', values = c('mediumpurple', 'gold3'))+
	scale_y_continuous('Net-to-gross ratio', limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2))+
	theme(axis.title.x = element_blank())

ggplot(data = data.frame(type = c(rep('Experiments', length(exp_trajs)),
				  rep('Ballistic movement', length(sim_trajs)),
				  rep('Empirical movement', length(model_trajs))), 
			 y = c(exp_trajs, sim_trajs, model_trajs)),
       aes(factor(type, levels = c('Experiments', 'Ballistic movement', 'Empirical movement')), y))+
	
	geom_boxplot(aes(fill = type), alpha = 0.6, show.legend = FALSE, outlier.shape = NA) +
	geom_jitter(width = 0.15, alpha = 0.25, size = 3)+
	scale_fill_manual('', values = c('gold3', 'brown4', 'mediumpurple'))+
	# scale_fill_manual('', values = c('mediumpurple', 'gold3', 'brown4'))+
	scale_y_continuous('Net-to-gross ratio', limits = c(0, 0.8), breaks = seq(0,0.8, 0.2))+
	theme(axis.title.x = element_blank())

global_eff <- function(path){
	files <- list.files(path)[grepl('food', list.files(path))]
	results <- lapply(files, function(i){
		f <- gsub('_food', '', i)
		d <- data.table(read_parquet(paste0(path, f)))
		m1 <- min(d[N >0, Frame])/2
		food <- data.table(read_parquet(paste0(path, i)))
		mint <- min(food[['t']]) - m1
		1/data.table(mint = mint, maxt = max(food[['t']])-min(food[['t']]))
		# 1/data.table(mint = mint, maxt = max(food[['t']]), t0 = 1/m1)
		# 1/data.table(mint = min(food[['t']]), maxt = max(food[['t']]))
	})
	rbindlist(results, idcol = TRUE)
}


path <- '/home/polfer/research/gits/AutomatAnts/results/2024/nest_avoid_explor/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('ballistic_\\d{1,2}', files))))

explor_trajs <- vapply(seq_along(files), function(i){
	food <- data.table(arrow::read_parquet(paste0(path, files[i], '_food.parquet')))[, .(node = revert_node(node)), by = 't']
	if(any(is.finite(food[['t']]))){
		mint <- food[which.min(t), .(Frame = round(t*2), node = node)]
		
		data <- data.table(arrow::read_parquet(paste0(path, files[i],
							      '_data.parquet')))[Frame <= mint[['Frame']],
							      		   c('Frame', 'id_out', 'pos')]
		
		data <- data[, .(node = parse_nodes(pos), id = parse_ids(id_out)), by = 'Frame']
		data <- data[id == unique(data[node == mint[['node']], id])]
		ref <- data.table(Frame = (min(data[['Frame']]):(max(data[['Frame']]))))
		dt <- merge(hex_sim, data, by = 'node')[, c('x', 'y', 'Frame')]
		final_data <- merge(dt, ref, all = TRUE, by = 'Frame')
		set(final_data, j = 'x', value = fillVec(final_data[['x']]))
		set(final_data, j = 'y', value = fillVec(final_data[['y']]))
		xy <- TrajFromCoords(final_data[, c('x', 'y', 'Frame')], fps = 2)
		print(i)
		TrajStraightness(xy)
	} else {0}
}, numeric(1)) 
explor_trajs[explor_trajs == 0] <- NA

ggplot(data = data.frame(type = c(rep('Experiments', length(exp_trajs)),
				  rep('Ballistic movement', length(sim_trajs)),
				  rep('Empirical movement', length(model_trajs)),
				  rep('Ballistic exploration', length(explor_trajs))), 
			 y = c(exp_trajs, sim_trajs, model_trajs, explor_trajs)),
       aes(factor(type, levels = c('Experiments', 'Ballistic movement', 'Empirical movement', 'Ballistic exploration')), y))+
	
	geom_boxplot(aes(fill = type), alpha = 0.6, show.legend = FALSE, outlier.shape = NA) +
	geom_jitter(width = 0.15, alpha = 0.25, size = 3)+
	scale_fill_manual('', values = c('gold3', 'brown4', 'mediumpurple', 'grey'))+
	# scale_fill_manual('', values = c('mediumpurple', 'gold3', 'brown4'))+
	scale_y_continuous('Net-to-gross ratio', limits = c(0, 0.8), breaks = seq(0,0.8, 0.2))+
	theme(axis.title.x = element_blank())


path <- '/home/polfer/research/gits/AutomatAnts/results/2024/nest_avoid/'
sim_effs <- global_eff(path)[is.finite(maxt)][, condition := 'Ballistic movement']
path <- '/home/polfer/research/gits/AutomatAnts/results/2024/default/'
model_effs <- global_eff(path)[is.finite(maxt)][, condition := 'Empirical movement']
path <- '/home/polfer/research/gits/AutomatAnts/results/2024/nest_avoid_explor/'
explor_effs <- global_eff(path)[, condition := 'Ballistic exploration']
det_global_eff <- rbindlist(lapply(det, function(i){
	x <- rbindlist(i@food)[['t']]
	# print(min(i@data[['Frame']]))
	mint <- min(x) - min(i@data[['Frame']])
	# 1/(data.table(mint = mint, maxt = max(x)-min(x), t0 = 1/min(i@data[['Frame']]))/2)
	1/(data.table(mint = mint, maxt = max(x)-min(x))/2)
	# 1/(data.table(mint = min(x), maxt = max(x))/2)
}))[, condition := 'Experiments']

# global eff
ggplot(data = melt(rbind(sim_effs[, -'.id'], model_effs[, -'.id'], det_global_eff), id.vars = 'condition'),
       aes(factor(condition, levels = c('Experiments', 
       				 'Ballistic movement', 'Empirical movement')), value, fill = condition)) +
	geom_boxplot(alpha = 0.6, show.legend = FALSE, outlier.shape = NA) +
	geom_jitter(width = 0.15, alpha = 0.25, size = 3, show.legend = FALSE)+
	facet_wrap(~ variable, scales = 'free')+
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')),
		   scales = 'free_y')+
	scale_x_discrete('') + 
	ylab(TeX('$Efficiency (s^{-1})$'))+
	scale_fill_manual('', values = c('mediumpurple', 'gold3', 'brown4'))+
	theme(legend.key.size = unit(30, 'pt'))

ggplot(data = melt(rbind(sim_effs[, -'.id'], explor_effs[, -'.id'], model_effs[, -'.id'], det_global_eff),
		   id.vars = 'condition'),
       aes(factor(condition, levels = c('Experiments', 
       				 'Ballistic movement', 'Empirical movement', 'Ballistic exploration')), value, fill = condition)) +
	geom_boxplot(alpha = 0.6, show.legend = FALSE, outlier.shape = NA) +
	geom_jitter(width = 0.15, alpha = 0.25, size = 3, show.legend = FALSE)+
	facet_wrap(~ variable, scales = 'free')+
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')),
		   scales = 'free_y')+
	scale_x_discrete('') + 
	ylab(TeX('$Efficiency (s^{-1})$'))+
	scale_fill_manual('', values = c( 'gold3', 'brown4', 'mediumpurple', 'grey'))+
	theme(legend.key.size = unit(30, 'pt'))

















foods <- rbindlist(lapply(det, function(i){
	d <- rbindlist(i@food)
	idx <- which.min(d[['t']])
	d[idx]
}))[, node := get_node(x, y), by = seq_along(det)]


ids <- lapply(seq_along(det), function(i){
	id <- setDT(det[[i]]@data)[Frame == foods[i, t] & node == foods[i, node]][['N_ind']]
	data <- det[[i]]@data[N_ind == id]
	Dt <- diff(data[, Frame])
	stops <- which(Dt != 1)
	stop_duration <- Dt[Dt != 1]
	list(data = data, stops = stops, stop_duration = stop_duration, interactions = sum(data[['Crossings']]))
})


###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################


##################################### DET EXPS ##########################################

##############################################################
##################### LR AND SR EXPLORERS DET ################
##############################################################


maxd <- 500 # 650
maxt <- 2400
min_pos <- 15

## calculate velocity, acceleration, straightness, distance travelled, time in arena...
LR_det <- lapply(seq_along(det), function(i){
	
	data <- setDT(det[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	data <- data[N_ind %in% data[, .N, by = 'N_ind'][N > 10, N_ind]]
	get_mot_matrix(data)
	
})
apply(Reduce(`+`, LR_det) / do.call('sum', LR_det), 2, function(i) i / sum(i))

SR_det <- lapply(seq_along(det), function(i){
	
	data <- setDT(det[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[!N_ind %in% inds]
	data <- data[N_ind %in% data[, .N, by = 'N_ind'][N > 10, N_ind]]
	get_mot_matrix(data)
	
})
apply(Reduce(`+`, SR_det) / do.call('sum', SR_det), 2, function(i) i / sum(i))



##############################################################
#################### LR AND SR EXPLORERS NF ################
##############################################################
maxd <- 500 # 650
maxt <- 2400
min_pos <- 15

## calculate velocity, acceleration, straightness, distance travelled, time in arena...
LR_nf <- lapply(seq_along(nf), function(i){
	
	data <- setDT(nf[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	data <- data[N_ind %in% data[, .N, by = 'N_ind'][N > 10, N_ind]]
	get_mot_matrix(data)

})
apply(Reduce(`+`, LR_nf) / do.call('sum', LR_nf), 2, function(i) i / sum(i))

SR_nf <- lapply(seq_along(nf), function(i){
	
	data <- setDT(nf[[i]]@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	set(data, j = 'd', value = dmatrix)
	data <- data[!N_ind %in% inds]
	data <- data[N_ind %in% data[, .N, by = 'N_ind'][N > 10, N_ind]]
	get_mot_matrix(data)
	
})
apply(Reduce(`+`, SR_nf) / do.call('sum', SR_nf), 2, function(i) i / sum(i))
	


	

colnames(exploration_trajs_det) <- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time')






M_det <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
full_det_result <- list()
for(e in seq_along(det)){
	exp <- det[[e]]
	tmax <- min(rbindlist(exp@food)[['t']])
	dt <- setDT(exp@data)[Frame <= tmax, c('N_ind', 'node')]
	Ns <- dt[, .N, N_ind]
	inds <- Ns[N >= min_pos, N_ind]
	result_list <- vector('list', length(inds))
	M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
	for(j in seq_along(inds)){
		df_ind <- exp@data[N_ind == inds[j], node]
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
					M_ind[as.character(result[ii-1]), as.character(result[ii])] <-
						M_ind[as.character(result[ii-1]), as.character(result[ii])]+1
					M_det[as.character(result[ii-1]), as.character(result[ii])] <-
						M_det[as.character(result[ii-1]), as.character(result[ii])]+1
				}
				
				
				# M_ind / sum(M_ind)
				result_list[[j]] <- apply(M_ind / sum(M_ind), 2, function(i) i / sum(i))
				full_det_result <- append(full_det_result, list(M_ind))
			}
		}
	}
	# M
	# apply(M / sum(M), 2, function(i) i / sum(i))
}

