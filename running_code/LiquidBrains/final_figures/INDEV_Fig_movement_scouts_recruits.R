source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/fit_functions.R')
load('~/research/gits/AnTracks/data/nf.RData')
load('~/research/gits/AnTracks/data/det.RData')

det <- lapply(det, class_ids)
nf <- lapply(nf, class_ids)

det_stats <- rbindlist(lapply(det, function(i){
	x <- i
	inds <- x@data[type != 'UNKNOWN', .N, by = 'N_ind'][N > 10, N_ind]
	x@data <- x@data[N_ind %in% inds]
	move_stats(x)
}))

nf_stats <- rbindlist(lapply(nf[-c(1, 2)], function(i){
	x <- i
	inds <- x@data[type != 'UNKNOWN', .N, by = 'N_ind'][N > 10, N_ind]
	x@data <- x@data[N_ind %in% inds]
	move_stats(x)
}))

stats <- rbindlist(list(det = det_stats, nf = nf_stats), idcol = 'exp')
stats <- na.omit(stats)

## 
stats <- stats[d > 1]
stats <- stats[!(type == 'Recruit' & exp == 'det' & d > 500)]

stats_mlt <- data.table(melt(stats[, c('exp', 'v', 'acc', 'd', 't', 'type')],
		  id.vars = c('exp', 'type')))
stats_mlt[['interaction']] <- apply(stats_mlt[, c('exp', 'type')], 1,
				    paste0, collapse = '', sep='')

## transformation of units 
stats_mlt[variable == 'v', value := value / 10] # to cm/s
stats_mlt[variable == 'acc', value := value / 10] # to cm/s^2
stats_mlt[variable == 'd', value := value / 10] # to cm
stats_mlt[variable == 't', value := value / 60] # to min

pairwise.wilcox.test(stats_mlt[variable == 'v', value], 
		     stats_mlt[variable == 'v', interaction], 
		     p.adjust.method = 'bonferroni')

subpanel_1 <- ggplot(data = stats_mlt[variable == 'v'],
		     aes(factor(type, levels = c('Scout', 'Recruit')),
		         value, fill = exp)) + 
	
	geom_boxplot(alpha = 0.6, outlier.shape = 1)+
	annotate('text', x = 0.81, y = 2.25, label = 'a', size = 12) +
	annotate('text', x = 1.19, y = 2.25, label = 'a', size = 12) +
	annotate('text', x = 1.81, y = 2.25, label = 'b', size = 12) +
	annotate('text', x = 2.19, y = 2.25, label = 'c', size = 12) +
	scale_fill_manual('', labels = c('Experimental', 'Control'), 
			  values = c('mediumpurple','gold3'))+
	xlab('') +
	scale_y_continuous(TeX("Velocity ($cm\\cdot s^{-1}$)"))+
	theme(legend.position = c(0.15, 0.8), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.justification = 'center', legend.title = element_blank(),
	      plot.title = element_text(size = 22))

pairwise.wilcox.test(stats_mlt[variable == 'acc', value], 
		     stats_mlt[variable == 'acc', interaction], 
		     p.adjust.method = 'bonferroni')
	
subpanel_2 <- ggplot(data = stats_mlt[variable == 'acc'],
		     aes(factor(type, levels = c('Scout', 'Recruit')),
		         value, fill = exp)) + 
	
	geom_boxplot(alpha = 0.6, outlier.shape = 1, show.legend = FALSE)+
	annotate('text', x = 0.81, y = 3, label = 'a', size = 12) +
	annotate('text', x = 1.19, y = 3, label = 'a', size = 12) +
	annotate('text', x = 1.81, y = 3, label = 'a', size = 12) +
	annotate('text', x = 2.19, y = 3, label = 'a', size = 12) +
	scale_fill_manual('', labels = c('Experimental', 'Control'), 
			  values = c('mediumpurple','gold3'))+
	xlab('') +
	scale_y_continuous(TeX("Acceleration ($cm\\cdot s^{-2}$)"))+
	theme(legend.position = c(0.15, 0.8), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.justification = 'center', legend.title = element_blank(),
	      plot.title = element_text(size = 22))+
	coord_cartesian(ylim = c(0, 5))

pairwise.wilcox.test(stats_mlt[variable == 'd', value], 
		     stats_mlt[variable == 'd', interaction], 
		     p.adjust.method = 'bonferroni')

subpanel_3 <- ggplot(data = stats_mlt[variable == 'd'],
		     aes(factor(type, levels = c('Scout', 'Recruit')),
		         value, fill = exp)) + 
	
	geom_boxplot(alpha = 0.6, outlier.shape = 1, show.legend = FALSE)+
	annotate('text', x = 0.81, y = 120, label = 'a', size = 12) +
	annotate('text', x = 1.19, y = 120, label = 'b', size = 12) +
	annotate('text', x = 1.81, y = 120, label = 'c', size = 12) +
	annotate('text', x = 2.19, y = 120, label = 'd', size = 12) +
	scale_fill_manual('', labels = c('Experimental', 'Control'), 
			  values = c('mediumpurple','gold3'))+
	xlab('') +
	scale_y_continuous(TeX('Maximum distance (cm)'))+
	theme(legend.position = c(0.15, 0.8), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.justification = 'center', legend.title = element_blank(),
	      plot.title = element_text(size = 22))+
	coord_cartesian(ylim = c(0, 125))


pairwise.wilcox.test(stats_mlt[variable == 't', value], 
		     stats_mlt[variable == 't', interaction], 
		     p.adjust.method = 'bonferroni')

subpanel_4 <- ggplot(data = stats_mlt[variable == 't'],
		     aes(factor(type, levels = c('Scout', 'Recruit')),
		         value, fill = exp)) + 
	
	geom_boxplot(alpha = 0.6, outlier.shape = 1, show.legend = FALSE)+
	annotate('text', x = 0.81, y = 15, label = 'a', size = 12) +
	annotate('text', x = 1.19, y = 15, label = 'b', size = 12) +
	annotate('text', x = 1.81, y = 15, label = 'c', size = 12) +
	annotate('text', x = 2.19, y = 15, label = 'c', size = 12) +
	scale_fill_manual('', labels = c('Experimental', 'Control'), 
			  values = c('mediumpurple','gold3'))+
	xlab('') +
	scale_y_continuous(TeX('Time (min)'))+
	theme(legend.position = c(0.15, 0.8), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.justification = 'center', legend.title = element_blank(),
	      plot.title = element_text(size = 22))+
	coord_cartesian(ylim = c(0, 20))

panel_a <- grid.arrange(subpanel_1 + theme(axis.line.x = element_blank(),
				axis.text.x = element_blank(),
				axis.ticks.x = element_blank()),
	     subpanel_2+ theme(axis.line.x = element_blank(),
	     		  axis.text.x = element_blank(),
	     		  axis.ticks.x = element_blank()),
	     subpanel_3, subpanel_4)

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

pls <- rbindlist(lapply(nf[-c(1, 2)], function(i){
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
	annotate('text', x = 0.81, y = 6, label = 'a', size = 12) +
	annotate('text', x = 1.19, y = 6, label = 'b', size = 12) +
	annotate('text', x = 1.81, y = 6, label = 'c', size = 12) +
	annotate('text', x = 2.19, y = 6, label = 'd', size = 12) +
	scale_fill_manual('', labels = c('Experimental', 'Control'), 
			  values = c('mediumpurple','gold3'))+
	xlab('') +
	scale_y_continuous(TeX('Directional persistence (steps)'))+
	coord_cartesian(ylim = c(0, 7))


det_straight <- rbindlist(lapply(det, get_straightness), idcol = 'exp')
nf_straight <- rbindlist(lapply(nf[-c(1, 2)], get_straightness), idcol = 'exp')

straightness <- rbindlist(list(det = det_straight, nf = nf_straight), idcol = 'condition')
straightness <- na.omit(straightness)


ggplot(data = straightness[, .(l = round(l/100)*10,
			       s = s, type = type, condition = condition)][, .(s = mean(s),
			       						sd = sd(s)),
			       					    by = c('l', 'type', 'condition')],
       aes(l, s, color = factor(condition, labels = c('Experimental', 'Control')),
           shape = factor(type, levels = c('Scout', 'Recruit'))), alpha = 0.5) + 
	geom_pointrange(aes(ymin = s - sd, ymax = s+sd),
			size = 1.5, alpha = 0.5) +
	# geom_path(data = pred_data, linewidth = 1.5)+
	scale_y_continuous('<Net-to-Gross ratio>', limits = c(0, 1), breaks = seq(0, 1, 0.2))+
	scale_x_continuous('Travelled distance (cm)', limits = c(0, 500))+
	scale_color_manual('', values = c('mediumpurple', 'gold3'), 
			   labels = c('Experimental', 'Control'))+
	scale_shape_manual('', values = c(21, 24))+
	theme(legend.title = element_blank(), 
	      legend.position = c(3/4, 3/4),
	      legend.box.background = element_rect(),
	      legend.background = element_blank(), plot.title = element_text(size = 22))+
	guides(shape = guide_legend(override.aes = list(shape = c(19, 17))))






setMethod('get_straightness', 'Experiment', function(.Object, ...){
	require(trajr)
	.Object <- class_ids(.Object, ...)
	data <- .Object@data[type != 'UNKNOWN']
	ids <- data[, .N, by = N_ind][N > 10, N_ind]
	
	rbindlist(lapply(seq_along(ids), function(i){
		sbst <- data[N_ind == ids[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		
		Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
		maxL <- max(Ls)
		if(maxL >= 100){
			idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
		} else {
			idxs <- nrow(xy)
		}
		straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
		df <- data.frame(l = Ls[idxs], s = straightness)
		df[['type']] <- sbst[, unique(type)]
		df
	}), idcol = 'ind')
})


data <- det[[1]]@data[type != 'UNKNOWN']
data <- data[N_ind %in% data[, .N, by = 'N_ind'][N > 10, N_ind]]
pls <- lapply(unique(data[['N_ind']]), function(i){
	a <- data[N_ind == i][c(1, diff(node)) != 0, node]
	if(length(a) > 3){
		dirs <- vapply(3:length(a), function(i){
			x <- list(hex[hex$node == a[[i-2]], c('x', 'y')],
				  hex[hex$node == a[[i-1]], c('x', 'y')],
				  hex[hex$node == a[[i]], c('x', 'y')])
			get_direction(x)
		}, numeric(1))
	}
	mean(find_pl(dirs))
})


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

ggplot(data = pls_det, aes(type, pl)) + geom_boxplot()


pls <- rbindlist(lapply(nf[-c(1, 2)], function(i){
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

ggplot(data = pls, aes(type, pl)) + geom_boxplot()


nds <- smple_data[['node']]


megatime <- sapply(1:1000, function(x){
	at <- Sys.time()
	nds_ind <- sapply(2:length(nds), function(i){
		if(nds[i] != nds[i-1]){
			return(nds[i])
		}
	})
	Sys.time() -at
})

gigatime <- c()
for(x in 1:1000){
	at <- Sys.time()
	data.table(nds = nds)[c(1, diff(nds)) != 0, nds]
	gigatime <- c(gigatime, Sys.time() -at)
}


nds_ind <- append(nds[1], nds_ind)
nds_ind <- nds_ind[lapply(nds_ind, length) > 0]




dirs_det <- rbindlist(lapply(det, function(i){
	data <- i@data[type != 'UNKNOWN']
	rbindlist(lapply(unique(data[['N_ind']]), function(x){
		a <- data[N_ind == x][c(1, diff(node)) != 0, node]
		m <- as.matrix(hex[a, c('x', 'y')])
		if(nrow(m) > 3){
			dirs <- get_directions_cpp(m)
			data.frame(dir = dirs, type = data[N_ind == x, unique(type)])
		} else {
			data.frame(dir = NA, type = data[N_ind == x, unique(type)])
		}
	}), idcol = 'N_ind')
}), idcol = 'exp')
dirs_det <- na.omit(dirs_det)
apply(table(dirs_det[, c('dir', 'type')]), 2, function(i) i / sum(i))

dirs_nf <- rbindlist(lapply(nf[-c(1, 2)], function(i){
	data <- i@data[type != 'UNKNOWN']
	rbindlist(lapply(unique(data[['N_ind']]), function(x){
		a <- data[N_ind == x][c(1, diff(node)) != 0, node]
		m <- as.matrix(hex[a, c('x', 'y')])
		if(nrow(m) > 3){
			dirs <- get_directions_cpp(m)
			data.frame(dir = dirs, type = data[N_ind == x, unique(type)])
		} else {
			data.frame(dir = NA, type = data[N_ind == x, unique(type)])
		}
	}), idcol = 'N_ind')
}), idcol = 'exp')
dirs_nf <- na.omit(dirs_nf)
apply(table(dirs_nf[, c('dir', 'type')]), 2, function(i) i / sum(i))













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

