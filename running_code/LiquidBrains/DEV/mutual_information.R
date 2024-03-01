#### LIBRARIES, DATA AND GENERIC FUNCTIONS ####
library(dtw)
library(dtwclust)
library(latex2exp)
library(arrow)
library(infotheo)
source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
load("/home/polfer/research/gits/AutomatAnts/data/zSeries.RData")

load("/home/polfer/research/gits/AutomatAnts/data/zSeries.RData")
load('~/research/gits/AnTracks/data/det.RData')
# load('~/research/gits/AnTracks/data/sto.RData')

rectify_nodes <- function(joly){
	apply(do.call('rbind', food_sim[['node']]), 1,
	      function(i) paste0('(', paste0(i, collapse = ', '), ')'))
}

parse_nodes <- function(nodes){
	unlist(strsplit(nodes, ';'))
}

parse_ids <- function(ids){
	as.integer(unlist(strsplit(ids, ',')))
}

MI <- function(x, tw = 2400, shift = 300){
	t0 <- tw %/% 2
	mint <- min(x[N > 0, Frame])
	sq <- seq(t0, max(x[['Frame']]) - tw, shift)
	result <- vector('list', length = length(sq))
	names(result) <- sq
	for(i in seq_along(sq)){
		upper_lim <- sq[i]+ t0 -1
		if(upper_lim > mint){
			current_sequence <- seq(sq[i] - t0, sq[i] + t0 -1, 1)
			dt <- x[Frame %in% current_sequence, .(node = parse_nodes(pos)), 
				by = 'Frame'][, .(N = .N), by = c('Frame', 'node')]
			M <- dt[CJ(Frame = current_sequence, node = node, unique = T),
				on=. (Frame, node)][is.na(N), N:=-1][N > 0, N:=1]
			m <- dcast(M, Frame ~ node, value.var = 'N')
			m$Frame <- NULL
			# result[[i]] <- mean(mutinformation(m))
			result[[i]] <- ncol(m)[-1]
		}
	}
	result
}

tMI <- function(x, tw = 2400, shift = 300){
	t0 <- tw %/% 2
	mint <- min(x[N > 0, Frame])
	sq <- seq(t0, max(x[['Frame']]) - tw, shift)
	result <- parallel::mclapply(seq_along(sq), function(i){
		upper_lim <- sq[i]+ t0 -1
		if(upper_lim > mint){
			current_sequence <- seq(sq[i] - t0, sq[i] + t0 -1, 1)
			dt <- x[Frame %in% current_sequence, .(node = parse_nodes(pos)), 
				by = 'Frame'][, .(N = .N), by = c('Frame', 'node')]
			M <- dt[CJ(Frame = current_sequence, node = node, unique = T),
				on=. (Frame, node)][is.na(N), N:=-1][N > 0, N:=1]
			m <- dcast(M, Frame ~ node, value.var = 'N')
			m$Frame <- NULL
			mean(mutinformation(m))
		} else {
			NA
		}
	}, mc.cores = 8L)
	names(result) <- sq
	result[lapply(result, is.na) == FALSE]
}

path <- '/home/polfer/research/gits/AutomatAnts/results/2024/beta_0.66/'
# path <- '/home/polfer/research/gits/AutomatAnts/results/2024/default/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('beta=0.66_\\d{1,2}', files))))
# files <- unique(unlist(regmatches(files, gregexpr('default_\\d{1,2}', files))))

# path <- '/home/polfer/research/gits/AutomatAnts/results/2024/beta_0.66/'
# files <- list.files(path)
# files <- unique(unlist(regmatches(files, gregexpr('beta=0.66_\\d{1,2}', files))))
pos <- data.table(read_parquet(paste0(path, files[4], '_positions.parquet')))




data <- data.table(read_parquet(paste0(path, files[1], '_data.parquet')))
t0 <- Sys.time()
r <- MI(data)
Sys.time() - t0

t0 <- Sys.time()
r2 <- MI(data, tw = 600, shift = 300)
Sys.time() - t0

t0 <- Sys.time()
r2 <- MI(data, tw = 200, shift = 200)
Sys.time() - t0

data <- data.table(read_parquet(paste0(path, files[2], '_data.parquet')))
t0 <- Sys.time()
r3 <- MI(data, tw = 600, shift = 200)
Sys.time() - t0

ggplot(data = rbindlist(lapply(r2, function(i) data.table(v = i)), idcol = T), aes(.id, v)) +
	geom_point(shape = 21, color = 'black', fill = NA) + geom_path()+
	scale_x_continuous('Time (min)', breaks = c(0, 50*120, 100*120, 150*120),
			   labels = c(0, 50, 100, 150))+
	scale_y_continuous('Average mutual information')

plot(rbindlist(lapply(r, function(i) data.table(v = i)), idcol = T), type = 'l')
plot(rbindlist(lapply(r2, function(i) data.table(v = i)), idcol = T), type = 'l')
plot(rbindlist(lapply(r3, function(i) data.table(v = i)), idcol = T), type = 'l')


plot(rbindlist(lapply(r, function(i) data.table(v = i)), idcol = T))
lines(rbindlist(lapply(r, function(i) data.table(v = i)), idcol = T))
	
path <- '/home/polfer/research/gits/AutomatAnts/results/2024/default/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('default_\\d{1,2}', files))))

data <- data.table(read_parquet(paste0(path, files[1], '_data.parquet')))

m$Frame <- NULL


dt <- data[, .(node = parse_nodes(pos), id = parse_ids(id_out)), by = 'Frame'][!is.na(node)]




global_eff <- function(path){
	files <- list.files(path)[grepl('food', list.files(path))]
	results <- lapply(files, function(i){
		food <- data.table(read_parquet(paste0(path, i)))
		1/data.table(mint = min(food[['t']]), maxt = max(food[['t']]))
	})
	rbindlist(results, idcol = TRUE)
}



eff <- function(path, filename, minpieces = 0){
	food <- data.table(read_parquet(paste0(path, filename ,'_food.parquet')))
	food <- unlist(food[!is.na(t), .(t = range(t))], use.names = FALSE) *2
	if(length(food)){
		dt <- data.table(read_parquet(paste0(path, filename, '_data.parquet')))
		
		
		filterdt <- unique(dt[Frame <= food[2], .(id = parse_ids(id_out)), by = 'Frame'])
		M <- dcast.data.table(data = filterdt, Frame ~ id, value.var = 'id', )
		.tmp <- rbindlist(lapply(2:ncol(M), function(i){
			idx <- which(!is.na(M[[i]]))
			data.frame(id = unique(M[[i]][idx]),Frame = M[[1]][idx], d=c(1, diff(idx)))
		}))
		.tmp[d > 1, 'd'] <- 0
		.tmp[, interval := cumsum(c(TRUE, diff(d) != 0)), by = id]
		.tmp[, interval := ifelse(interval %% 2 == 0, interval + 1, interval), by = id]
		.tmp_result <- .tmp[, .(start_frame = min(Frame), end_frame = max(Frame)), by = .(id, interval)]
		
		.tp1 <- .tmp_result[start_frame <= food[1]][, end_frame := ifelse(end_frame < food[1], end_frame, food[1])]
		tp1_value <- sum(.tp1[, end_frame] - .tp1[, start_frame])/2
		.tp2 <- .tmp_result[start_frame <= food[2] & end_frame >= food[1]][, start_frame := ifelse(start_frame < food[1], food[1], start_frame)]
		.tp2 <- .tp2[, end_frame := ifelse(end_frame < food[2], end_frame, food[2])]
		tp2_value <- sum(.tp2[, end_frame] - .tp2[, start_frame])/2
		data.table(tp1 = tp1_value, tp2 = tp2_value)
	}
}

global_sim_effs <- global_eff(path)[is.finite(maxt)][, condition := 'Simulations']
det_global_eff <- rbindlist(lapply(det, function(i){
	x <- rbindlist(i@food)[['t']]
	1/(data.table(mint = min(x), maxt = max(x))/2)
	}))[, condition := 'Experiments']

ggplot(data = melt(rbind(global_sim_effs[, -'.id'], det_global_eff), id.vars = 'condition'),
       aes(condition, value, fill = condition)) + geom_boxplot(alpha = 0.6) +
	facet_wrap(~ variable, scales = 'free')+
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')),
		   scales = 'free_y')+
	scale_x_discrete('') + 
	ylab(TeX('$Efficiency (s^{-1})$'))+
	scale_fill_manual('', values = c('mediumpurple', 'gold3'))+
	theme(legend.key.size = unit(30, 'pt'))

collective_sim_effs <- lapply(seq_along(files), function(i){
	do.call('gc', args = list(verbose = FALSE))
	eff(path, files[i])
})

det_tps <- lapply(det, get_eff)
sto_tps <- lapply(sto, get_eff)
tps <- cbind(condition = 'exp', rbindlist(l = list(rbindlist(det_tps),rbindlist(sto_tps))))

coll_dt <- rbindlist(collective_sim_effs)

coll_dt[['condition']] <- 'Simulations'


ggplot(data = melt(rbindlist(list(coll_dt, cbind(condition = 'Experiments', rbindlist(det_tps))), use.names = TRUE), id.vars = 'condition'), 
       aes(condition, 1/value, fill = condition))+
	geom_boxplot(alpha = 0.6)+ 
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')),
		   scales = 'free_y')+
	scale_x_discrete('') + 
	ylab(TeX('$Efficiency (s^{-1})$'))+
	scale_fill_manual('', values = c('mediumpurple', 'gold3'))+
	theme(legend.key.size = unit(30, 'pt'))

ggplot(data = rbind(melt(rbind(global_sim_effs[, -'.id'], det_global_eff), id.vars = 'condition'),
       data.table(melt(rbindlist(list(coll_dt, 
				  cbind(condition = 'Experiments',
				        rbindlist(det_tps))), use.names = TRUE), 
		   id.vars = 'condition'))[, value := 1/value]), 
       aes(condition, value, fill = condition))+
	geom_boxplot(alpha = 0.6) +
	facet_wrap(~ factor(variable, labels = c('Exploration (collective)', 'Exploitation (collective)', 
						 'Exploration (individual)', 'Exploitation (individual)')),
		   scales = 'free_y')+
	scale_x_discrete('') + 
	ylab(TeX('$Efficiency (s^{-1})$'))+
	scale_fill_manual('', values = c('mediumpurple', 'gold3'))+
	theme(legend.key.size = unit(30, 'pt'))

# ggplot(data = melt(rbindlist(list(coll_dt, tps), use.names = TRUE), id.vars = 'condition'), 
#        aes(condition, 1/value, fill = condition))+
# 	geom_boxplot(alpha = 0.6)+ 
# 	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')),
# 		   scales = 'free_y')+
# 	scale_x_discrete('') + 
# 	ylab(TeX('$Efficiency (s^{-1})$'))+
# 	scale_fill_manual('', values = c('mediumpurple', 'gold3', 'deepskyblue4'))+
# 	theme(legend.key.size = unit(30, 'pt'))




rectify_nodes <- function(n){
	print(n)
	n <- rbindlist(lapply(n[['node']], function(i) data.table(i)))
	print(n)
	apply(n, 1,
	      function(i) paste0('(', paste0(i, collapse = ', '), ')'))
}

revert_node <- function(n){
	n <- do.call('rbind', n)
	apply(n, 1, function(i) paste0('(', paste0(i, collapse = ', '), ')'))
}


afunc <- function(path, filename, filter = NA){
	do.call('gc', args = list(verbose = FALSE))
	food_sim <- data.table(read_parquet(paste0(path, filename,'_food.parquet')))
	food_sim <- data.table(node = revert_node(food_sim[['node']]), t = food_sim[['t']])
	# food_sim <- food_sim[, .(node = revert_node(node), t = t)]
	dt <- data.table(read_parquet(paste0(path, filename,'_data.parquet')))[, c('Frame', 'pos', 'id_out')]
	if(is.numeric(filter)){
		dt <- dt[Frame <= filter]
	}
	dt_complete <- unique(dt[, .(node = parse_nodes(pos), id = parse_ids(id_out)), by = 'Frame'][!is.na(node)])
	
	ids_filter <- dt_complete[node %in% food_sim[['node']], .(id = unique(id))][['id']]
	dt_filtered <- dt_complete[id %in% ids_filter]
	
	tmp_1 <- dt_filtered[, .(mint = min(Frame)), by = c('node', 'id')][node %in% food_sim[['node']]]
	tmp_2 <- dt_filtered[, .(Frame = Frame), by = c('node','id')][node == '(0, 22)']
	
	times_to_food(tmp_1, tmp_2)
}

times_to_food <- function(x, y){
	
	unlist(lapply(1:nrow(x), function(i){
		maxt <- x[i, mint]
		f <- max(y[id == x[i, id] & Frame < maxt & node == '(0, 22)', Frame])
		maxt - f
	}))
}

t0 <- Sys.time()
times_simulation <- c()
for(i in files){
	times_simulation <- c(times_simulation, afunc(path, i))
}
Sys.time() - t0

t0 <- Sys.time()
times_simulation_tp1 <- c()
for(i in files){
	times_simulation_tp1 <- c(times_simulation_tp1, afunc(path, i, filter = 2990))
}
Sys.time() - t0

t0 <- Sys.time()
times_simulation_tp2 <- c()
for(i in files){
	times_simulation_tp2 <- c(times_simulation_tp2, afunc(path, i, filter = 6525))
}
Sys.time() - t0







get_times <- function(exp){
	rbindlist(exp@food)[['t']]
}

get_foodnodes <- function(exp, ref = hex[hex$y > 1000, ]){
	get_node(rbindlist(exp@food)[, 1:2], xy = ref)
}

get_ids <- function(exp, ...){
	if(!is.data.table(exp@data)){
		setDT(exp@data)
	}
	fnds <- get_foodnodes(exp, ...)
	exp@data[node %in% fnds, .(id = unique(N_ind))][['id']]
}
tracklets_to_food <- function(exp, filter = NA){
	nest_nodes <- c(662, 608, 634, 635)
	fnds <- get_foodnodes(exp)
	dt <- setDT(exp@data)
	ids <- get_ids(exp)
	sbst_exp <- dt[N_ind %in% ids]
	
	if(is.numeric(filter)){
		sbst_exp <- sbst_exp[Frame <= filter]
	}
	cross_ids <- sbst_exp[, .(x = length(unique(Crossings)) == 1), by = 'N_ind']
	filtered_ids <- sbst_exp[N_ind %in% cross_ids[x == TRUE, N_ind]]
	result <- c()
	for(i in ids){
		visited_nodes <- filtered_ids[N_ind == i & node %in% fnds,
					      .(visited_nodes = unique(node))][['visited_nodes']]
		if(length(visited_nodes) > 0){
			tmp_origin <- filtered_ids[N_ind == i, .(mint = min(Frame)), by = node][node %in% nest_nodes]
			tmp_destiny <- filtered_ids[N_ind == i, .(mint = min(Frame)), by = node][node %in% visited_nodes]
			tmp <- tmp_destiny[which.min(mint), mint] - tmp_origin[which.min(mint), mint]
			result <- c(result, tmp)
		}
	}
	result[result > 0]
}

apply(do.call('rbind', lapply(det, function(i) range(rbindlist(i@food)[['t']]))), 2, max)
tdet_tp1 <- unlist(lapply(det, tracklets_to_food, filter = 2990))
tdet_tp2 <- unlist(lapply(det, tracklets_to_food, filter = 6525))
tdet <- unlist(lapply(det, tracklets_to_food))

# tsto <- unlist(lapply(sto, tracklets_to_food, filter = FALSE))

# tdet_noFil <- unlist(lapply(det, tracklets_to_food, filter = FALSE))
# tsto_noFil <- unlist(lapply(sto, tracklets_to_food, filter = FALSE))

a <- ggplot(data = data.frame(x = c(times_simulation[is.finite(times_simulation)] / 120,
			       tdet/120),
			 type = c(rep('Simulation', sum(is.finite(times_simulation))),
			 	 rep('Determinist', length(tdet)))), 
       aes(x, fill = type)) +
	geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 35, 0.5),
		       color = 'black', alpha = 0.5, position = 'identity') + 
	geom_density(alpha = 0.5)+
	scale_x_continuous('Time (min)', breaks = seq(0, 10, 2.5),
			   limits = c(0, 10))+

	geom_vline(xintercept = c(median(times_simulation[is.finite(times_simulation)] / 120),
				  median(tdet/120)),
		   linewidth = 1, linetype = 2, color = c('mediumpurple', 'gold3'))+
	scale_y_continuous('Density', breaks = seq(0, 10, 0.15)) +
	scale_fill_manual('',values = c('gold3', 'mediumpurple'))+
	annotation_custom(
		ggplotGrob(ggplot(data = data.frame(x = c(times_simulation[is.finite(times_simulation)] / 120,
							  tdet/120),
							  type = c(rep('Simulation', sum(is.finite(times_simulation))),
							  	 rep('Determinist', length(tdet)))))+

			   	geom_violin(aes(type, x, fill = type), trim = FALSE,
			   		    alpha = 0.5, show.legend = FALSE)+
			   	theme(
			   		axis.ticks.x = element_blank())+xlab('')+
			   	scale_y_continuous('Time (min)', breaks = seq(0, 150, 15))+
			   	scale_fill_manual(values = c('gold3', 'mediumpurple'))),
		xmin = 4, xmax = 7.5, ymin = 0.2, ymax = 0.65)+
	geom_rect(aes(xmin = 4, xmax = 7.5, ymin = 0.2, ymax = 0.69), fill = NA, color = 'black',
		  linewidth = 0.5)+
	theme(legend.position = c(0.85, 0.5), 
	      legend.background = element_rect(fill = NA, colour = 'black'),
	      legend.title = element_blank())


b <- ggplot(data = data.frame(x = c(times_simulation_tp1[is.finite(times_simulation_tp1)] / 120,
			       tdet_tp1/120),
			 type = c(rep('Simulation', sum(is.finite(times_simulation_tp1))),
			 	 rep('Determinist', length(tdet_tp1)))), 
       aes(x, fill = type)) +
	geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 35, 0.5),
		       color = 'black', alpha = 0.5, position = 'identity') + 
	geom_density(alpha = 0.5)+
	scale_x_continuous('Time (min)', breaks = seq(0, 10, 2.5),
			   limits = c(0, 10))+
	
	geom_vline(xintercept = c(median(times_simulation_tp1[is.finite(times_simulation_tp1)] / 120),
				  median(tdet_tp1/120)),
		   linewidth = 1, linetype = 2, color = c('mediumpurple', 'gold3'))+
	scale_y_continuous('Density', breaks = seq(0, 10, 0.15)) +
	scale_fill_manual('',values = c('gold3', 'mediumpurple'))+
	annotation_custom(
		ggplotGrob(ggplot(data = data.frame(x = c(times_simulation_tp1[is.finite(times_simulation_tp1)] / 120,
							  tdet_tp1/120),
							  type = c(rep('Simulation', sum(is.finite(times_simulation_tp1))),
							  	 rep('Determinist', length(tdet_tp1)))))+
			   	
			   	geom_violin(aes(type, x, fill = type), trim = FALSE,
			   		    alpha = 0.5, show.legend = FALSE)+
			   	theme(
			   		axis.ticks.x = element_blank())+xlab('')+
			   	scale_y_continuous('Time (min)', breaks = seq(0, 150, 2.5))+
			   	scale_fill_manual(values = c('gold3', 'mediumpurple'))),
		xmin = 4, xmax = 7.5, ymin = 0.2, ymax = 0.65)+
	geom_rect(aes(xmin = 4, xmax = 7.5, ymin = 0.2, ymax = 0.69), fill = NA, color = 'black',
		  linewidth = 0.5)+
	theme(legend.position = c(0.85, 0.5), 
	      legend.background = element_rect(fill = NA, colour = 'black'),
	      legend.title = element_blank())


c <- ggplot(data = data.frame(x = c(times_simulation_tp2[is.finite(times_simulation_tp2)] / 120,
			       tdet_tp2/120),
			 type = c(rep('Simulation', sum(is.finite(times_simulation_tp2))),
			 	 rep('Determinist', length(tdet_tp2)))), 
       aes(x, fill = type)) +
	geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 35, 0.5),
		       color = 'black', alpha = 0.5, position = 'identity') + 
	geom_density(alpha = 0.5)+
	scale_x_continuous('Time (min)', breaks = seq(0, 10, 2.5),
			   limits = c(0, 10))+
	
	geom_vline(xintercept = c(median(times_simulation_tp2[is.finite(times_simulation_tp2)] / 120),
				  median(tdet_tp2/120)),
		   linewidth = 1, linetype = 2, color = c('mediumpurple', 'gold3'))+
	scale_y_continuous('Density', breaks = seq(0, 10, 0.15)) +
	scale_fill_manual('',values = c('gold3', 'mediumpurple'))+
	annotation_custom(
		ggplotGrob(ggplot(data = data.frame(x = c(times_simulation_tp2[is.finite(times_simulation_tp2)] / 120,
							  tdet_tp2/120),
							  type = c(rep('Simulation', sum(is.finite(times_simulation_tp2))),
							  	 rep('Determinist', length(tdet_tp2)))))+
			   	
			   	geom_violin(aes(type, x, fill = type), trim = FALSE,
			   		    alpha = 0.5, show.legend = FALSE)+
			   	theme(
			   		axis.ticks.x = element_blank())+xlab('')+
			   	scale_y_continuous('Time (min)', breaks = seq(0, 150, 5))+
			   	scale_fill_manual(values = c('gold3', 'mediumpurple'))),
		xmin = 4, xmax = 7.5, ymin = 0.2, ymax = 0.65)+
	geom_rect(aes(xmin = 4, xmax = 7.5, ymin = 0.2, ymax = 0.69), fill = NA, color = 'black',
		  linewidth = 0.5)+
	theme(legend.position = c(0.85, 0.5), 
	      legend.background = element_rect(fill = NA, colour = 'black'),
	      legend.title = element_blank())


png('/home/polfer/research/gits/AnTracks/plots/distrib_times_splitted.png', width = 6000, height = 6000, res = 300)
ggarrange(b, c, a, labels = c('TP1', 'TP2', 'FULL'), ncol = 1)
dev.off()
# ggplot(data = data.frame(x = c(times_simulation[is.finite(times_simulation)] / 120,
# 			       c(tdet_noFil, tsto_noFil)/120),
# 			 type = c(rep('Simulation', sum(is.finite(times_simulation))),
# 			 	 rep('Experiment', length(c(tdet_noFil, tsto_noFil))))), 
#        aes(x, fill = type)) +
# 	geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 35, 0.15),
# 		       color = 'black', alpha = 0.5, position = 'identity') + 
# 	geom_density(alpha = 0.5)+
# 	scale_x_continuous('Time (min)', breaks = seq(0, 10, 2.5),
# 			   limits = c(0, 10))+
# 	# geom_vline(xintercept = c(median(times_simulation[is.finite(times_simulation)] / 120),
# 	# 			  median( c(tdet_noFil, tsto_noFil)/120)),
# 	# 	   linewidth = 1.2, linetype = 2)+
# 	geom_vline(xintercept = c(median(times_simulation[is.finite(times_simulation)] / 120),
# 				  median( c(tdet_noFil, tsto_noFil)/120)),
# 		   linewidth = 1, linetype = 2, color = c('mediumpurple', 'gold3'))+
# 	scale_y_continuous('Density', breaks = seq(0, 10, 0.15)) +
# 	scale_fill_manual('',values = c('gold3', 'mediumpurple'))+
# 	annotation_custom(
# 		ggplotGrob(ggplot(data = data.frame(x = c(times_simulation[is.finite(times_simulation)] / 120,
# 							  c(tdet_noFil, tsto_noFil)/120),
# 							  type = c(rep('Simulation', sum(is.finite(times_simulation))),
# 							  	 rep('Experiment', length(c(tdet_noFil, tsto_noFil))))))+
# 			   	# geom_boxplot(aes(type, x, fill = type), width = 0.5,
# 			   	# 	     alpha = 0.5,outlier.shape = 1, show.legend = FALSE)+
# 			   	geom_violin(aes(type, x, fill = type), trim = FALSE,
# 			   		    alpha = 0.5, show.legend = FALSE)+
# 			   	theme(
# 			   		axis.ticks.x = element_blank())+xlab('')+
# 			   	scale_y_continuous('Time (min)', breaks = seq(0, 100, 10))+
# 			   	scale_fill_manual(values = c('gold3', 'mediumpurple'))),
# 		xmin = 4, xmax = 7.5, ymin = 0.2, ymax = 0.65)+
# 	geom_rect(aes(xmin = 4, xmax = 7.5, ymin = 0.2, ymax = 0.69), fill = NA, color = 'black',
# 		  linewidth = 0.5)+
# 	theme(legend.position = c(0.85, 0.5), 
# 	      legend.background = element_rect(fill = NA, colour = 'black'),
# 	      legend.title = element_blank())











ndtw <- function(x, y, ...) {
	dtw::dtw(x, y, step.pattern = asymmetric,
		 distance.only = TRUE, ...)$normalizedDistance
}

process_data <- function(path, 
			 vars = c('Frame', 'N')){
	ref <- data.table(Frame = 1:21600)
	t <- movingAverage(ref[['Frame']], 60, 0)
	
	files <- list.files(path)
	files <- files[!grepl('food', files) & !grepl('position', files) &
		       	!grepl('data', files) & !grepl('keys', files)]
	
	a <- lapply(files, function(i){
		
		y <- merge(ref, data.table(read_parquet(paste0(path, i)))[, ..vars], all = TRUE, by = 'Frame')
		y[['N']] <- fillVec(y[['N']])
		n <- movingAverage(y[['N']], 60L, 0L)
		z <- zscore(n)
		data.table(Frame = t, N = n, Z = z)
		
	})
	names(a) <- files
	a
}

plot_clusters <- function(data_list, clusters){
	kprop <- 100* (table(clusters) / length(clusters))
	n <- names(kprop)
	xPlot <- rbindlist(lapply(seq_along(data_list), function(i){
		data.table(Frame = data_list[[i]][['Frame']], N = data_list[[i]][['N']], k = clusters[i], exp = i)
	}))
	
	ggplot(data = xPlot, aes(Frame, N)) + geom_path(aes(group = exp), color = 'grey70', alpha = 0.75)+
		geom_path(data = xPlot[, .(N = mean(N), exp = exp), by = c('Frame', 'k')], linewidth = 1)+
		facet_wrap(~factor(k, levels = c('norm', 'flat', 'late', 'low'),
				   labels = c(paste0('Experimental pattern (', 
				   		  ifelse('norm' %in% n, kprop[['norm']], 0),
				   		  '%)'), 
				   	   paste0('Plateau (', 
				   	          ifelse('flat' %in% n, kprop[['flat']], 0),
				   	          '%)'), 
				   	   paste0('Late start (', 
				   	          ifelse('late' %in% n, kprop[['late']], 0),
				   	          '%)'), 
				   	   paste0('Low activity (', 
				   	          ifelse('low' %in% n, kprop[['low']], 0),
				   	          '%)'))), 
			   nrow = 1) +
		ylab('Activity (number of ants in arena)') + 
		scale_x_continuous('Time (min)', breaks = seq(0, 21600, 30 * 120), labels = seq(0, 180, 30))+
		theme(strip.text = element_text(size = 16, margin = margin(t = 5, b = 5, unit = 'pt')),
		      aspect.ratio = 0.75)
}



# path <- '/home/polfer/research/gits/AutomatAnts/results/2024/default/'
path <- '/home/polfer/research/gits/AutomatAnts/results/2024/beta_0.5/'
path <- '/home/polfer/research/gits/AutomatAnts/results/2024/beta_0.66/'
path <- '/home/polfer/research/gits/AutomatAnts/results/2024/beta_0.7/'
path <- '/home/polfer/research/gits/AutomatAnts/results/2024/beta_0.8/'

path <- '/home/polfer/research/gits/AutomatAnts/results/2024/Jij_0.01/'
# path <- '/home/polfer/research/gits/AutomatAnts/results/2024/beta_1'

knames <- c('norm', 'late', 'flat', 'low')
data <- process_data(path)
distances <- rbindlist(lapply(data, function(i){
	data.table(t(sapply(zSeries, function(x){
		ndtw(i[['Z']], x)
	})))
}))
colnames(distances) <- knames
k <- knames[apply(distances, 1, which.min)]

plot_clusters(data, k)



files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('Jij_0.01_\\d{1,2}', files))))
full_data <- rbindlist(lapply(files, function(i){
	data.table(read_parquet(paste0(path, i,'.parquet')))
}), idcol = TRUE)[order(Frame)]
full_data_avg <- full_data[, .(N = mean(N)), by = 'Frame']
filtered_avg <- full_data[.id %in% which(k == 'norm'), .(N = mean(N)), by = 'Frame']


det_avg <- rbindlist(lapply(det, function(i){
	data.table(i@data)[, .N, by = 'Frame']
}))[order(Frame)][, .(N = mean(N)), by = 'Frame'][N == 1 & Frame < 500, N := 0]

det_avg <- rbindlist(lapply(det, function(i){
	setDT(i@data)
	x <- merge(i@data[, .N, by = 'Frame'], data.table(Frame = 1:21600), all = TRUE, by = 'Frame')
	if(i@date != det[[7]]@date){
		x[is.na(N), 'N'] <- 0		
	}
	x
}))[, .(N = mean(N, na.rm = T)), by = 'Frame']
# det_avg[['condition']] <- 'det'
# det_avg[['.id']] <- -1

data_peak <- rbindlist(data, idcol = TRUE)
# data_peak <- rbindlist(data[lapply(k, function(i) i == 'norm') == TRUE], idcol = TRUE)
# data_peak <- rbindlist(data[lapply(k, function(i) i %in% c('norm', 'flat')) == TRUE], idcol = TRUE)

filtered_peak <- rbindlist(data[lapply(k, function(i) i == 'norm') == TRUE], idcol = TRUE)

# data_avg <- data_peak[, .(N = mean(N, na.rm = TRUE)), by = 'Frame']

det_N <- rbindlist(lapply(det, function(i){
	setDT(i@data)
	x <- merge(i@data[, .N, by = 'Frame'], data.table(Frame = 1:21600), all = TRUE, by = 'Frame')
	if(i@date != det[[7]]@date){
		x[is.na(N), 'N'] <- 0		
	}
	x[, .(N = movingAverage(N, t = 60)[1:695], Frame = seq(31, 21545, 31))]
}), idcol = TRUE) 


plot_clusters(data[lapply(k, function(i) i) == 'norm'], k[k == 'norm'])+
	geom_path(data = det_avg, color = 'black', linewidth = 2, alpha = 0.5)+
	geom_path(data = det_avg, color = 'grey50', linewidth = 1.25, alpha = 0.5)



ggplot(data = filtered_peak, aes(Frame, N)) +
	geom_path(aes(group = .id), color = 'grey80', alpha = 0.4, linewidth = 0.8)+
	geom_path(data = det_N, aes(group = .id), color = 'grey30', alpha = 0.5, linewidth = 0.8)+
	
	geom_path(data = filtered_avg, color = 'black', linewidth = 3, alpha = 0.5)+
	geom_path(data = filtered_avg, aes(color = 'Simulations'), linewidth = 2, alpha = 0.5)+
	geom_path(data = det_avg, color = 'black', linewidth = 3, alpha = 0.5)+
	geom_path(data = det_avg, aes(color = 'Experiment'), linewidth = 2, alpha = 0.5)+
	# geom_path(data = rbindlist(list(det_avg, full_data_avg), idcol = TRUE),
	# 	  aes(color = factor(.id, labels = c('Simulations', 'Empirical data'))), linewidth = 2)+
	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50)*120, labels = seq(0, 150, 50))+
	scale_y_continuous('Occupancy')+
	scale_color_manual('',values = c('grey30', 'grey80'))

ggplot(data = data_peak, aes(Frame, N)) +
	geom_path(aes(group = .id), color = 'grey80', alpha = 0.4, linewidth = 0.8)+
	geom_path(data = det_N, aes(group = .id), color = 'grey30', alpha = 0.5, linewidth = 0.8)+

	geom_path(data = full_data_avg, color = 'black', linewidth = 3, alpha = 0.5)+
	geom_path(data = full_data_avg, aes(color = 'Simulations'), linewidth = 2, alpha = 0.5)+
	geom_path(data = det_avg, color = 'black', linewidth = 3, alpha = 0.5)+
	geom_path(data = det_avg, aes(color = 'Experiment'), linewidth = 2, alpha = 0.5)+
	# geom_path(data = rbindlist(list(det_avg, full_data_avg), idcol = TRUE),
	# 	  aes(color = factor(.id, labels = c('Simulations', 'Empirical data'))), linewidth = 2)+
	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50)*120, labels = seq(0, 150, 50))+
	scale_y_continuous('Occupancy')+
	scale_color_manual('',values = c('grey30', 'grey80'))




















path <- '/home/polfer/research/gits/AutomatAnts/results/2024/beta_0.66/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('beta=0.66_\\d{1,2}', files))))

sim_effs <- 1/(rbindlist(lapply(seq_along(files), function(i){
	food <- read_parquet(paste0(path, files[i], '_food.parquet'))[['t']]
	data.table(tmin = min(food), tmax = max(food))
})))

det_effs <- 1/(rbindlist(lapply(det, function(i){
	data.table(tmin = min(rbindlist(i@food)[['t']]), tmax = max(rbindlist(i@food)[['t']]))
	}))/2)

sto_effs <- 1/(rbindlist(lapply(sto, function(i){
	data.table(tmin = min(rbindlist(i@food)[['t']]), tmax = max(rbindlist(i@food)[['t']]))
}))/2)

ggplot(data = melt(rbindlist(list(sims = sim_effs, exps = rbind(det_effs, sto_effs)), idcol = TRUE),
		   id.vars = '.id'), aes(.id, value))+
	geom_boxplot() + 
	facet_wrap(~factor(variable, labels = c('Exploration', 'Exploitation')), scales = 'free')













