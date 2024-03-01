#### LIBRARIES, DATA AND GENERIC FUNCTIONS ####
library(dtw)
library(dtwclust)
library(latex2exp)
library(arrow)
library(infotheo)
source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
load("/home/polfer/research/gits/AutomatAnts/data/zSeries.RData")

load('~/research/gits/AnTracks/data/det.RData')

#### FUNCTIONS ####
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


# example_data <- data.table(read_parquet(paste0(path, files[1], '_data.parquet')))[, c('T', 'id_out', 'pos')]
# processed_data <- example_data[, .(id = parse_ids(id_out), pos = parse_nodes(pos)), by = 'T']
# xy <- setDT(apply(do.call('rbind', strsplit(gsub('[()]', '', processed_data[['pos']]), ', ')), 2, as.numeric))
# colnames(xy) <- c('x', 'y')
# xy[['id']] <- processed_data[['id']]
# xy[['t']] <- processed_data[['T']]
# 
# hist(,
#      breaks = 1000,
#      xlim = c(0, 1000))
# 
# ggplot(data = xy[, .(v = c(0, sqrt(diff(x)^2 + diff(y)^2) / diff(t))), 
# 				by = 'id'],
#        aes(v))
# 
# hist(det[[1]]@data[, .(v = c(0, sqrt(diff(Xmm)^2 + diff(Ymm)^2) / diff(Time_sec))/50), by = 'N_ind'][['v']], 
#      breaks = seq(0, 20, 0.005), xlim = c(0, 1))
# 
# ggplot(data = data.table(v = unlist(lapply(det, function(i){
# 	i@data[, .(v = c(0, sqrt(diff(Xmm)^2 + diff(Ymm)^2) / diff(Time_sec))/50), by = 'N_ind'][['v']]
# }))), aes(v)) + geom_histogram(breaks = seq(0, 1, 0.005), fill = 'grey80', color = 'black')+
# 	scale_x_continuous('Speed (nodes / s)') + ylab('Frequency')

path <- '/home/polfer/research/gits/AutomatAnts/results/2024/Jij_0.01/'
# path <- '/home/polfer/research/gits/AutomatAnts/results/2024/beta_exploration/beta_0.7/'
files <- list.files(path)
# files <- unique(unlist(regmatches(files, gregexpr('beta_0.7_\\d{1,2}', files))))
files <- unique(unlist(regmatches(files, gregexpr('Jij_0.01_\\d{1,2}', files))))


#### SPATIAL DYNAMICS
pos <- data.table(read_parquet(paste0(path, files[1], '_positions.parquet')))
food <- data.table(read_parquet(paste0(path, files[1], '_food.parquet')))[, node := revert_node(node)]
food <- as.data.frame(merge(food, hex_sim, by = 'node', sort = FALSE)[, c('x', 'y')])
food <- list(GP1 = food[1:6, ], GP2 = food[7:12, ])
draw_nodes(pos, 
	   add = draw_hexagons(sim_edges, linewidth = 1.5,
	   		    add = geom_foodpatches(food = food, fill = 'mediumpurple', alpha = 0.8)), 
	   z = pos[['z']], size = 3, show.legend = FALSE) + 
	theme(aspect.ratio = 0.5)+
	scale_fill_viridis(option = 'C')
	# scico::scale_fill_scico()

#### POPULATION DYNAMICS ####

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

full_data <- rbindlist(lapply(files, function(i){
	data.table(read_parquet(paste0(path, i,'.parquet')))
}), idcol = TRUE)[order(Frame)]
full_data_avg <- full_data[, .(N = mean(N)), by = 'Frame']
filtered_avg <- full_data[.id %in% which(k == 'norm'), .(N = mean(N)), by = 'Frame']

foods <- rbindlist(lapply(files, function(i){
	data.table(read_parquet(paste0(path, i,'_food.parquet')))[, c('t')]
}), idcol = TRUE)[, .(mint = min(t), maxt = max(t)), by = '.id'][is.finite(mint),
								 lapply(.SD, mean), .SDcols = c('mint', 'maxt')]
foods_det <- rbindlist(lapply(det, function(i){
	a <- rbindlist(i@food)[['t']]
	data.table(mint = min(a)/2, maxt = max(a)/2)
}))[, lapply(.SD, mean)]

det_avg <- rbindlist(lapply(det, function(i){
	setDT(i@data)
	x <- merge(i@data[, .N, by = 'Frame'], data.table(Frame = 1:21600), all = TRUE, by = 'Frame')
	if(i@date != det[[7]]@date){
		x[is.na(N), 'N'] <- 0		
	}
	x
}))[, .(N = mean(N, na.rm = T)), by = 'Frame']

data_peak <- rbindlist(data, idcol = TRUE)
filtered_peak <- rbindlist(data[lapply(k, function(i) i == 'norm') == TRUE], idcol = TRUE)

det_N <- rbindlist(lapply(det, function(i){
	setDT(i@data)
	x <- merge(i@data[, .N, by = 'Frame'], data.table(Frame = 1:21600), all = TRUE, by = 'Frame')
	if(i@date != det[[7]]@date){
		x[is.na(N), 'N'] <- 0		
	}
	x[, .(N = movingAverage(N, t = 60)[1:695], Frame = seq(31, 21545, 31))]
}), idcol = TRUE) 

# filtered by "exp-pattern"
# ggplot(data = filtered_peak, aes(Frame, N)) +
# 	geom_path(aes(group = .id), color = 'grey80', alpha = 0.4, linewidth = 0.8)+
# 	geom_path(data = det_N, aes(group = .id), color = 'grey30', alpha = 0.5, linewidth = 0.8)+
# 	
# 	geom_path(data = filtered_avg, color = 'black', linewidth = 3, alpha = 0.5)+
# 	geom_path(data = filtered_avg, aes(color = 'Simulations'), linewidth = 2, alpha = 0.5)+
# 	geom_path(data = det_avg, color = 'black', linewidth = 3, alpha = 0.5)+
# 	geom_path(data = det_avg, aes(color = 'Experiment'), linewidth = 2, alpha = 0.5)+
# 	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50)*120, labels = seq(0, 150, 50))+
# 	scale_y_continuous('Occupancy')+
# 	scale_color_manual('',values = c('grey30', 'grey80'))+
# 	geom_vline(xintercept = unlist(foods)*2, linetype = 2, linewidth = 1)+
# 	geom_vline(xintercept = unlist(foods_det)*2, linetype = 3, linewidth = 1)

# all simulations
ggplot(data = data_peak, aes(Frame, N)) +
	geom_path(aes(group = .id), color = 'grey80', alpha = 0.4, linewidth = 0.8)+
	geom_path(data = det_N, aes(group = .id), color = 'grey30', alpha = 0.5, linewidth = 0.8)+
	
	geom_path(data = full_data_avg, color = 'black', linewidth = 3, alpha = 0.5)+
	geom_path(data = full_data_avg, aes(color = 'Simulations'), linewidth = 2, alpha = 0.5)+
	geom_path(data = det_avg, color = 'black', linewidth = 3, alpha = 0.5)+
	geom_path(data = det_avg, aes(color = 'Experiment'), linewidth = 2, alpha = 0.5)+
	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50)*120, labels = seq(0, 150, 50))+
	scale_y_continuous('Occupancy')+
	scale_color_manual('',values = c('grey30', 'grey80'))+
	geom_vline(xintercept = unlist(foods)*2, linetype = 2, linewidth = 1)+
	geom_vline(xintercept = unlist(foods_det)*2, linetype = 3, linewidth = 1)


#### FORAGING EFFICIENCY ####
global_sim_effs <- global_eff(path)[is.finite(maxt)][, condition := 'Simulations']
det_global_eff <- rbindlist(lapply(det, function(i){
	x <- rbindlist(i@food)[['t']]
	1/(data.table(mint = min(x), maxt = max(x))/2)
}))[, condition := 'Experiments']

# global eff
ggplot(data = melt(rbind(global_sim_effs[, -'.id'], det_global_eff), id.vars = 'condition'),
       aes(condition, value, fill = condition)) + geom_boxplot(alpha = 0.6) +
	facet_wrap(~ variable, scales = 'free')+
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')),
		   scales = 'free_y')+
	scale_x_discrete('') + 
	ylab(TeX('$Efficiency (s^{-1})$'))+
	scale_fill_manual('', values = c('mediumpurple', 'gold3'))+
	theme(legend.key.size = unit(30, 'pt'))


# collective efficiency
collective_sim_effs <- lapply(seq_along(files), function(i){
	do.call('gc', args = list(verbose = FALSE))
	eff(path, files[i])
})

det_tps <- lapply(det, get_eff)
# sto_tps <- lapply(sto, get_eff)
# tps <- cbind(condition = 'exp', rbindlist(l = list(rbindlist(det_tps),rbindlist(sto_tps))))

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

