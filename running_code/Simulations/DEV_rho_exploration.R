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


det_global_eff <- rbindlist(lapply(det, function(i){
	x <- rbindlist(i@food)[['t']]
	mint <- min(x) - min(i@data[['Frame']])
	1/(data.table(mint = mint, maxt = max(x)-min(x))/2)
}))[, condition := 'Experiments']


####### RHO = 0.5 AND RECRUITMENT #######

path_rho <- '/home/polfer/research/gits/AutomatAnts/results/2024/rho/wRec/rho_0.5/'
f_rho <- list.files(path_rho)
files_rho <- f_rho[grepl('data', f_rho)]

eff_rho_0.5 <- global_eff(path_rho)[, condition := 'Simulations']

#### efficiency rho = 0.5 ####
ggplot(data = melt(rbind(eff_rho_0.5[, -'.id'], det_global_eff), id.vars = 'condition'),
       aes(factor(condition, labels = c('Experiments', 'rho')), value, fill = condition)) +
	geom_boxplot(alpha = 0.6, show.legend = FALSE, outlier.shape = NA) +
	geom_jitter(width = 0.15, alpha = 0.25, size = 3, show.legend = FALSE)+
	facet_wrap(~ variable, scales = 'free')+
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')),
		   scales = 'free_y')+
	scale_x_discrete('') + 
	ylab(TeX('$Efficiency (s^{-1})$'))+
	scale_fill_manual('', values = c('mediumpurple', 'gold3'))+
	theme(legend.key.size = unit(30, 'pt'))


f_rho <- list.files(path_rho)
files_rho <- f_rho[!grepl('data', f_rho) & !grepl('food', f_rho) & !grepl('keys', f_rho) &
		   	!grepl('positions', f_rho)]

data <- process_data(path_rho)

full_data <- rbindlist(lapply(files_rho, function(i){
	data.table(read_parquet(paste0(path_rho, i)))
}), idcol = TRUE)[order(Frame)]
full_data_avg <- full_data[, .(N = mean(N)), by = 'Frame']


foods <- rbindlist(lapply(f_rho[grepl('food', f_rho)], function(i){
	data.table(read_parquet(paste0(path_rho, i)))[, c('t')]
}), idcol = TRUE)[, .(mint = min(t), maxt = max(t)), by = '.id'][is.finite(mint),
								 lapply(.SD, median), .SDcols = c('mint', 'maxt')]


# collective efficiency
coll_eff_rho_0.5 <- lapply(gsub('_data.parquet','',files_rho), function(i){
	do.call('gc', args = list(verbose = FALSE))
	eff(path_rho, i)
})

det_tps <- lapply(det, get_eff)

coll_rho_0.5 <- rbindlist(coll_eff_rho_0.5)
coll_rho_0.5[['condition']] <- 'Simulations'

ggplot(data = rbind(melt(rbind(eff_rho_0.5[, -'.id'], det_global_eff), id.vars = 'condition'),
			 data.table(melt(rbindlist(list(coll_rho_0.5, 
			 			       cbind(condition = 'Experiments',
			 			             rbindlist(det_tps))), use.names = TRUE), 
			 		id.vars = 'condition'))[, value := 1/value]), 
			 aes(condition, value, fill = condition))+
	
	geom_boxplot(alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +
	geom_jitter(alpha = 0.25, width = 0.15, size = 3, show.legend = FALSE)+
	facet_wrap(~ factor(variable, labels = c('Exploration time', 'Exploitation time', 
						 'Exploration efficiency', 'Exploitation efficiency')),
		   scales = 'free_y')+
	scale_x_discrete('') + 
	ylab(TeX('$Efficiency (s^{-1})$'))+
	scale_fill_manual('', values = c('mediumpurple', 'gold3'))+
	theme(legend.key.size = unit(30, 'pt'))


data_peak <- rbindlist(data, idcol = TRUE)
#### N rho = 0.5 ####
ggplot(data = data_peak, aes(Frame, N)) +
	geom_path(aes(group = .id), color = 'grey80', alpha = 0.2, linewidth = 0.8)+
	geom_path(data = det_N, aes(group = .id), color = 'grey30', alpha = 0.5, linewidth = 0.8)+
	
	geom_path(data = full_data_avg, color = 'black', linewidth = 3, alpha = 0.5)+
	geom_path(data = full_data_avg, aes(color = 'Simulations'), linewidth = 2, alpha = 0.5)+
	geom_path(data = det_avg, color = 'black', linewidth = 3, alpha = 0.5)+
	geom_path(data = det_avg, aes(color = 'Experiments'), linewidth = 2, alpha = 0.5)+
	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50)*120, labels = seq(0, 150, 50))+
	scale_y_continuous('Occupancy')+
	scale_color_manual('',values = c('grey30', 'grey80'))+
	geom_vline(xintercept = unlist(foods)*2, linetype = 2, linewidth = 1)+
	geom_vline(xintercept = unlist(foods_det)*2, linetype = 3, linewidth = 1)+
	theme(legend.position = c(0.8, 0.85),
	      legend.background = element_rect(fill = 'white', color = 'black'),
	      legend.title = element_blank())

#### positions rho = 0.5 ####
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

food <- data.table(read_parquet(paste0(path_rho, f_rho[grepl('food', f_rho)][1])))[, node := revert_node(node)]
food <- as.data.frame(merge(food, hex_sim, by = 'node', sort = FALSE)[, c('x', 'y')])
food <- list(GP1 = food[1:6, ], GP2 = food[7:12, ])

pos_det <- merge(hex, rbindlist(lapply(det, function(i){
	i@data[, .N, by = 'node']
}), idcol = TRUE)[, .(N = sum(N)), by = 'node'], by = 'node')
files_rho <- f_rho[grepl('data', f_rho)]

pos_data <- rbindlist(lapply(files_rho, function(i){ ## <--- takes a couple of minutes and a bit of RAM
	do.call('gc', args = list(verbose = FALSE))
	data.table(read_parquet(paste0(path_rho, i)))[, .(node = parse_nodes(pos))][, .N, by = 'node']
}))

sum_pos_data <- merge(hex_sim, pos_data[, .(N = sum(N)), by = 'node'], by = 'node')

draw_nodes(sum_pos_data, 
	   add = draw_hexagons(sim_edges, linewidth = 2, color = 'black',
	   		    add = draw_hexagons(sim_edges, linewidth = 2, color = 'black', 
	   		    		    add = geom_foodpatches(food = food, fill = 'mediumpurple', alpha = 0.9))), 
	   z = rank(sum_pos_data[['N']]), size = 4, show.legend = FALSE) + 
	theme(aspect.ratio = 0.5)+
	scale_fill_viridis(option = 'plasma')

draw_nodes(pos_det, 
	   add = draw_hexagons(edges, linewidth = 2, color = 'black',
	   		    add = draw_hexagons(edges, linewidth = 2, color = 'black', 
	   		    		    add = geom_foodpatches(fill = 'mediumpurple', alpha = 0.9))), 
	   z = rank(pos_det[['N']]), size = 4) + 
	scale_fill_viridis('Occupancy',option = 'plasma', breaks = c(1, 620),
			   labels = c('Low', 'High')) +
	theme(legend.position = 'bottom', aspect.ratio = 0.5,
	      legend.margin = margin(t = -20, unit = 'pt'))+
	guides(fill = guide_colorbar(title.position = 'top', barwidth = 15, title.hjust = 0.5))
####### RHO EXPLORATION #######
path_rho <- '/home/polfer/research/gits/AutomatAnts/results/2024/rho/wRec/rho_0.9/'
f_rho <- list.files(path_rho)
files_rho <- f_rho[grepl('data', f_rho)]

eff_rho_0.1 <- global_eff(path_rho)[, condition := 'rho_0.1']
eff_rho_0.5 <- global_eff(path_rho)[, condition := 'rho_0.5']
eff_rho_0.9 <- global_eff(path_rho)[, condition := 'rho_0.9']

ggplot(data = melt(rbind(eff_rho_0.9[, -'.id'], det_global_eff), id.vars = 'condition'),
       aes(factor(condition, labels = c('Experiments', 'rho')), value, fill = condition)) +
	geom_boxplot(alpha = 0.6, show.legend = FALSE, outlier.shape = NA) +
	geom_jitter(width = 0.15, alpha = 0.25, size = 3, show.legend = FALSE)+
	facet_wrap(~ variable, scales = 'free')+
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')),
		   scales = 'free_y')+
	scale_x_discrete('', labels = c('Experiments', TeX('$\\rho = 0.9$'))) + 
	ylab(TeX('$Efficiency (s^{-1})$'))+
	scale_fill_manual('', values = c('mediumpurple', 'gold3', 'brown4'))+
	theme(legend.key.size = unit(30, 'pt'))




path_rho <- '/home/polfer/research/gits/AutomatAnts/results/2024/rho/wRec/rho_0.9/'
f_rho <- list.files(path_rho)
files_rho <- f_rho[!grepl('data', f_rho) & !grepl('food', f_rho) & !grepl('keys', f_rho) &
		   	!grepl('positions', f_rho)]

data <- process_data(path_rho)

full_data <- rbindlist(lapply(files_rho, function(i){
	data.table(read_parquet(paste0(path_rho, i)))
}), idcol = TRUE)[order(Frame)]
full_data_avg <- full_data[, .(N = mean(N)), by = 'Frame']


foods <- rbindlist(lapply(f_rho[grepl('food', f_rho)], function(i){
	data.table(read_parquet(paste0(path_rho, i)))[, c('t')]
}), idcol = TRUE)[, .(mint = min(t), maxt = max(t)), by = '.id'][is.finite(mint),
								 lapply(.SD, median), .SDcols = c('mint', 'maxt')]


data_peak <- rbindlist(data, idcol = TRUE)



ggplot(data = data_peak, aes(Frame, N)) +
	geom_path(aes(group = .id), color = 'grey80', alpha = 0.2, linewidth = 0.8)+
	geom_path(data = det_N, aes(group = .id), color = 'grey30', alpha = 0.5, linewidth = 0.8)+
	
	geom_path(data = full_data_avg, color = 'black', linewidth = 3, alpha = 0.5)+
	geom_path(data = full_data_avg, aes(color = 'Simulations'), linewidth = 2, alpha = 0.5)+
	geom_path(data = det_avg, color = 'black', linewidth = 3, alpha = 0.5)+
	geom_path(data = det_avg, aes(color = 'Experiments'), linewidth = 2, alpha = 0.5)+
	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50)*120, labels = seq(0, 150, 50))+
	scale_y_continuous('Occupancy')+
	scale_color_manual('',values = c('grey30', 'grey80'))+
	geom_vline(xintercept = unlist(foods)*2, linetype = 2, linewidth = 1)+
	geom_vline(xintercept = unlist(foods_det)*2, linetype = 3, linewidth = 1)+
	theme(legend.position = c(0.8, 0.85),
	      legend.background = element_rect(fill = 'white', color = 'black'),
	      legend.title = element_blank())

#### WITHOUT RECRUITMENT ####

foods_det <- rbindlist(lapply(det, function(i){
	a <- rbindlist(i@food)[['t']]
	data.table(mint = min(a)/2, maxt = max(a)/2)
}))[, lapply(.SD, median)] # <--- median !

det_avg <- rbindlist(lapply(det, function(i){
	setDT(i@data)
	x <- merge(i@data[, .N, by = 'Frame'], data.table(Frame = 1:21600), all = TRUE, by = 'Frame')
	if(i@date != det[[7]]@date){
		x[is.na(N), 'N'] <- 0		
	}
	x
}))[, .(N = mean(N, na.rm = T)), by = 'Frame']

det_N <- rbindlist(lapply(det, function(i){
	setDT(i@data)
	x <- merge(i@data[, .N, by = 'Frame'], data.table(Frame = 1:21600), all = TRUE, by = 'Frame')
	if(i@date != det[[7]]@date){
		x[is.na(N), 'N'] <- 0		
	}
	x[, .(N = movingAverage(N, t = 60)[1:695], Frame = seq(31, 21545, 31))]
}), idcol = TRUE)

path_rho <- '/home/polfer/research/gits/AutomatAnts/results/2024/rho/woutRec/rho_0.1/'
f_rho <- list.files(path_rho)
files_rho <- f_rho[grepl('data', f_rho)]

global_eff <- function(path){
	files <- list.files(path)[grepl('food', list.files(path))]
	results <- lapply(files, function(i){
		f <- gsub('_food', '', i)
		d <- data.table(read_parquet(paste0(path, f)))
		m1 <- min(d[N >0, Frame])/2
		food <- data.table(read_parquet(paste0(path, i)))
		mint <- min(food[['t']]) - m1
		1/data.table(mint = mint, maxt = max(food[['t']])-min(food[['t']]))
	})
	rbindlist(results, idcol = TRUE)
}

eff_rho_0.1 <- global_eff(path_rho)[, condition := 'rho_0.1']
eff_rho_0.5 <- global_eff(path_rho)[, condition := 'rho_0.5']
eff_rho_0.9 <- global_eff(path_rho)[, condition := 'rho_0.9']

det_global_eff <- rbindlist(lapply(det, function(i){
	x <- rbindlist(i@food)[['t']]
	mint <- min(x) - min(i@data[['Frame']])
	1/(data.table(mint = mint, maxt = max(x)-min(x))/2)
}))[, condition := 'Experiments']

ggplot(data = melt(rbind(eff_rho_0.1[, -'.id'], det_global_eff), id.vars = 'condition'),
       aes(factor(condition, labels = c('Experiments', 'rho')), value, fill = condition)) +
	geom_boxplot(alpha = 0.6, show.legend = FALSE, outlier.shape = NA) +
	geom_jitter(width = 0.15, alpha = 0.25, size = 3, show.legend = FALSE)+
	facet_wrap(~ variable, scales = 'free')+
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')),
		   scales = 'free_y')+
	scale_x_discrete('', labels = c('Experiments', TeX('$\\rho = 0.9$'))) + 
	ylab(TeX('$Efficiency (s^{-1})$'))+
	scale_fill_manual('', values = c('mediumpurple', 'gold3', 'brown4'))+
	theme(legend.key.size = unit(30, 'pt'))




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

library(dtwclust)



path_rho <- '/home/polfer/research/gits/AutomatAnts/results/2024/rho/woutRec/rho_0.5/'
f_rho <- list.files(path_rho)
files_rho <- f_rho[!grepl('data', f_rho) & !grepl('food', f_rho) & !grepl('keys', f_rho) &
		   	!grepl('positions', f_rho)]

data <- process_data(path_rho)

full_data <- rbindlist(lapply(files_rho, function(i){
	data.table(read_parquet(paste0(path_rho, i)))
}), idcol = TRUE)[order(Frame)]
full_data_avg <- full_data[, .(N = mean(N)), by = 'Frame']


foods <- rbindlist(lapply(f_rho[grepl('food', f_rho)], function(i){
	data.table(read_parquet(paste0(path_rho, i)))[, c('t')]
}), idcol = TRUE)[, .(mint = min(t), maxt = max(t)), by = '.id'][is.finite(mint),
								 lapply(.SD, median), .SDcols = c('mint', 'maxt')]


data_peak <- rbindlist(data, idcol = TRUE)



ggplot(data = data_peak, aes(Frame, N)) +
	geom_path(aes(group = .id), color = 'grey80', alpha = 0.2, linewidth = 0.8)+
	geom_path(data = det_N, aes(group = .id), color = 'grey30', alpha = 0.5, linewidth = 0.8)+
	
	geom_path(data = full_data_avg, color = 'black', linewidth = 3, alpha = 0.5)+
	geom_path(data = full_data_avg, aes(color = 'Simulations'), linewidth = 2, alpha = 0.5)+
	geom_path(data = det_avg, color = 'black', linewidth = 3, alpha = 0.5)+
	geom_path(data = det_avg, aes(color = 'Experiments'), linewidth = 2, alpha = 0.5)+
	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50)*120, labels = seq(0, 150, 50))+
	scale_y_continuous('Occupancy')+
	scale_color_manual('',values = c('grey30', 'grey80'))+
	geom_vline(xintercept = unlist(foods)*2, linetype = 2, linewidth = 1)+
	geom_vline(xintercept = unlist(foods_det)*2, linetype = 3, linewidth = 1)+
	theme(legend.position = c(0.8, 0.85),
	      legend.background = element_rect(fill = 'white', color = 'black'),
	      legend.title = element_blank())
