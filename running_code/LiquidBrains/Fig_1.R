#### LIBRARIES, DATA AND GENERIC FUNCTIONS ####
library(dtw)
library(dtwclust)
library(latex2exp)
library(arrow)
# library(infotheo)
source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')

load('~/research/gits/AnTracks/data/det.RData')


## EXPERIMENTS
# foods_det <- rbindlist(lapply(det, function(i){
# 	a <- rbindlist(i@food)[['t']]
# 	data.table(mint = min(a)/2, maxt = max(a)/2)
# }))[, lapply(.SD, mean)] # <--- mean !
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

path <- '/home/polfer/research/gits/AutomatAnts/results/2024/default/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('default_\\d{1,2}', files))))


#### SPATIAL DYNAMICS
pos <- data.table(read_parquet(paste0(path, files[3], '_positions.parquet')))
food <- data.table(read_parquet(paste0(path, files[1], '_food.parquet')))[, node := revert_node(node)]
food <- as.data.frame(merge(food, hex_sim, by = 'node', sort = FALSE)[, c('x', 'y')])
food <- list(GP1 = food[1:6, ], GP2 = food[7:12, ])
c <- draw_nodes(pos, 
	   add = draw_hexagons(sim_edges, linewidth = 2, color = 'black',
	   		    add = draw_hexagons(sim_edges, linewidth = 2, color = 'black', 
	   		    		    add = geom_foodpatches(food = food, fill = 'mediumpurple', alpha = 0.9))), 
	   z = pos[['z']], size = 4, show.legend = FALSE) + 
	theme(aspect.ratio = 0.5)+
	scale_fill_viridis(option = 'plasma')

pos_det <- merge(hex, rbindlist(lapply(det, function(i){
	i@data[, .N, by = 'node']
}), idcol = TRUE)[, .(N = sum(N)), by = 'node'], by = 'node')

c2 <- draw_nodes(c1, 
	   add = draw_hexagons(edges, linewidth = 2, color = 'black',
	   		    add = draw_hexagons(edges, linewidth = 2, color = 'black', 
	   		    		    add = geom_foodpatches(fill = 'mediumpurple', alpha = 0.9))), 
	   z = rank(c1[['N']]), size = 4) + 
	scale_fill_viridis('Occupancy',option = 'plasma', breaks = c(1, 620),
			   labels = c('Low', 'High')) +
	theme(legend.position = 'bottom', aspect.ratio = 0.5,
	      legend.margin = margin(t = -20, unit = 'pt'))+
	guides(fill = guide_colorbar(title.position = 'top', barwidth = 15, title.hjust = 0.5))
	
	
# scico::scale_fill_scico()

#### POPULATION DYNAMICS ####

data <- process_data(path)

full_data <- rbindlist(lapply(files, function(i){
	data.table(read_parquet(paste0(path, i,'.parquet')))
}), idcol = TRUE)[order(Frame)]
full_data_avg <- full_data[, .(N = mean(N)), by = 'Frame']

foods <- rbindlist(lapply(files, function(i){
	data.table(read_parquet(paste0(path, i,'_food.parquet')))[, c('t')]
}), idcol = TRUE)[, .(mint = min(t), maxt = max(t)), by = '.id'][is.finite(mint),
								 lapply(.SD, mean), .SDcols = c('mint', 'maxt')]
foods <- rbindlist(lapply(files, function(i){
	data.table(read_parquet(paste0(path, i,'_food.parquet')))[, c('t')]
}), idcol = TRUE)[, .(mint = min(t), maxt = max(t)), by = '.id'][is.finite(mint),
								 lapply(.SD, median), .SDcols = c('mint', 'maxt')]


data_peak <- rbindlist(data, idcol = TRUE)
 
# comparison of population dynamics
a <- ggplot(data = data_peak, aes(Frame, N)) +
	geom_path(aes(group = .id), color = 'grey80', alpha = 0.2, linewidth = 0.8)+
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

a <- ggplot(data = data_peak, aes(Frame, N)) +
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
	theme(legend.position = c(0.8, 0.75),
	      legend.background = element_rect(fill = 'white', color = 'black'),
	      legend.title = element_blank())


#### FORAGING EFFICIENCY ####
global_sim_effs <- global_eff(path)[is.finite(maxt)][, condition := 'Simulations']
det_global_eff <- rbindlist(lapply(det, function(i){
	x <- rbindlist(i@food)[['t']]
	# print(min(i@data[['Frame']]))
	mint <- min(x) - min(i@data[['Frame']])
	# 1/(data.table(mint = mint, maxt = max(x)-min(x), t0 = 1/min(i@data[['Frame']]))/2)
	1/(data.table(mint = mint, maxt = max(x)-min(x))/2)
	# 1/(data.table(mint = min(x), maxt = max(x))/2)
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

coll_dt <- rbindlist(collective_sim_effs)
coll_dt[['condition']] <- 'Simulations'

# ggplot(data = melt(rbindlist(list(coll_dt, cbind(condition = 'Experiments', rbindlist(det_tps))), use.names = TRUE), id.vars = 'condition'), 
#        aes(condition, 1/value, fill = condition))+
# 	geom_boxplot(alpha = 0.6)+ 
# 	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')),
# 		   scales = 'free_y')+
# 	scale_x_discrete('') + 
# 	ylab(TeX('$Efficiency (s^{-1})$'))+
# 	scale_fill_manual('', values = c('mediumpurple', 'gold3'))+
# 	theme(legend.key.size = unit(30, 'pt'))

b <- ggplot(data = rbind(melt(rbind(global_sim_effs[, -'.id'], det_global_eff), id.vars = 'condition'),
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

b <- ggplot(data = rbind(melt(rbind(global_sim_effs[, -'.id'], det_global_eff), id.vars = 'condition'),
			 data.table(melt(rbindlist(list(coll_dt, 
			 			       cbind(condition = 'Experiments',
			 			             rbindlist(det_tps))), use.names = TRUE), 
			 		id.vars = 'condition'))[, value := 1/value]), 
			 aes(condition, value, fill = condition))+
	geom_boxplot(alpha = 0.6, show.legend = FALSE) +
	facet_wrap(~ factor(variable, labels = c('Exploration time', 'Exploitation time', 
						 'Exploration efficiency', 'Exploitation efficiency')),
		   scales = 'free_y')+
	scale_x_discrete('') + 
	ylab(TeX('$Efficiency (s^{-1})$'))+
	scale_fill_manual('', values = c('mediumpurple', 'gold3'))+
	theme(legend.key.size = unit(30, 'pt'))


# grid.arrange(grobs = list(a, b, c),
# 	     layout_matrix = rbind(c(rep(1, 10), rep(2, 10)), 
# 	     		      c(rep(3, 8),NA, NA, rep(2, 10))))
grid.arrange(grobs = list(a, b, c),
	     layout_matrix = rbind(c(1, 2), c(3, 2)))

# grid.arrange(grobs = list(a + theme(aspect.ratio = 0.5), b, c2, c),
# 	     layout_matrix = rbind(c(1, 2), c(1, 2), c(3, 2), 
# 	     		      c(3, NA), c(4, NA), c(4, NA)))


png('/home/polfer/research/gits/AnTracks/plots/summary_default.png', 8000, 6000, res = 550)
grid.arrange(grobs = list(a + theme(aspect.ratio = 0.5), b, 
			  ggarrange(c2, c, ncol = 1, common.legend = TRUE, legend = 'bottom')),
			  layout_matrix = rbind(c(1, 2), c(1, 2), c(3, 2), 
			  		      c(3, NA), c(3, NA), c(3, NA)))
dev.off()


########################################################
############ FALTA MUTUAL INFORMATION !!!!! ############
########################################################



