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
	
# ids[[1]][['data']][ids[[1]][['stops']]+1][, d := ids[[1]][['stop_duration']]]
# 
# ind_5 <- det[[1]]@data[N_ind == 5]
# ind_6 <- det[[1]]@data[N_ind == 6]
# 
# Dt <- diff(ind_5[, Frame])
# stops <- which(Dt != 1)
# stop_duration <- Dt[Dt != 1]
# 
# draw_hexagons(add = geom_foodpatches()) +
# 	geom_point(data = ind_5, aes(Xmm, Ymm, size = 1/Frame, color = Frame), show.legend = FALSE)+
# 	scale_color_viridis_c()+
# 	scale_size_continuous(range = c(3, 9))+
# 	geom_point(data = ind_5[stops +1][, d := stop_duration], aes(Xmm, Ymm, size = d), shape = 13, 
# 		   show.legend = FALSE)
# 
# 
# 
# 
# det[[1]]@data[N_ind == 5 & Crossings == 1]
# det[[1]]@data[Frame == 1601]
# geom_foodpatches(add = draw_hexagons(add = ggplot(data = ind_5, aes(Xmm, Ymm))+
# 				     	geom_point(aes(size = 1/Frame, color = Frame), show.legend = FALSE)+
# 				     	scale_color_viridis_c()))
# geom_foodpatches(add = draw_hexagons(add = ggplot(data = ind_6, aes(Xmm, Ymm))+
# 				     	geom_point(aes(size = 1/Frame, color = Frame), show.legend = FALSE)+
# 				     	scale_color_viridis_c()))
