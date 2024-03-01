source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
library(arrow)
library(latex2exp)

path <- '/home/polfer/research/gits/AutomatAnts/results/2024/uniform_rec/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('uniform_\\d{1,2}', files))))

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

afunc <- function(path, filename, filter = FALSE){
	do.call('gc', args = list(verbose = FALSE))
	food_sim <- data.table(read_parquet(paste0(path, filename,'_food.parquet')))
	food_sim <- food_sim[, .(node = rectify_nodes(node), t = t)]
	dt <- data.table(read_parquet(paste0(path, filename,'_data.parquet')))[, c('Frame', 'pos', 'id_out')]
	if(filter){
		dt <- dt[Frame <= max(food_sim[['t']] * 2)]
	}
	dt_complete <- dt[, .(node = parse_nodes(pos), id = parse_ids(id_out)), by = 'Frame'][!is.na(node)]
	
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



path2 <- '/home/polfer/research/gits/AutomatAnts/results/2024/uniform_noRec/'
files2 <- list.files(path2)
files2 <- unique(unlist(regmatches(files2, gregexpr('uniform_\\d{1,2}', files2))))
t0 <- Sys.time()
times_simulation2 <- c()
for(i in files2){
	times_simulation2 <- c(times_simulation2, afunc(path2, i))
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
tracklets_to_food <- function(exp, filter = TRUE){
	nest_nodes <- c(662, 608, 634, 635)
	fnds <- get_foodnodes(exp)
	dt <- setDT(exp@data)
	ids <- get_ids(exp)
	sbst_exp <- dt[N_ind %in% ids]

	if(filter){
		sbst_exp <- sbst_exp[Frame <= max(rbindlist(exp@food)[['t']])]
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
load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/sto.RData')


tdet_noFil <- unlist(lapply(det, tracklets_to_food, filter = FALSE))
tsto_noFil <- unlist(lapply(sto, tracklets_to_food, filter = FALSE))



ggplot(data = data.frame(x = c(times_simulation[is.finite(times_simulation)] / 120,
			       times_simulation2[is.finite(times_simulation2)] / 120),
			 type = c(rep('With recruitment', sum(is.finite(times_simulation))),
			 	 rep('Without recruitment', sum(is.finite(times_simulation2))))), 
       aes(x, fill = type)) +
	geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 35, 0.15),
		       color = 'black', alpha = 0.5, position = 'identity') + 
	geom_density(alpha = 0.5)+
	scale_x_continuous('Time (min)', breaks = seq(0, 10, 2.5),
			   limits = c(0, 10))+
	# geom_vline(xintercept = c(median(times_simulation[is.finite(times_simulation)] / 120),
	# 			  median( c(tdet_noFil, tsto_noFil)/120)),
	# 	   linewidth = 1.2, linetype = 2)+
	geom_vline(xintercept = c(median(times_simulation2[is.finite(times_simulation2)] / 120),
				  median( times_simulation[is.finite(times_simulation)] / 120)),
		   linewidth = 1, linetype = 2, color = c('mediumpurple', 'gold3'))+
	scale_y_continuous('Density', breaks = seq(0, 10, 0.15)) +
	scale_fill_manual('',values = c('gold3', 'mediumpurple'))+
	annotation_custom(
		ggplotGrob(ggplot(data = data.frame(x = c(times_simulation[is.finite(times_simulation)] / 120,
							  times_simulation2[is.finite(times_simulation2)] / 120),
							  type = c(rep('With recruitment',
							  	     sum(is.finite(times_simulation))),
							  	 rep('Without recruitment', 
							  	     sum(is.finite(times_simulation2))))))+

			   	geom_violin(aes(type, x, fill = type), trim = FALSE,
			   		    alpha = 0.5, show.legend = FALSE)+
			   	theme(
			   		axis.ticks.x = element_blank())+xlab('')+
			   	scale_y_continuous('Time (min)', breaks = seq(0, 100, 10))+
			   	scale_fill_manual(values = c('gold3', 'mediumpurple'))),
		xmin = 4, xmax = 7.5, ymin = 0.2, ymax = 0.65)+
	geom_rect(aes(xmin = 4, xmax = 7.5, ymin = 0.2, ymax = 0.69), fill = NA, color = 'black',
		  linewidth = 0.5)+
	theme(legend.position = c(0.85, 0.5), 
	      legend.background = element_rect(fill = NA, colour = 'black'),
	      legend.title = element_blank())


ggplot(data = data.frame(x = c(times_simulation[is.finite(times_simulation)] / 120,
			       c(tdet_noFil, tsto_noFil)/120),
			 type = c(rep('Simulation', sum(is.finite(times_simulation))),
			 	 rep('Experiment', length(c(tdet_noFil, tsto_noFil))))), 
       aes(x, fill = type)) +
	geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 35, 0.15),
		       color = 'black', alpha = 0.5, position = 'identity') + 
	geom_density(alpha = 0.5)+
	scale_x_continuous('Time (min)', breaks = seq(0, 10, 2.5),
			   limits = c(0, 10))+
	# geom_vline(xintercept = c(median(times_simulation[is.finite(times_simulation)] / 120),
	# 			  median( c(tdet_noFil, tsto_noFil)/120)),
	# 	   linewidth = 1.2, linetype = 2)+
	geom_vline(xintercept = c(median(times_simulation[is.finite(times_simulation)] / 120),
				  median( c(tdet_noFil, tsto_noFil)/120)),
		   linewidth = 1, linetype = 2, color = c('mediumpurple', 'gold3'))+
	scale_y_continuous('Density', breaks = seq(0, 10, 0.15)) +
	scale_fill_manual('',values = c('gold3', 'mediumpurple'))+
	annotation_custom(
		ggplotGrob(ggplot(data = data.frame(x = c(times_simulation[is.finite(times_simulation)] / 120,
							  c(tdet_noFil, tsto_noFil)/120),
							  type = c(rep('Simulation', sum(is.finite(times_simulation))),
							  	 rep('Experiment', length(c(tdet_noFil, tsto_noFil))))))+
			   	# geom_boxplot(aes(type, x, fill = type), width = 0.5,
			   	# 	     alpha = 0.5,outlier.shape = 1, show.legend = FALSE)+
			   	geom_violin(aes(type, x, fill = type), trim = FALSE,
			   		    alpha = 0.5, show.legend = FALSE)+
			   	theme(
			   		axis.ticks.x = element_blank())+xlab('')+
			   	scale_y_continuous('Time (min)', breaks = seq(0, 100, 10))+
			   	scale_fill_manual(values = c('gold3', 'mediumpurple'))),
		xmin = 4, xmax = 7.5, ymin = 0.2, ymax = 0.65)+
	geom_rect(aes(xmin = 4, xmax = 7.5, ymin = 0.2, ymax = 0.69), fill = NA, color = 'black',
		  linewidth = 0.5)+
	theme(legend.position = c(0.85, 0.5), 
	      legend.background = element_rect(fill = NA, colour = 'black'),
	      legend.title = element_blank())


ggplot(data = data.frame(x = c(times_simulation2[is.finite(times_simulation2)] / 120,
			       c(tdet_noFil, tsto_noFil)/120),
			 type = c(rep('Simulation', sum(is.finite(times_simulation2))),
			 	 rep('Experiment', length(c(tdet_noFil, tsto_noFil))))), 
       aes(x, fill = type)) +
	geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 35, 0.15),
		       color = 'black', alpha = 0.5, position = 'identity') + 
	geom_density(alpha = 0.5)+
	scale_x_continuous('Time (min)', breaks = seq(0, 10, 2.5),
			   limits = c(0, 10))+
	# geom_vline(xintercept = c(median(times_simulation2[is.finite(times_simulation2)] / 120),
	# 			  median( c(tdet_noFil, tsto_noFil)/120)),
	# 	   linewidth = 1.2, linetype = 2)+
	geom_vline(xintercept = c(median(times_simulation2[is.finite(times_simulation2)] / 120),
				  median( c(tdet_noFil, tsto_noFil)/120)),
		   linewidth = 1, linetype = 2, color = c('mediumpurple', 'gold3'))+
	scale_y_continuous('Density', breaks = seq(0, 10, 0.15)) +
	scale_fill_manual('',values = c('gold3', 'mediumpurple'))+
	annotation_custom(
		ggplotGrob(ggplot(data = data.frame(x = c(times_simulation2[is.finite(times_simulation2)] / 120,
							  c(tdet_noFil, tsto_noFil)/120),
							  type = c(rep('Simulation', sum(is.finite(times_simulation2))),
							  	 rep('Experiment', length(c(tdet_noFil, tsto_noFil))))))+
			   	# geom_boxplot(aes(type, x, fill = type), width = 0.5,
			   	# 	     alpha = 0.5,outlier.shape = 1, show.legend = FALSE)+
			   	geom_violin(aes(type, x, fill = type), trim = FALSE,
			   		    alpha = 0.5, show.legend = FALSE)+
			   	theme(
			   		axis.ticks.x = element_blank())+xlab('')+
			   	scale_y_continuous('Time (min)', breaks = seq(0, 100, 10))+
			   	scale_fill_manual(values = c('gold3', 'mediumpurple'))),
		xmin = 4, xmax = 7.5, ymin = 0.2, ymax = 0.65)+
	geom_rect(aes(xmin = 4, xmax = 7.5, ymin = 0.2, ymax = 0.69), fill = NA, color = 'black',
		  linewidth = 0.5)+
	theme(legend.position = c(0.85, 0.5), 
	      legend.background = element_rect(fill = NA, colour = 'black'),
	      legend.title = element_blank())




ggplot(data = data.frame(x = c(times_simulation[is.finite(times_simulation)] / 120,
			       times_simulation2[is.finite(times_simulation2)] / 120,
			       c(tdet_noFil, tsto_noFil)/120),
			 type = c(rep('With recruitment', sum(is.finite(times_simulation))),
			 	 rep('Without recruitment', sum(is.finite(times_simulation2))),
			 	 rep('Experiment', length(c(tdet_noFil, tsto_noFil))))), 
       aes(x, fill = type)) +
	geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 35, 0.15),
		       color = 'black', alpha = 0.5, position = 'identity') + 
	geom_density(alpha = 0.5)+
	scale_x_continuous('Time (min)', breaks = seq(0, 10, 2.5),
			   limits = c(0, 10))+

	geom_vline(xintercept = c(median( c(tdet_noFil, tsto_noFil)/120),
				  median(times_simulation2[is.finite(times_simulation2)] / 120),
				  median(times_simulation[is.finite(times_simulation)] / 120)),
		   linewidth = 1, linetype = 2, color = c('deepskyblue4', 'mediumpurple', 'gold3'))+
	scale_y_continuous('Density', breaks = seq(0, 10, 0.15)) +
	scale_fill_manual('',values = c('deepskyblue4', 'gold3', 'mediumpurple'))+
	annotation_custom(
		ggplotGrob(ggplot(data = data.frame(x = c(times_simulation[is.finite(times_simulation)] / 120,
							 times_simulation2[is.finite(times_simulation2)] / 120,
							 c(tdet_noFil, tsto_noFil)/120),
							 type = c(rep('With recruitment', sum(is.finite(times_simulation))),
							 	 rep('Without recruitment', sum(is.finite(times_simulation2))),
							 	 rep('Experiment', length(c(tdet_noFil, tsto_noFil))))))+

			   	geom_violin(aes(type, x, fill = type), trim = FALSE,
			   		    alpha = 0.5, show.legend = FALSE)+
			   	theme(
			   		axis.ticks.x = element_blank())+xlab('')+
			   	scale_y_continuous('Time (min)', breaks = seq(0, 100, 10))+
			   	scale_fill_manual(values = c('deepskyblue4', 'gold3', 'mediumpurple'))),
		xmin = 3, xmax = 7.5, ymin = 0.4, ymax = 0.9)+
	geom_rect(aes(xmin = 3, xmax = 7.5, ymin = 0.4, ymax = 0.97), fill = NA, color = 'black',
		  linewidth = 0.5)+
	theme(legend.position = c(0.85, 0.5), 
	      legend.background = element_rect(fill = NA, colour = 'black'),
	      legend.title = element_blank())
