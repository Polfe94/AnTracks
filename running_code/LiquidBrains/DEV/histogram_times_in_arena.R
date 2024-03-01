source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
library(arrow)
library(latex2exp)

path <- '/home/polfer/research/gits/AutomatAnts/results/with_recruitment/parameters/uniform/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('uniform_\\d{1,2}', files))))


afunc <- function(path, filename, filter = FALSE){

	food_sim <- data.table(read.csv(paste0(path, filename,'_food.csv')))
	dt <- data.table(read.csv(paste0(path, filename,'.csv')))[, c('Frame', 'pos')]
	if(filter){
		dt <- dt[Frame <= max(food_sim[['t']] * 2)]
	}
	dt_complete <- dt[, .(node = parse_nodes(pos), id = parse_ids(pos)), by = 'Frame'][!is.na(node)]
	
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
times_simulation_filtered <- c()
for(i in files){
	times_simulation_filtered <- c(times_simulation_filtered, afunc(path, i, TRUE))
}
Sys.time() - t0

ggplot(data = data.frame(x = times_simulation[is.finite(times_simulation)] / 120), aes(x)) +
	geom_histogram(breaks = seq(0, 35, 0.3), color = 'black', fill = 'gold3', alpha = 0.5) + 
	scale_x_continuous('Time (min)', breaks = seq(0, 40, 5),
							limits = c(0, 35))+
	geom_vline(xintercept = median(times_simulation[is.finite(times_simulation)] / 120),
linewidth = 1, linetype = 2)+
	ylab('Frequency')


load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/sto.RData')


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
	ids <- get_ids(exp)#dt[node %in% det_food_nodes, .(id = unique(N_ind))][['id']]
	sbst_exp <- dt[N_ind %in% ids]
	
	# if(filter){
	# 	cross_ids <- sbst_exp[, .(x = length(unique(Crossings)) == 1), by = 'N_ind']
	# 	filtered_ids <- sbst_exp[N_ind %in% cross_ids[x == TRUE, N_ind]]
	# } else {
	# 	filtered_ids <- sbst_exp
	# }
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
	print(paste0('Values below 0: ', sum(result <= 0)))
	result[result > 0]
}

tdet <- unlist(lapply(det, tracklets_to_food))
tsto <- unlist(lapply(sto, tracklets_to_food))

tdet_noFil <- unlist(lapply(det, tracklets_to_food, filter = FALSE))
tsto_noFil <- unlist(lapply(sto, tracklets_to_food, filter = FALSE))

ggplot(data = data.frame(v = c(tdet, tsto)/120, type = c(rep('det', length(tdet)),
							   rep('sto', length(tsto)))),
       # aes(v, fill = type)) + geom_histogram(alpha = 0.5, color = 'black', 
       # 				      bins = 100) + 
       aes(v, fill = type)) + geom_histogram(alpha = 0.5, color = 'black', 
       				      breaks = seq(0, 35, 0.3)) + 
	scale_fill_manual('', values = c('gold3', 'mediumpurple'), 
			  labels = c('Determinist', 'Stochastic'))+
	ylab('Frequency') +
	geom_vline(xintercept = median(c(tdet, tsto)/120), linetype = 2, linewidth = 1)+
	# geom_vline(xintercept = c(median(tdet/120),median(tsto/120)), linetype = 2, linewidth = 1)+
	scale_x_continuous('Time (min)', breaks = seq(0, 35, 5))+
	ggtitle('With filter (free tracklets without interactions)')

ggplot(data = data.frame(v = c(tdet_noFil, tsto_noFil)/120, type = c(rep('det', length(tdet_noFil)),
							 rep('sto', length(tsto_noFil)))),
       aes(v, fill = type)) + geom_histogram(alpha = 0.5, color = 'black', bins = 100) + 
	scale_fill_manual('', values = c('gold3', 'mediumpurple'), 
			  labels = c('Determinist', 'Stochastic'))+
	ylab('Frequency') +
	geom_vline(xintercept = median(c(tdet_noFil, tsto_noFil)/120), linetype = 2, linewidth = 1)+
	# geom_vline(xintercept = c(median(tdet_noFil/120),median(tsto_noFil/120)), linetype = 2, linewidth = 1)+
	scale_x_continuous('Time (min)', breaks = seq(0, 35, 5))+
	ggtitle('Without filter (full tracks with interactions)')

ggplot(data = data.frame(v = c(c(tdet, tsto)/120 ,c(tdet_noFil, tsto_noFil)/120),
			 type = c(rep('det', length(tdet)),
			 	 rep('sto', length(tsto)),
			 c(rep('det', length(tdet_noFil)),
			 	 rep('sto', length(tsto_noFil)))),
			 filter = c(rep('Filtered', length(c(tdet, tsto))),
			 	   rep('Unfiltered', length(c(tdet_noFil, tsto_noFil))))),
       aes(v, fill = type))+
	geom_histogram(alpha = 0.5, color = 'black', bins = 100) + 
	scale_fill_manual('', values = c('gold3', 'mediumpurple'), 
			  labels = c('Determinist', 'Stochastic'))+
	ylab('Frequency') +
	geom_vline(xintercept = median(c(tdet_noFil, tsto_noFil)/120), linetype = 2, linewidth = 1)+
	# geom_vline(xintercept = c(median(tdet_noFil/120),median(tsto_noFil/120)), linetype = 2, linewidth = 1)+
	scale_x_continuous('Time (min)', breaks = seq(0, 35, 5))+
	facet_wrap(~ filter)

ggplot(data = data.frame(v = c(c(tdet, tsto)/120 ,c(tdet_noFil, tsto_noFil)/120),
			 type = c(rep('det', length(tdet)),
			 	 rep('sto', length(tsto)),
			 	 c(rep('det', length(tdet_noFil)),
			 	   rep('sto', length(tsto_noFil)))),
			 filter = c(rep('Filtered', length(c(tdet, tsto))),
			 	   rep('Unfiltered', length(c(tdet_noFil, tsto_noFil))))),
       aes(x = v, fill = type))+
	geom_histogram(alpha = 0.5, color = 'black', bins = 100, aes(y = after_stat(density))) + 
	scale_fill_manual('', values = c('gold3', 'mediumpurple'), 
			  labels = c('Determinist', 'Stochastic'))+
	ylab('Frequency') +
	geom_vline(xintercept = median(c(tdet_noFil, tsto_noFil)/120), linetype = 2, linewidth = 1)+
	# geom_vline(xintercept = c(median(tdet_noFil/120),median(tsto_noFil/120)), linetype = 2, linewidth = 1)+
	scale_x_continuous('Time (min)', breaks = seq(0, 35, 5))+
	facet_wrap(~ filter)

ggplot(data = data.frame(v = c(c(tdet, tsto)/120 ,c(tdet_noFil, tsto_noFil)/120),
			 type = c(rep('det', length(tdet)),
			 	 rep('sto', length(tsto)),
			 	 c(rep('det', length(tdet_noFil)),
			 	   rep('sto', length(tsto_noFil)))),
			 filter = c(rep('Filtered', length(c(tdet, tsto))),
			 	   rep('Unfiltered', length(c(tdet_noFil, tsto_noFil))))),
       aes(v))+
	geom_histogram(alpha = 0.5, color = 'black', bins = 100) + 
	# scale_fill_manual('', values = c('gold3', 'mediumpurple'), 
	# 		  labels = c('Determinist', 'Stochastic'))+
	ylab('Frequency') +
	geom_vline(xintercept = median(c(tdet_noFil, tsto_noFil)/120), linetype = 2, linewidth = 1)+
	# geom_vline(xintercept = c(median(tdet_noFil/120),median(tsto_noFil/120)), linetype = 2, linewidth = 1)+
	scale_x_continuous('Time (min)', breaks = seq(0, 35, 5))+
	facet_wrap(~ filter)
	



## Estimated q95_SIM = 8.9125; q95_EXP = 8.625

aplot <- ggplot(data = data.frame(x = c(times_simulation[is.finite(times_simulation)] / 120,
			 c(tdet_noFil, tsto_noFil)/120),
			 type = c(rep('Simulation', sum(is.finite(times_simulation))),
			 	 rep('Experiment', length(c(tdet_noFil, tsto_noFil))))), 
       aes(x, fill = type)) +
	geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 35, 0.4),
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
	scale_y_continuous('Relative frequency', breaks = seq(0, 1, 0.15)) +
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



## Filtered (from t0 to last exploitation time)

bplot <- ggplot(data = data.frame(x = c(times_simulation_filtered[is.finite(times_simulation_filtered)] / 120,
			       c(tdet, tsto)/120),
			 type = c(rep('Simulation', sum(is.finite(times_simulation_filtered))),
			 	 rep('Experiment', length(c(tdet, tsto))))), 
       aes(x, fill = type)) +
	geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 35, 0.4),
		       color = 'black', alpha = 0.5, position = 'identity') + 
	geom_density(alpha = 0.5)+
	scale_x_continuous('Time (min)', breaks = seq(0, 10, 2.5),
			   limits = c(0, 10))+
	# geom_vline(xintercept = c(median(times_simulation_filtered[is.finite(times_simulation_filtered)] / 120),
	# 			  median( c(tdet, tsto)/120)),
	# 	   linewidth = 1.2, linetype = 2)+
	geom_vline(xintercept = c(median(times_simulation_filtered[is.finite(times_simulation_filtered)] / 120),
				  median( c(tdet, tsto)/120)),
		   linewidth = 1, linetype = 2, color = c('mediumpurple', 'gold3'))+
	scale_y_continuous('Relative frequency', breaks = seq(0, 1, 0.15)) +
	scale_fill_manual('',values = c('gold3', 'mediumpurple'))+
	annotation_custom(
		ggplotGrob(ggplot(data = data.frame(x = c(times_simulation_filtered[is.finite(times_simulation_filtered)] / 120,
							  c(tdet, tsto)/120),
							  type = c(rep('Simulation', sum(is.finite(times_simulation_filtered))),
							  	 rep('Experiment', length(c(tdet, tsto))))))+
			   	# geom_boxplot(aes(type, x, fill = type), width = 0.5,
			   	# 	     alpha = 0.5,outlier.shape = 1, show.legend = FALSE)+
			   	geom_violin(aes(type, x, fill = type),trim = FALSE,
			   		    alpha = 0.5, show.legend = FALSE)+
			   	theme(
			   		axis.ticks.x = element_blank())+xlab('')+
			   	scale_y_continuous('Time (min)', breaks = seq(0, 20, 2.5))+
			   	scale_fill_manual(values = c('gold3', 'mediumpurple'))),
		xmin = 4, xmax = 7.5, ymin = 0.2, ymax = 0.65)+
	geom_rect(aes(xmin = 4, xmax = 7.5, ymin = 0.2, ymax = 0.69), fill = NA, color = 'black',
		  linewidth = 0.5)+
	theme(legend.position = c(0.85, 0.5), 
	      legend.background = element_rect(fill = NA, colour = 'black'),
	      legend.title = element_blank())



ggarrange(aplot, bplot, ncol = 1, nrow = 2)




path2 <- '/home/polfer/research/gits/AutomatAnts/results/without_recruitment/parameters/uniform/'
files2 <- list.files(path2)
files2 <- unique(unlist(regmatches(files2, gregexpr('uniform_\\d{1,2}', files2))))

t0 <- Sys.time()
times_noRec <- c()
for(i in files2){
	times_noRec <- c(times_noRec, afunc(path2, i))
}
Sys.time() - t0



cplot <- ggplot(data = data.frame(x = c(times_noRec[is.finite(times_noRec)] / 120,
					c(tdet_noFil, tsto_noFil)/120),
					type = c(rep('Simulation', sum(is.finite(times_noRec))),
						 rep('Experiment', length(c(tdet_noFil, tsto_noFil))))), 
		aes(x, fill = type)) +
	geom_histogram(aes(y = after_stat(density)), breaks = seq(0, 35, 0.4),
		       color = 'black', alpha = 0.5, position = 'identity') + 
	geom_density(alpha = 0.5)+
	scale_x_continuous('Time (min)', breaks = seq(0, 10, 2.5),
			   limits = c(0, 10))+
	# geom_vline(xintercept = c(median(times_noRec[is.finite(times_noRec)] / 120),
	# 			  median( c(tdet_noFil, tsto_noFil)/120)),
	# 	   linewidth = 1.2, linetype = 2)+
	geom_vline(xintercept = c(median(times_noRec[is.finite(times_noRec)] / 120),
				  median( c(tdet_noFil, tsto_noFil)/120)),
		   linewidth = 1, linetype = 2, color = c('mediumpurple', 'gold3'))+
	scale_y_continuous('Relative frequency', breaks = seq(0, 1, 0.15)) +
	scale_fill_manual('',values = c('gold3', 'mediumpurple'))+
	annotation_custom(
		ggplotGrob(ggplot(data = data.frame(x = c(times_noRec[is.finite(times_noRec)] / 120,
							  c(tdet_noFil, tsto_noFil)/120),
							  type = c(rep('Simulation', sum(is.finite(times_noRec))),
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

cplot


ggarrange(aplot, cplot, ncol = 1)


############### INDIVIDUAL CHECKS ################


dt <- data.table(read.csv(paste0(path, 'uniform_10.csv')))[, c('Frame', 'pos')]
ids <- gsub('(', '', unlist(regmatches(dt[['pos']], gregexpr('\\(\\d{1,2},', dt[['pos']]))))
dt_filtered <- dt[!is.na(parse_nodes(dt[['pos']]))]
food_sim <- data.table(read.csv(paste0(path, 'uniform_10_food.csv')))

dt_complete <- rbindlist(lapply(seq_len(nrow(dt)), function(i){
	dt[i, .(Frame = Frame, node = parse_nodes(pos), id = parse_ids(pos))]
}))[!is.na(node)]

ids_filter <- dt_complete[node %in% food_sim[['node']], .(id = unique(id))][['id']]
dt_filtered <- dt_complete[id %in% ids_filter]

lalala <- dt_filtered[, .(mint = min(Frame)), by = c('node', 'id')][node %in% food_sim[['node']]]
lelele <- dt_filtered[, .(Frame = Frame), by = c('node','id')][node == '(0, 22)']

dt <- setDT(det[[1]]@data)

## criteria to fulfill: starts in nodes 662, 608, 634 or 635
nest_nodes <- c(662, 608, 634, 635)
det_food_nodes <- get_node(rbindlist(det[[1]]@food)[, 1:2], xy = hex[hex$y > 1000, ])
det_food_times <- get_times(det[[1]])#rbindlist(det[[1]]@food)[['t']]

ids <- get_ids(det[[1]])#dt[node %in% det_food_nodes, .(id = unique(N_ind))][['id']]
sbst_exp <- dt[N_ind %in% ids]

cross_ids <- sbst_exp[, .(x = length(unique(Crossings)) == 1), by = 'N_ind'] # 86
filtered_ids <- sbst_exp[Crossings != 0, .(N_ind = N_ind, 
					   node = node,
					   Frame = Frame)] # 64
filtered_ids <- sbst_exp[N_ind %in% cross_ids[x == TRUE, N_ind]]

visited_nodes <- filtered_ids[N_ind == 17 & node %in% det_food_nodes, .(visited_nodes = unique(node))][['visited_nodes']]
filtered_ids[N_ind == 17, .(maxt = max(Frame), mint = min(Frame)), by = node][node %in% c(662, 371)]
(3555-3086)/120 ## time in minutes

times_ind <- filtered_ids[N_ind == 17, .(mint = min(Frame), maxt = max(Frame)), by = node][node %in% c(662, visited_nodes)]


# example trajectory from one individual (ind = 17)
draw_hexagons(edges = edges) +
	geom_point(data = rbind(cbind(label = 'Complete', filtered_ids[N_ind == 17]),
				cbind(label = 'Path to food', filtered_ids[N_ind == 17 & Frame < 3555]),
				cbind(label = 'Path to nest', filtered_ids[N_ind == 17 & Frame > 3555])),
				aes(Xmm, Ymm, size = label, color = label))+
	scale_size_manual('', values = c(5, 3, 1.5))+ scale_color_manual('', values = c('black', 'purple', 'gold3'))+
	guides(size = guide_legend(override.aes = list(size = 5)))

# check 2
# times_ind <- filtered_ids[N_ind == 630, .(mint = min(Frame), maxt = max(Frame)), by = node][node %in% c(608, visited_nodes)]
# draw_hexagons(edges = edges) +
# 	geom_point(data = rbind(cbind(label = 'Complete', filtered_ids[N_ind == 630]),
# 				cbind(label = 'Path to food', filtered_ids[N_ind == 630 & Frame < 21111]),
# 				cbind(label = 'Path to nest', filtered_ids[N_ind == 630 & Frame > 21111])),
# 				aes(Xmm, Ymm, size = label, color = label))+
# 	scale_size_manual('', values = c(5, 3, 1.5))+ scale_color_manual('', values = c('black', 'purple', 'gold3'))+
# 	guides(size = guide_legend(override.aes = list(size = 5)))


# some trajectory examples...
draw_hexagons(edges = edges) + geom_point(data = filtered_ids[N_ind == 17], aes(Xmm, Ymm))
draw_hexagons(edges = edges) + geom_point(data = filtered_ids[N_ind == 630], aes(Xmm, Ymm))
draw_hexagons(edges = edges) + geom_point(data = filtered_ids[N_ind == 621], aes(Xmm, Ymm))
