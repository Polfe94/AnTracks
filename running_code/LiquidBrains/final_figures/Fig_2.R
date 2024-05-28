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

global_eff <- function(path){
	files <- list.files(path)[grepl('food', list.files(path))]
	results <- lapply(files, function(i){
		f <- gsub('_food', '', i)
		d <- data.table(read_parquet(paste0(path, f)))
		m1 <- min(d[N >0, Frame])/2
		food <- data.table(read_parquet(paste0(path, i)))
		1/data.table(mint = min(food[['t']]) - m1,
			     maxt = max(food[['t']])-min(food[['t']]))
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

pos_data <- rbindlist(lapply(files, function(i){ ## <--- takes a couple of minutes and a bit of RAM
	do.call('gc', args = list(verbose = FALSE))
	data.table(read_parquet(paste0(path, i, '_data.parquet')))[, .(node = parse_nodes(pos))][, .N, by = 'node']
}))

sum_pos_data <- merge(hex_sim, pos_data[, .(N = sum(N)), by = 'node'], by = 'node')
c <- draw_nodes(sum_pos_data, 
		add = draw_hexagons(sim_edges, linewidth = 2, color = 'black',
				    add = draw_hexagons(sim_edges, linewidth = 2, color = 'black', 
				    		    add = geom_foodpatches(food = food, fill = 'mediumpurple', alpha = 0.9))), 
		z = rank(sum_pos_data[['N']]), size = 4, show.legend = FALSE) + 
	theme(aspect.ratio = 0.5 , plot.margin = margin(t = -5, b = -15), plot.title = element_text(size = 22))+
	scale_fill_viridis(option = 'plasma') +
	ggtitle('','Model')

# POS <- rbindlist(lapply(seq_along(files), function(i){
# 	data.table(read_parquet(paste0(path, files[i], '_positions.parquet')))
# }))
# POS_trimmed <- rbindlist(lapply(seq_along(files), function(i){
# 	data.table(read_parquet(paste0(path, files[i], '_positions.parquet')))
# }))

# c <- draw_nodes(pos, 
# 	   add = draw_hexagons(sim_edges, linewidth = 2, color = 'black',
# 	   		    add = draw_hexagons(sim_edges, linewidth = 2, color = 'black', 
# 	   		    		    add = geom_foodpatches(food = food, fill = 'mediumpurple', alpha = 0.9))), 
# 	   z = pos[['z']], size = 4, show.legend = FALSE) + 
# 	theme(aspect.ratio = 0.5)+
# 	scale_fill_viridis(option = 'plasma')

pos_det <- merge(hex, rbindlist(lapply(det, function(i){
	i@data[, .N, by = 'node']
}), idcol = TRUE)[, .(N = sum(N)), by = 'node'], by = 'node')

# pos_det_trimmed <- data.table(merge(hex[hex$y > 1000, ], rbindlist(lapply(det, function(i){
# 	x <- min(rbindlist(i@food)[['t']])
# 	i@data[Frame <= x, .N, by = 'node']
# }), idcol = TRUE)[, .(N = sum(N)), by = 'node'], by = 'node', all = TRUE))[is.na(N), N := 0]
# 
# draw_nodes(pos_det_trimmed, 
# 	   add = draw_hexagons(edges, linewidth = 2, color = 'black',
# 	   		    add = draw_hexagons(edges, linewidth = 2, color = 'black', 
# 	   		    		    add = geom_foodpatches(fill = 'mediumpurple', alpha = 0.9))), 
# 	   z = rank(pos_det_trimmed[['N']]), size = 4) + 
# 	scale_fill_viridis('Occupancy',option = 'plasma', breaks = range(rank(pos_det_trimmed[['N']])),
# 			   labels = c('Low', 'High')) +
# 	theme(legend.position = 'bottom', aspect.ratio = 0.5,
# 	      legend.margin = margin(t = -20, unit = 'pt'))+
# 	guides(fill = guide_colorbar(title.position = 'top', barwidth = 15, title.hjust = 0.5))

c2 <- draw_nodes(pos_det, 
		 add = draw_hexagons(edges, linewidth = 2, color = 'black',
		 		    add = draw_hexagons(edges, linewidth = 2, color = 'black', 
		 		    		    add = geom_foodpatches(fill = 'mediumpurple', alpha = 0.9))), 
		 z = rank(pos_det[['N']]), size = 4) + 
	scale_fill_viridis('Density',option = 'plasma', breaks = c(1, 620),
			   labels = c('Low', 'High')) +
	theme(legend.position = 'bottom', aspect.ratio = 0.5,
	      plot.margin = margin(t = 0, b = -20), plot.title = element_text(size = 22))+
	guides(fill = guide_colorbar(title.position = 'top', barwidth = 15, title.hjust = 0.5))+
	ggtitle('B', 'DET')


# scico::scale_fill_scico()

#### POPULATION DYNAMICS ####

data <- process_data(path)

full_data <- rbindlist(lapply(files, function(i){
	data.table(read_parquet(paste0(path, i,'.parquet')))
}), idcol = TRUE)[order(Frame)]
full_data_avg <- full_data[Frame < 175*120, .(N = mean(N)), by = 'Frame'] ## filtered

foods <- rbindlist(lapply(files, function(i){
	food <- data.table(read_parquet(paste0(path, i,'_food.parquet')))[['t']]
	data <- data.table(read_parquet(paste0(path, i,'.parquet')))
	mint <- data[N > 0, min(Frame)] /2
	r <- range(food) - mint
	data.table(mint = r[1], maxt = r[2])
	
}), idcol = TRUE)[is.finite(mint),lapply(.SD, median), .SDcols = c('mint', 'maxt')]


data_peak <- rbindlist(data, idcol = TRUE)

# comparison of population dynamics
a <- ggplot(data = data_peak[Frame < 175*120], aes(Frame, N)) +
	geom_path(aes(group = .id), color = 'grey80', alpha = 0.2, linewidth = 0.8)+
	geom_path(data = det_N[Frame < 175*120], aes(group = .id), color = 'grey30', alpha = 0.5, linewidth = 0.8)+
	
	geom_path(data = full_data_avg, color = 'black', linewidth = 3, alpha = 0.5)+
	geom_path(data = full_data_avg, aes(color = 'Simulations'), linewidth = 2, alpha = 0.5)+
	geom_path(data = det_avg[Frame < 175*120], color = 'black', linewidth = 3, alpha = 0.5)+
	geom_path(data = det_avg[Frame < 175*120], aes(color = 'Experiments'), linewidth = 2, alpha = 0.5)+
	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50)*120, labels = seq(0, 150, 50))+
	scale_y_continuous('N (number of ants in arena)')+
	scale_color_manual('',values = c('grey30', 'grey80'), labels = c('DET', 'Model'))+
	geom_vline(xintercept = unlist(foods)*2, linetype = 2, linewidth = 1)+
	geom_vline(xintercept = unlist(foods_det)*2, linetype = 3, linewidth = 1)+
	theme(legend.position = c(0.8, 0.85),
	      legend.background = element_rect(fill = 'white', color = 'black'),
	      legend.title = element_blank(),plot.title = element_text(size = 22))+
	ggtitle('A')


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



round(wilcox.test(formula = value ~ condition,data = data.table(melt(rbindlist(list(coll_dt, 
			       cbind(condition = 'Experiments',
			             rbindlist(det_tps))), use.names = TRUE), 
		id.vars = 'condition'))[variable == 'tp1'])$p.value, 3)

data_b <- data.table(rbind(melt(rbind(global_sim_effs[, -'.id'], det_global_eff), id.vars = 'condition'),
		data.table(melt(rbindlist(list(coll_dt, 
					       cbind(condition = 'Experiments',
					             rbindlist(det_tps))), use.names = TRUE), 
				id.vars = 'condition'))[, value := 1/value]))
vars <- c('mint', 'maxt', 'tp1', 'tp2')
pvals <- data.table(value = round(vapply(vars,
					 function(i){wilcox.test(formula = value ~ condition, 
					 			data = data_b[variable == i])$p.value},
					 numeric(1)),
				  3), 
		    variable = vars)
pvals[['y']] <- 0.8*vapply(vars, function(i) max(data_b[variable == i, value], na.rm = T), numeric(1))
pvals[['condition']] <- NA
pvals[1:2, 'y'] <- pvals[2:1, y]

lbl_function <- function(label){
	l <- paste0('p = ', label)
	ifelse(nchar(label) == 5, l, paste0(l, '0'))
}


b <- ggplot(data = data_b, aes(condition, value, fill = condition))+
	
	geom_boxplot(alpha = 0.6, outlier.shape = NA, show.legend = FALSE) +
	geom_jitter(alpha = 0.25, width = 0.15, size = 3, show.legend = FALSE)+
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation', 
						 'Exploration (per capita)', 'Exploitation (per capita)')),
		   scales = 'free_y')+
	scale_x_discrete('', labels = c('DET', 'Model')) + 
	ylab(TeX('Search efficiency ($s^{-1}$)'))+
	geom_segment(data = pvals, aes(x = 'Experiments', xend = 'Simulations', y = y, yend = y,
				       group = factor(variable, levels = vars)))+
	geom_text(data = pvals, aes(x = 1.5, y = y*1.1, 
				    group = factor(variable, levels = vars), 
				    label = lbl_function(value)), size = 5) +
	scale_fill_manual('', values = c('mediumpurple', 'brown4'))+
	# scale_fill_manual('', values = c('mediumpurple', 'gold3'))+
	theme(legend.key.size = unit(30, 'pt'),
	      plot.title = element_text(size = 22))+
	ggtitle('C')

png('/home/polfer/research/gits/AnTracks/plots/figures_LiquidBrains/Fig_2.png', 8000, 8000, res = 600)
grid.arrange(grobs = list(a + theme(aspect.ratio = 0.5), b, 
			  ggarrange(c2, c, ncol = 1, common.legend = TRUE, legend = 'bottom')),
			  layout_matrix = rbind(rep(1, 10), rep(1, 10), 
			  		      c(rep(3, 4), rep(2, 6)), 
			  		      c(rep(3, 4), rep(2, 6))))
dev.off()

# png('/home/polfer/research/gits/AnTracks/plots/summary_default.png', 8000, 6000, res = 550)
# grid.arrange(grobs = list(a + theme(aspect.ratio = 0.5), b, 
# 			  ggarrange(c2, c, ncol = 1, common.legend = TRUE, legend = 'bottom')),
# 			  layout_matrix = rbind(c(1, 2), c(1, 2), c(3, 2), 
# 			  		      c(3, NA), c(3, NA), c(3, NA)))
# dev.off()


########################################################
############ FALTA MUTUAL INFORMATION !!!!! ############
########################################################

load("/home/polfer/research/gits/AnTracks/results/MI_det.RData")
load("/home/polfer/research/gits/AnTracks/results/MI_default.RData")
mi_sim <- rbindlist(lapply(seq_along(mutual_info_result), function(i){
	x <- mutual_info_result[[i]]
	n <- names(x)
	rbindlist(lapply(seq_along(x), function(ii){
		data.frame(t = as.numeric(n)[ii], x = as.numeric(x[[ii]]))
	}))
}), idcol = TRUE)

mi_exp <- rbindlist(lapply(seq_along(det_mi), function(i){
	
	x <- det_mi[[i]]
	n <- names(x)
	rbindlist(lapply(seq_along(x), function(ii){
		data.frame(t = as.numeric(n)[ii], x = as.numeric(x[[ii]]))
	}))
	
	
}), idcol = TRUE)


# ggplot(data = mi_sim, aes(t, x, group = t)) + geom_boxplot()
# ggplot(data = mi_sim[, .(x = mean(x)), by = 't'], aes(t, x))+
# 	geom_point(size = 3, shape = 21) +
# 	geom_path()+
# 	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50)*120,
# 			   labels = seq(0, 150, 50))+
# 	scale_y_continuous('<MI>')
# ggplot(data = mi_exp[, .(x = mean(x)), by = 't'], aes(t, x))+
# 	geom_point(size = 3, shape = 21) +
# 	geom_path()+
# 	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50)*120,
# 			   labels = seq(0, 150, 50))+
# 	scale_y_continuous('<MI>')


d <- ggplot(data = mi_sim[, .(x = mean(x)), by = 't'][!is.na(x), .(x = norm_range(x, 0, 1), t = t)], aes(t, x), alpha = 0.5)+
	geom_point(size = 3, shape = 21, aes(color = 'Simulations')) +
	geom_path(aes(color = 'Simulations'))+
	geom_point(size = 3, shape = 21, 
		   data = mi_exp[, .(x = mean(x)), by = 't'][!is.na(x), .(x = norm_range(x, 0, 1), t = t)], aes(color = 'Experiments')) +
	geom_path(data = mi_exp[, .(x = mean(x)), by = 't'][!is.na(x), .(x = norm_range(x, 0, 1), t = t)], aes(color = 'Experiments'))+
	geom_vline(xintercept = unlist(foods)*2, linetype = 2, linewidth = 1)+
	geom_vline(xintercept = unlist(foods_det)*2, linetype = 3, linewidth = 1)+
	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50)*120,
			   labels = seq(0, 150, 50))+
	scale_y_continuous(TeX('$<MI_{sim}>$'), labels = seq(0, 1, length.out = 5)*
			   	round(max(mi_sim[, .(x = mean(x)), by = 't'][['x']], na.rm = TRUE), 4),
			   sec.axis = sec_axis(trans = function(i) i, 
			   		    labels = seq(0, 1, length.out = 5)*
			   		    	round(max(mi_exp[, .(x = mean(x)), by = 't'][['x']], na.rm = TRUE), 4),
			   		    name = TeX('$<MI_{exp}>$'))) +
	scale_color_manual('', values = c('mediumpurple', 'gold3'))+
	theme(legend.position = c(0.825, 0.85), 
	      legend.background = element_rect(color = 'black', fill = NA),
	      legend.title = element_blank())



# png('/home/polfer/research/gits/AnTracks/plots/summary_default_prov.png', 8000, 6000, res = 550)
# grid.arrange(grobs = list(a + theme(aspect.ratio = 0.5, plot.title = element_text(face = 'bold'))+
# 			  	ggtitle('A'), b+ggtitle('B')+theme(plot.title = element_text(face = 'bold')), 
# 			  ggarrange(c2, c, ncol = 1, common.legend = TRUE,
# 			  	  legend = 'bottom', labels = c('C'), hjust = -8.25, vjust = 0.5 ), 
# 			  d + ggtitle('D')+theme(plot.title = element_text(face = 'bold'))),
# 			  layout_matrix = rbind(c(1, 2), c(1, 2), c(3, 2), 
# 			  		      c(3, 2), c(3, 4), c(3, 4)))
# dev.off()
# 
# 
# 
# png('/home/polfer/research/gits/AnTracks/plots/summary_default_prov.png', 8000, 6000, res = 550)
# grid.arrange(grobs = list(a + theme(aspect.ratio = 0.5, plot.title = element_text(face = 'bold'))+
# 			  	ggtitle('A'), b+ggtitle('B')+theme(plot.title = element_text(face = 'bold')), 
# 			  ggarrange(c2, c, ncol = 1, common.legend = TRUE,
# 			  	  legend = 'bottom', labels = c('C'), hjust = -8.25, vjust = 0.5 ), 
# 			  d + ggtitle('D')+theme(plot.title = element_text(face = 'bold'))),
# 			  layout_matrix = rbind(c(1, 2), c(1, 2), c(3, 2), 
# 			  		      c(3, 2), c(3, 4), c(3, 4)))
# dev.off()


png('/home/polfer/research/gits/AnTracks/plots/summary_default_prov.png', 8000, 6000, res = 550)
grid.arrange(grobs = list(a + theme(aspect.ratio = 0.5, plot.title = element_text(face = 'bold'))+
			  	ggtitle('A'), b+ggtitle('B')+theme(plot.title = element_text(face = 'bold')), 
			  ggarrange(c2, c, ncol = 1, common.legend = TRUE,
			  	  legend = 'bottom', labels = c('C'), hjust = -8.25, vjust = 0.5 ), 
			  d + ggtitle('D')+theme(plot.title = element_text(face = 'bold'))),
			  layout_matrix = rbind(c(1, 2), c(1, 2), c(3, 2), 
			  		      c(3, 2), c(3, 4), c(3, 4)))
dev.off()


ggsave('/home/polfer/research/gits/AnTracks/plots/summary_default.pdf', 
       plot = grid.arrange(grobs = list(a + theme(aspect.ratio = 0.5, plot.title = element_text(face = 'bold'))+
       				 	ggtitle('A'), b+ggtitle('B')+theme(plot.title = element_text(face = 'bold')), 
       				 ggarrange(c2, c, ncol = 1, common.legend = TRUE,
       				 	  legend = 'bottom', labels = c('C'), hjust = -8.25, vjust = 0.5 ), 
       				 d + ggtitle('D')+theme(plot.title = element_text(face = 'bold'))),
       				 layout_matrix = rbind(c(1, 2), c(1, 2), c(3, 2), 
       				 		      c(3, 2), c(3, 4), c(3, 4))), device = 'pdf',
       dpi = 300, width = 4400, height = 3300, unit = 'px')


