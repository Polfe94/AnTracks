#### LIBRARIES, DATA AND GENERIC FUNCTIONS ####
library(dtw)
library(dtwclust)
library(latex2exp)
library(arrow)
# library(infotheo)
source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')

load('~/research/gits/AnTracks/data/det.RData')

path <- '/home/polfer/research/gits/AutomatAnts/results/2024/default/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('default_\\d{1,2}', files))))

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

#### SPATIAL DYNAMICS
pos <- data.table(read_parquet(paste0(path, files[3], '_positions.parquet')))
data <- data.table(read_parquet(paste0(path, files[3], '_data.parquet')))
food <- data.table(read_parquet(paste0(path, files[3], '_food.parquet')))[, node := revert_node(node)]

foods <- lapply(seq_along(files), function(i){
	data.table(read_parquet(paste0(path, files[i], '_food.parquet')))[, node := revert_node(node)][, t := round(t * 2)]
})

ts <- round(food[['t']]*2)
dt <- data[, .(node = parse_nodes(pos), id = parse_ids(id_out)), by = c('Frame')]
dt_filtered <- dt[Frame <= min(ts)]

data_full <- merge(hex_sim, dt[, .N, by = 'node'], all = TRUE)[is.na(N), N := 0]
data_show <- merge(hex_sim, dt_filtered[, .N, by = 'node'], all = TRUE)[is.na(N), N := 0]


foodpatches <- as.data.frame(merge(food, hex_sim, by = 'node', sort = FALSE)[, c('x', 'y')])
foodpatches <- list(GP1 = foodpatches[1:6, ], GP2 = foodpatches[7:12, ])

draw_nodes(data_show,
	   add = draw_hexagons(sim_edges, linewidth = 2, color = 'black',
	   		    add = draw_hexagons(sim_edges, linewidth = 2, color = 'black',
	   		    		    add = geom_foodpatches(food = foodpatches, fill = 'mediumpurple', alpha = 0.9))),
	   z = rank(data_show[['N']]), size = 4) +
	scale_fill_viridis('Occupancy',option = 'plasma', breaks = range(rank(data_show[['N']])),
			   labels = c('Low', 'High')) +
	theme(legend.position = 'bottom', aspect.ratio = 0.5,
	      legend.margin = margin(t = -20, unit = 'pt'))+
	guides(fill = guide_colorbar(title.position = 'top', barwidth = 15, title.hjust = 0.5))

draw_nodes(data_full,
	   add = draw_hexagons(sim_edges, linewidth = 2, color = 'black',
	   		    add = draw_hexagons(sim_edges, linewidth = 2, color = 'black',
	   		    		    add = geom_foodpatches(food = foodpatches, fill = 'mediumpurple', alpha = 0.9))),
	   z = rank(data_full[['N']]), size = 4) +
	scale_fill_viridis('Occupancy',option = 'plasma', breaks = range(rank(data_full[['N']])),
			   labels = c('Low', 'High')) +
	theme(legend.position = 'bottom', aspect.ratio = 0.5,
	      legend.margin = margin(t = -20, unit = 'pt'))+
	guides(fill = guide_colorbar(title.position = 'top', barwidth = 15, title.hjust = 0.5))


# POS <- rbindlist(lapply(seq_along(files), function(i){
# 	data.table(read_parquet(paste0(path, files[i], '_positions.parquet')))
# }))
# POS_trimmed <- rbindlist(lapply(seq_along(files), function(i){
# 	data.table(read_parquet(paste0(path, files[i], '_positions.parquet')))

bigdata <- vector('list', length(files))
for(i in seq_along(bigdata)){
	bigdata[[i]] <- data.table(read_parquet(paste0(path, files[i], '_data.parquet')))[N > 0, c('pos','id_out', 'Frame')]
}

bigdata <- rbindlist(bigdata, idcol = TRUE)[pos != '']

foods <- rbindlist(foods, idcol = TRUE)

pos_data <- bigdata[, .(node = parse_nodes(pos)), by = c('.id', 'Frame')]

explor_pos <- rbindlist(lapply(seq_len(max(foods[, .id])), function(i){
pos_data[.id == i & Frame <=  min(foods[.id == i, t]), .N, by = 'node']}), idcol =TRUE)  


#### EXPLORATION SIMULATIONS

load('~/research/gits/AnTracks/results/explor_pos.RData')
load('~/research/gits/AnTracks/results/explor_pos_random.RData')
load('~/research/gits/AnTracks/results/explor_pos_ballistic.RData')

food <- as.data.frame(merge(food, hex_sim, by = 'node', sort = FALSE)[, c('x', 'y')])
food <- list(GP1 = food[1:6, ], GP2 = food[7:12, ])

explor_pos <- merge(hex_sim, explor_pos[, .(N = sum(N)), by = 'node'], all = TRUE)[is.na(N), N := 0]
draw_nodes(explor_pos,
	   add = draw_hexagons(sim_edges, linewidth = 2, color = 'black',
	   		    add = draw_hexagons(sim_edges, linewidth = 2, color = 'black',
	   		    		    add = geom_foodpatches(food = food, fill = 'mediumpurple', alpha = 0.9))),
	   z = rank(explor_pos[['N']]), size = 4) +
	scale_fill_viridis('Occupancy',option = 'plasma', breaks = range(rank(explor_pos[['N']])),
			   labels = c('Low', 'High')) +
	theme(legend.position = 'bottom', aspect.ratio = 0.5,
	      legend.margin = margin(t = -20, unit = 'pt'))+
	guides(fill = guide_colorbar(title.position = 'top', barwidth = 15, title.hjust = 0.5))


#### FULL SIMULATIONS

load('~/research/gits/AnTracks/results/full_pos.RData')

full_pos <- merge(hex_sim, full_pos[, .(N = sum(N)), by = 'node'], all = TRUE)[is.na(N), N := 0]
draw_nodes(full_pos,
	   add = draw_hexagons(sim_edges, linewidth = 2, color = 'black',
	   		    add = draw_hexagons(sim_edges, linewidth = 2, color = 'black',
	   		    		    add = geom_foodpatches(food = foodpatches, fill = 'mediumpurple', alpha = 0.9))),
	   z = rank(full_pos[['N']]), size = 4) +
	scale_fill_viridis('Occupancy',option = 'plasma', breaks = range(rank(full_pos[['N']])),
			   labels = c('Low', 'High')) +
	theme(legend.position = 'bottom', aspect.ratio = 0.5,
	      legend.margin = margin(t = -20, unit = 'pt'))+
	guides(fill = guide_colorbar(title.position = 'top', barwidth = 15, title.hjust = 0.5))



### DETERMINISTS 

pos_det <- merge(hex, rbindlist(lapply(det, function(i){
	i@data[, .N, by = 'node']
}), idcol = TRUE)[, .(N = sum(N)), by = 'node'], by = 'node')

pos_det_trimmed <- data.table(merge(hex[hex$y > 1000, ], rbindlist(lapply(det, function(i){
	x <- min(rbindlist(i@food)[['t']])
	i@data[Frame <= x, .N, by = 'node']
}), idcol = TRUE)[, .(N = sum(N)), by = 'node'], by = 'node', all = TRUE))[is.na(N), N := 0]

draw_nodes(pos_det_trimmed,
	   add = draw_hexagons(edges, linewidth = 2, color = 'black',
	   		    add = draw_hexagons(edges, linewidth = 2, color = 'black',
	   		    		    add = geom_foodpatches(fill = 'mediumpurple', alpha = 0.9))),
	   z = rank(pos_det_trimmed[['N']]), size = 4) +
	scale_fill_viridis('Occupancy',option = 'plasma', breaks = range(rank(pos_det_trimmed[['N']])),
			   labels = c('Low', 'High')) +
	theme(legend.position = 'bottom', aspect.ratio = 0.5,
	      legend.margin = margin(t = -20, unit = 'pt'))+
	guides(fill = guide_colorbar(title.position = 'top', barwidth = 15, title.hjust = 0.5))




### NO FOOD
load('~/research/gits/AnTracks/data/nf.RData')

pos_nf <- data.table(merge(hex[hex$y > 1000, ], rbindlist(lapply(nf, function(i){
	setDT(i@data)[, .N, by = 'node']
}), idcol = TRUE)[, .(N = sum(N)), by = 'node'], by = 'node', all = TRUE))[is.na(N), N := 0]

# pos_nf_trimmed <- data.table(merge(hex[hex$y > 1000, ], rbindlist(lapply(nf, function(i){
# 	x <- min(rbindlist(i@food)[['t']])
# 	i@data[Frame <= x, .N, by = 'node']
# }), idcol = TRUE)[, .(N = sum(N)), by = 'node'], by = 'node', all = TRUE))[is.na(N), N := 0]

draw_nodes(pos_nf,
	   add = draw_hexagons(edges, linewidth = 2, color = 'black',
	   		    add = draw_hexagons(edges, linewidth = 2, color = 'black',
	   		    		    add = geom_foodpatches(fill = 'mediumpurple', alpha = 0.9))),
	   z = rank(pos_nf[['N']]), size = 4) +
	scale_fill_viridis('Occupancy',option = 'plasma', breaks = range(rank(pos_nf[['N']])),
			   labels = c('Low', 'High')) +
	theme(legend.position = 'bottom', aspect.ratio = 0.5,
	      legend.margin = margin(t = -20, unit = 'pt'))+
	guides(fill = guide_colorbar(title.position = 'top', barwidth = 15, title.hjust = 0.5))
