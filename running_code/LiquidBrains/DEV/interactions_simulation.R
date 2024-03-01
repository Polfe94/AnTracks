#### LIBRARIES, DATA AND GENERIC FUNCTIONS ####
library(dtw)
library(dtwclust)
library(latex2exp)
library(arrow)
library(infotheo)
source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')

parse_nodes <- function(nodes){
	unlist(strsplit(nodes, ';'))
}

parse_ids <- function(ids){
	as.integer(unlist(strsplit(ids, ',')))
}

path <- '/home/polfer/research/gits/AutomatAnts/results/2024/Jij_0.01/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('Jij_0.01_\\d{1,2}', files))))
sim_data <- sim_data <- data.table(read_parquet(paste0(path, files[3], '_data.parquet')))

## spatial
sum_sim <- sim_data[, .(node = parse_nodes(pos), id = parse_ids(id_out)), by = 'Frame'][order(Frame)]
duplicates <- sum_sim[, .N, by = c('Frame', 'node', 'id')]
ggplot(data = duplicates, aes('', N)) + geom_boxplot()+
	ylab('Number of ants in a single node') + xlab('') +theme(axis.ticks.x = element_blank())
interactions <- merge(hex_sim[, c('node', 'x', 'y')],
		      duplicates[, .(z = sum(N)), by = 'node'], by = 'node', all = TRUE)[is.na(z), z:= 0]
interactions[z >= quantile(z, p = 0.95), 'z'] <- quantile(interactions[['z']], p = 0.95)
draw_nodes(interactions, z = rank(interactions[['z']]), size = 4) +
	scale_fill_viridis(option = 'C')

## temporal
tI <- duplicates[N > 1, .(I = sum(N)), by = c('Frame')]
rbindlist(list(tI[, .(Frame = Frame, N = I, norm = norm_range(I, 0, 1))], 
	       sim_data[, .(Frame = Frame, N, norm = norm_range(N, 0, 1))]), idcol = TRUE)

ggplot(data = rbindlist(list(tI[, .(Frame = Frame, N = I, norm = norm_range(I, 0, 1))], 
			     sim_data[, .(Frame = Frame, N, norm = norm_range(N, 0, 1))]), idcol = TRUE),
       aes(Frame, norm))  +geom_path(aes(color = factor(.id, labels = c('Interactions', 'N'))), alpha = 0.75)+
       	scale_color_manual('', values = c('mediumpurple', 'gold3'))+
	scale_x_continuous('Time (min)', breaks = seq(0, 150, 50)*120,
			   labels = seq(0, 150, 50))+
	scale_y_continuous('N', sec.axis = sec_axis(trans = ~. *1,
			   		    name = 'Interactions'))
	
			    # 			   		    breaks = seq(0, 60, length.out = 9),
			    # 			   		    labels = seq(0, 2000, length.out = 9),
			    # 			   		    name = 'Interactions (cumulative sum)'))+


## Si
