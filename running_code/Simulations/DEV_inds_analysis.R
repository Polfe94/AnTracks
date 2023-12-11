source('~/research/gits/AnTracks/src/Experiment.R')
library(arrow)
dt <- data.table(read_parquet('/home/polfer/research/gits/AutomatAnts/results/some_data.parquet'))[, -1]
keys <- data.table(read.csv('/home/polfer/research/gits/AutomatAnts/results/some_keys.csv'))[, -1]

parse_ids <- function(x){
	as.numeric(regmatches(x, gregexpr('\\d{1,2}', x))[[1]])
}

get_missing_ids <- function(x){
	x <- unique(unlist(x))
	if(length(x)) (0:99)[-(x+1)] else 0:99
}

ids <- dt[, .(id = parse_ids(id_out)), by = 'Frame']
ids2 <- dt[, .(id = list(parse_ids(id_out))), by = 'Frame']
ids_nest <- ids2[, .(id = get_missing_ids(id)), by = 'Frame']

Out <- merge(ids, keys, by = 'id', all = TRUE)
set(Out, j = 'type', value = 'Arena')
In <- merge(ids_nest, keys, by = 'id', all = TRUE)
set(In, j = 'type', value = 'Nest')

df <- rbind(Out, In)
df_means <- df[, .(g = mean(g)), by = c('Frame', 'type')]
# ggplot(data = df, aes(g, fill = factor(type))) +
# 	geom_histogram(alpha = 0.4, color = 'black',breaks = seq(0, 1, 0.02), position = 'dodge',
# 		       aes(y = after_stat(density)))+
# 	scale_fill_manual('', values = c('gold4', 'dodgerblue3'))+
# 	geom_density(alpha = 0.2)

# ggplot(data = df, aes(g, fill = factor(type))) +
# 	
# 	scale_fill_manual('', values = c('gold3', 'dodgerblue4'))+
# 	geom_density(alpha = 0.4) + 
# 	scale_x_continuous('Sensitivity (g)', breaks = seq(0, 1, 0.25))+
# 	scale_y_continuous('Relative frequency', breaks = function(i){
# 		seq(0, round(max(i)), length.out = 6)
# 	})

ggplot(data = df, aes(g, fill = factor(type))) +
	
	scale_fill_manual('', values = c('gold3', 'dodgerblue4'))+
	geom_histogram(alpha = 0.4, aes(y = after_stat(density)), position = 'identity',
		       breaks = seq(0, 1, 0.04), color = 'black') + 
	scale_x_continuous('Sensitivity (g)', breaks = seq(0, 1, 0.25))+
	scale_y_continuous('Relative frequency', breaks = function(i){
		seq(0, round(max(i)), length.out = 6)
	})+
	geom_density(alpha = 0.2)
# ggplot(data = df_means, aes(g, fill = factor(type))) +
# 	
# 	scale_fill_manual('', values = c('gold3', 'dodgerblue4'))+
# 	geom_density(alpha = 0.4) + 
# 	scale_x_continuous('Sensitivity (g)', breaks = seq(0, 1, 0.25))+
# 	scale_y_continuous('Relative frequency', breaks = function(i){
# 		seq(0, round(max(i)), length.out = 6)
# 	})

ggplot(data = df_means, aes(g, fill = factor(type))) +
	
	scale_fill_manual('', values = c('gold3', 'dodgerblue4'))+
	geom_histogram(alpha = 0.4, aes(y = after_stat(density)), position = 'identity',
		       breaks = seq(0, 1, 0.02), color = 'black') + 
	scale_x_continuous('Sensitivity (g)', breaks = seq(0, 1, 0.25))+
	scale_y_continuous('Relative frequency', breaks = function(i){
		seq(0, round(max(i)), length.out = 6)
	})+
	geom_density(alpha = 0.2, aes(color = factor(type)), show.legend = FALSE)+
	scale_color_manual('', values = c('gold3', 'dodgerblue4'))
 
# ids2 <- merge(ids, keys, by = 'id', all = TRUE)
# set(ids2, j = 'Frame', value = round(ids2[['T']] *2))
# 
# ids_unique <- ids2[, .(g = unique(g), id = unique(id)), by = 'Frame']
# 
# hist(ids2[['g']])
# hist(ids_unique[['g']])
# 
# 
# nest_ids <- merge(keys, ids_unique, by = 'id', all.x = TRUE)
# nest_ids <- CJ(frame = unique(ids_unique[['Frame']]), id = keys[['id']], g = keys[['g']])
# 

