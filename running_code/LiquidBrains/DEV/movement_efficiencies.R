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

path <- '/home/polfer/research/gits/AutomatAnts/results/2024/modified_mot/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('modified_mot_\\d{1,2}', files))))

path <- '/home/polfer/research/gits/AutomatAnts/results/2024/random/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('random_\\d{1,2}', files))))

path <- '/home/polfer/research/gits/AutomatAnts/results/2024/nest_avoid/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('ballistic_\\d{1,2}', files))))

path <- '/home/polfer/research/gits/AutomatAnts/results/2024/nest_avoid_bias_1.5/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('ballistic_\\d{1,2}', files))))

foods <- rbindlist(lapply(files, function(i){
	data.table(read_parquet(paste0(path, i,'_food.parquet')))[, c('t')]
}), idcol = TRUE)[, .(mint = min(t), maxt = max(t)), by = '.id'][is.finite(mint),
								 lapply(.SD, median), .SDcols = c('mint', 'maxt')]

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

ggplot(data = rbind(melt(rbind(global_sim_effs[, -'.id'], det_global_eff), id.vars = 'condition'),
		    data.table(melt(rbindlist(list(coll_dt, 
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