source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
library(arrow)
library(latex2exp)
path <- '/home/polfer/research/gits/AutomatAnts/results/rec_noRec_nf/rec/uniform/'
dt <- data.table(read_parquet(paste0(path, 'uniform_0_data.parquet')))
n <- data.table(read_parquet(paste0(path, 'uniform_0.parquet')))
pos <- data.table(read_parquet(paste0(path, 'uniform_0_positions.parquet')))
food <- data.table(read_parquet(paste0(path, 'uniform_0_food.parquet')))
keys <- data.table(read_parquet(paste0(path, 'uniform_0_keys.parquet')))

load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/sto.RData')


parse_ids <- function(x){
	as.numeric(strsplit(x, ',')[[1]])
}

# parse_ids <- function(x){
# 	as.numeric(unlist(strsplit(x, ',')))
# }

get_time <- function(x, threshold = 10){
	if(is.unsorted(x)){
		x <- sort(x)
	}
	y <- diff(x)
	y <- y[y <= threshold]
	sum(y)
}

ind_efficiency <- function(path, filename){
	if(file.exists(paste0(path, filename, '_data.parquet'))){
		dt <- data.table(read_parquet(paste0(path, filename, '_data.parquet')))
		food <- data.table(read_parquet(paste0(path, filename ,'_food.parquet')))
		food <- unlist(food[!is.na(t), .(t = range(t))], use.names = FALSE)
		if(all(is.na(food))){
			food <- rep(10800, 2)
		} else if(food[1] == food[2]){
			food[2] <- 10800
		}
		dt <- dt[id_out != '' & `T` <= food[2], .(id = parse_ids(id_out)), by = 'T']
		tp1 <- sum(dt[`T` <= food[1], .(time = get_time(`T`)), by = 'id'][['time']])
		tp2 <- sum(dt[, .(time = get_time(`T`)), by = 'id'][['time']])
		data.table(tp1 = tp1, tp2 = tp2)
	} else {
		warning(paste0(filename, ' does not exist'))
		}
}

ind_efficiency_v2 <- function(path, filename, minpieces = 12){
	if(file.exists(paste0(path, filename, '_data.parquet'))){
		food <- data.table(read_parquet(paste0(path, filename ,'_food.parquet')))
		if(length(food[!is.na(t), t]) > (minpieces-1)){
			food <- unlist(food[!is.na(t), .(t = range(t))], use.names = FALSE)
			dt <- data.table(read_parquet(paste0(path, filename, '_data.parquet')))
			dt <- dt[id_out != '' & `T` <= food[2], .(id = parse_ids(id_out)), by = 'T']
			tp1 <- sum(dt[`T` < food[1], .(time = get_time(`T`)), by = 'id'][['time']])
			tp2 <- sum(dt[`T` > food[1], .(time = get_time(`T`)), by = 'id'][['time']])
			data.table(tp1 = tp1, tp2 = tp2)
		}

	} else {
		warning(paste0(filename, ' does not exist'))
	}
}

## INDIVIDUAL TESTING
draw_hexagons(edges = sim_edges, color = 'grey50', linewidth = 1.2)+
	geom_point(data = pos, aes(x, y, fill = z), shape = 21,size = 4)+
	scale_fill_viridis_c(option = 'C') + 
	theme(legend.position = 'none', aspect.ratio = 0.5)

ggplot(data = n, aes(Frame, N)) + geom_path() + 
	geom_vline(xintercept = range(food[['t']]), linetype = 2, linewidth = 1) +
	scale_x_continuous('Time (min)', limits = c(0, 21600), 
			   breaks = seq(0, 21600, length.out = 7),
			   labels = seq(0, 180, length.out = 7))+
	ylab('Number of ants in arena')


filterdt <- dt[id_out != '' & `T` <= max(food[['t']]), .(id = parse_ids(id_out)), by = 'Frame']
M <- dcast(data = filterdt, Frame ~ id)
.tmp <- rbindlist(lapply(2:ncol(M), function(i){
	idx <- which(!is.na(M[[i]]))
	data.frame(id = unique(M[[i]][idx]),Frame = M[[1]][idx], d=c(1, diff(idx)))
	}))
.tmp[d > 1, 'd'] <- 0
.tmp[, interval := cumsum(c(TRUE, diff(d) != 0)), by = id]
.tmp[, interval := ifelse(interval %% 2 == 0, interval + 1, interval), by = id]
.tmp_result <- .tmp[, .(start_frame = min(Frame), end_frame = max(Frame)), by = .(id, interval)]
f <- range(food[['t']])*2
.tp1 <- .tmp_result[start_frame <= f[1]][, end_frame := ifelse(end_frame < f[1], end_frame, f[1])]
tp1_value <- sum(.tp1[, end_frame] - .tp1[, start_frame])/2
.tp2 <- .tmp_result[start_frame >= f[1] & start_frame <= f[2]][, end_frame := ifelse(end_frame < f[2], end_frame, f[2])]
tp2_value <- sum(.tp2[, end_frame] - .tp2[, start_frame])/2


dt_ids <- dt[id_out != '' & `T` <= max(food[['t']]), .(id = parse_ids(id_out), Frame = Frame)]
dt_tp1 <- sum(filterdt[`T` <= min(food[['t']]), .(time = get_time(`T`)), by = 'id'][['time']])
dt_tp2 <- sum(filterdt[, .(time = get_time(`T`)), by = 'id'][['time']])



files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('.*uniform_\\d{1,2}', files))))

t0 <- Sys.time()
result_list <- lapply(seq_along(files), function(x){
	do.call('gc', args = list(verbose = FALSE))
	ind_efficiency(path, files[x])
})
Sys.time()- t0
names(result_list) <- files
L <- result_list[lapply(result_list,length) > 1]
# L <- lapply(L, function(i) data.table(t(i)))
result_dt <- rbindlist(L, idcol = TRUE)

t0 <- Sys.time()
result_list_v2 <- lapply(seq_along(files), function(x){
	do.call('gc', args = list(verbose = FALSE))
	ind_efficiency_v2(path, files[x])
})
Sys.time()- t0
names(result_list_v2) <- files
Lv2 <- result_list_v2[lapply(result_list_v2,length) > 1]
result_dt_v2 <- rbindlist(Lv2, idcol = TRUE)

# ggplot(data = melt(result_dt, id.vars = '.id'), aes('', 1/value))+
# 	geom_boxplot()+ 
# 	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')),
# 		   scales = 'free_y')+
# 	xlab('') + ylab(TeX('$Efficiency (s^{-1})$'))




###########
# filterdt <- dt[id_out != '' & `T` <= max(food[['t']]), .(id = parse_ids(id_out)), by = 'T']
# dt_tp1 <- sum(filterdt[`T` <= min(food[['t']]), .(time = get_time(`T`)), by = 'id'][['time']])
# dt_tp2 <- sum(filterdt[, .(time = get_time(`T`)), by = 'id'][['time']])


det_tps <- lapply(seq_along(det), function(i){
	
	food <- range(do.call('rbind', det[[i]]@food)[['t']])
	dt <- setDT(det[[i]]@data)
	filterdt <- dt[Frame <= food[2], .(id = N_ind), by = 'Frame']
	tp1 <- sum(filterdt[Frame <= food[1], .(time = get_time(Frame)), by = 'id'][['time']])/2
	tp2 <- sum(filterdt[, .(time = get_time(Frame)), by = 'id'][['time']])/2
	data.table(condition = 'exp', tp1 = tp1, tp2 = tp2)
})

sto_tps <- lapply(seq_along(sto), function(i){
	food <- range(do.call('rbind', sto[[i]]@food)[['time']][['t']])
	dt <- setDT(sto[[i]]@data)
	filterdt <- dt[Frame <= food[2], .(id = N_ind), by = 'Frame']
	tp1 <- sum(filterdt[Frame <= food[1], .(time = get_time(Frame)), by = 'id'][['time']])/2
	tp2 <- sum(filterdt[, .(time = get_time(Frame)), by = 'id'][['time']])/2
	data.table(condition = 'exp', tp1 = tp1, tp2 = tp2)
})

tps <- rbindlist(l = list(rbindlist(det_tps),rbindlist(sto_tps)))
colnames(result_dt)[1] <- 'condition'
result_dt[['condition']] <- 'sim'



ggplot(data = melt(rbindlist(list(result_dt, tps)), id.vars = 'condition'), 
       aes(condition, 1/value))+
	geom_boxplot()+ 
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')),
		   scales = 'free_y')+
	xlab('') + ylab(TeX('$Efficiency (s^{-1})$'))




det_tps_v2 <- lapply(seq_along(det), function(i){
	
	food <- range(do.call('rbind', det[[i]]@food)[['t']])
	dt <- setDT(det[[i]]@data)
	filterdt <- dt[Frame <= food[2], .(id = N_ind), by = 'Frame']
	tp1 <- sum(filterdt[Frame <= food[1], .(time = get_time(Frame)), by = 'id'][['time']])/2
	tp2 <- sum(filterdt[Frame > food[1], .(time = get_time(Frame)), by = 'id'][['time']])/2
	data.table(condition = 'exp', tp1 = tp1, tp2 = tp2)
})

sto_tps_v2 <- lapply(seq_along(sto), function(i){
	food <- range(do.call('rbind', sto[[i]]@food)[['t']])
	dt <- setDT(sto[[i]]@data)
	filterdt <- dt[Frame <= food[2], .(id = N_ind), by = 'Frame']
	tp1 <- sum(filterdt[Frame <= food[1], .(time = get_time(Frame)), by = 'id'][['time']])/2
	tp2 <- sum(filterdt[Frame > food[1], .(time = get_time(Frame)), by = 'id'][['time']])/2
	data.table(condition = 'exp', tp1 = tp1, tp2 = tp2)
})

tps_v2 <- rbindlist(l = list(rbindlist(det_tps_v2),rbindlist(sto_tps_v2)))
colnames(result_dt_v2)[1] <- 'condition'
result_dt_v2[['condition']] <- 'sim'


ggplot(data = melt(rbindlist(list(result_dt_v2, tps_v2)), id.vars = 'condition'), 
       aes(condition, 1/value))+
	geom_boxplot()+ 
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')),
		   scales = 'free_y')+
	xlab('') + ylab(TeX('$Efficiency (s^{-1})$'))
