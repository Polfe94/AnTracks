source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
library(arrow)
library(latex2exp)
load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/sto.RData')

path <- '/home/polfer/research/gits/AutomatAnts/results/rec_noRec_nf/rec/uniform/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('.*uniform_\\d{1,2}', files))))

parse_ids <- function(x){
	as.numeric(strsplit(x, ',')[[1]])
}


# get_time <- function(x, threshold = 10){
# 	if(is.unsorted(x)){
# 		x <- sort(x)
# 	}
# 	y <- diff(x)
# 	y <- y[y <= threshold]
# 	sum(y)
# }

ind_efficiency_v2 <- function(path, filename){

	food <- data.table(read_parquet(paste0(path, filename ,'_food.parquet')))
	food <- unlist(food[!is.na(t), .(t = range(t))], use.names = FALSE)
	dt <- data.table(read_parquet(paste0(path, filename, '_data.parquet')))
	dt <- dt[id_out != '' & `T` <= food[2], .(id = parse_ids(id_out)), by = 'T']
	tp1 <- sum(dt[`T` < food[1], .(time = get_time(`T`)), by = 'id'][['time']])
	tp2 <- sum(dt[`T` > food[1], .(time = get_time(`T`)), by = 'id'][['time']])
	data.table(tp1 = tp1, tp2 = tp2)
	
}

eff <- function(path, filename, minpieces = 12){
	food <- data.table(read_parquet(paste0(path, filename ,'_food.parquet')))
	food <- unlist(food[!is.na(t), .(t = range(t))], use.names = FALSE) *2
	if(length(food)){
		dt <- data.table(read_parquet(paste0(path, filename, '_data.parquet')))
		
		
		filterdt <- dt[Frame <= food[2], .(id = parse_ids(id_out)), by = 'Frame']
		M <- dcast(data = filterdt, Frame ~ id, value.var = 'id')
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

t0 <- Sys.time()
eff(path, 'uniform_0')
Sys.time()-t0

t0 <- Sys.time()
ind_efficiency_v2(path, 'uniform_0')
Sys.time()-t0



t0 <- Sys.time()
eff(path, 'uniform_1')
Sys.time()-t0

t0 <- Sys.time()
ind_efficiency_v2(path, 'uniform_1')
Sys.time()-t0








t0 <- Sys.time()
result_list_v2 <- lapply(seq_along(files), function(x){
	do.call('gc', args = list(verbose = FALSE))
	ind_efficiency_v2(path, files[x])
})
Sys.time()- t0
names(result_list_v2) <- files
Lv2 <- result_list_v2[lapply(result_list_v2,length) > 1]
result_dt_v2 <- rbindlist(Lv2, idcol = TRUE)





t0 <- Sys.time()
result_list <- lapply(seq_along(files), function(x){
	do.call('gc', args = list(verbose = FALSE))
	eff(path, files[x])
})
Sys.time()- t0
names(result_list) <- files
L <- result_list[lapply(result_list,length) > 1]
result_dt <- rbindlist(L, idcol = TRUE)





det_tps <- lapply(det, get_eff)
sto_tps <- lapply(sto, get_eff)
tps <- cbind(condition = 'exp', rbindlist(l = list(rbindlist(det_tps),rbindlist(sto_tps))))

colnames(result_dt)[1] <- 'condition'
result_dt[['condition']] <- 'sim_1'
# colnames(result_dt_v2)[1] <- 'condition'
# result_dt_v2[['condition']] <- 'sim_2'

ggplot(data = melt(rbindlist(list(result_dt, tps)), id.vars = 'condition'), 
       aes(condition, 1/value))+
	geom_boxplot()+ 
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')),
		   scales = 'free_y')+
	xlab('') + ylab(TeX('$Efficiency (s^{-1})$'))
# ggplot(data = melt(rbindlist(list(result_dt,result_dt_v2, tps)), id.vars = 'condition'), 
#        aes(condition, 1/value))+
# 	geom_boxplot()+ 
# 	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')),
# 		   scales = 'free_y')+
# 	xlab('') + ylab(TeX('$Efficiency (s^{-1})$'))
