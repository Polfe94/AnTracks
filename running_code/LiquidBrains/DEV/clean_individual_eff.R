source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
library(arrow)
library(latex2exp)
load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/sto.RData')


path <- '/home/polfer/research/gits/AutomatAnts/results/2024/uniform_rec/'
files <- list.files(path)
files <- unique(unlist(regmatches(files, gregexpr('uniform_\\d{1,2}', files))))

path2 <- '/home/polfer/research/gits/AutomatAnts/results/2024/uniform_noRec/'
files2 <- list.files(path2)
files2 <- unique(unlist(regmatches(files2, gregexpr('uniform_\\d{1,2}', files2))))

parse_ids <- function(x){
	as.numeric(strsplit(x, ',')[[1]])
}


get_time <- function(x, threshold = 10){
	if(is.unsorted(x)){
		x <- sort(x)
	}
	y <- diff(x)
	y <- y[y <= threshold]
	sum(y)
}


eff <- function(path, filename, minpieces = 0){
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
result_list <- lapply(seq_along(files), function(x){
	do.call('gc', args = list(verbose = FALSE))
	eff(path, files[x])
})
Sys.time()- t0
names(result_list) <- files
L <- result_list[lapply(result_list,length) > 1]
result_dt <- rbindlist(L, idcol = TRUE)

t0 <- Sys.time()
result_list2 <- lapply(seq_along(files2), function(x){
	do.call('gc', args = list(verbose = FALSE))
	eff(path2, files2[x])
})
Sys.time()- t0
names(result_list2) <- files2
L2 <- result_list2[lapply(result_list2,length) > 1]
result_dt2 <- rbindlist(L2, idcol = TRUE)



det_tps <- lapply(det, get_eff)
sto_tps <- lapply(sto, get_eff)
tps <- cbind(condition = 'exp', rbindlist(l = list(rbindlist(det_tps),rbindlist(sto_tps))))

colnames(result_dt)[1] <- 'condition'
result_dt[['condition']] <- 'sim_wRec'
colnames(result_dt2)[1] <- 'condition'
result_dt2[['condition']] <- 'sim_woutRec'



ggplot(data = melt(rbindlist(list(result_dt, result_dt2, tps)), id.vars = 'condition'), 
       aes(condition, 1/value, fill = condition))+
	geom_boxplot(alpha = 0.6)+ 
	facet_wrap(~ factor(variable, labels = c('Exploration', 'Exploitation')),
		   scales = 'free_y')+
	scale_x_discrete('', labels = c('Experiment', 'Without recruitment', 'With recruitment')) + 
	ylab(TeX('$Efficiency (s^{-1})$'))+
	scale_fill_manual('', values = c('mediumpurple', 'gold3', 'deepskyblue4'),
			  labels = c('Experiment', 'Without recruitment', 'With recruitment'))+
	theme(legend.key.size = unit(30, 'pt'))
