path_static <- '~/research/AutomatAnts/results/2024/hetero_model/rho/static/'
path_shutdown <- '~/research/AutomatAnts/results/2024/hetero_model/rho/shutdown/'
f <- list.files(path_static)
files <- unique(unlist(regmatches(f, gregexpr('rho_\\d{1}\\.\\d{1,3}_\\d{1,2}', f))))

library(arrow)
library(data.table)
library(parallel)

nCores <- 20L

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


# collective efficiency
collective_static <- mclapply(seq_along(files), function(i){
	do.call('gc', args = list(verbose = FALSE))
	df <- eff(path_static, files[i])
	df[['rho']] <- round(as.numeric(strsplit(files[i], '_')[[1]][2]), 2)
	df
}, mc.cores = nCores)

collective_shutdown <- mclapply(seq_along(files), function(i){
	do.call('gc', args = list(verbose = FALSE))
	df <- eff(path_shutdown, files[i])
	df[['rho']] <- round(as.numeric(strsplit(files[i], '_')[[1]][2]), 2)
	df
}, mc.cores = nCores)

save(collective_shutdown, file = '/home/usuaris/pol.fernandez/research/AnTracks/results/collective_shutdown.RData')
save(collective_static, file = '/home/usuaris/pol.fernandez/research/AnTracks/results/collective_static.RData')