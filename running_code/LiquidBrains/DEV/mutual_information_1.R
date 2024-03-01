library(data.table)
library(arrow)
library(infotheo)

path <- '/home/polfer/research/gits/AutomatAnts/results/2024/scenarios/g0.75_J_0.01/'
x <- data.table(read_parquet(paste0(path, 'g0.75_J_0.01_0_data.parquet')))

x[['Frame']] <- round(x[['T']])

parse_nodes <- function(nodes){
	unlist(strsplit(nodes, ';'))
}

parse_ids <- function(ids){
	as.integer(unlist(strsplit(ids, ',')))
}

MI <- function(x, tw = 2400, shift = 300){
	t0 <- tw %/% 2
	mint <- min(x[N > 0, Frame])
	sq <- seq(t0, max(x[['Frame']]) - tw, shift)
	result <- vector('list', length = length(sq))
	names(result) <- sq
	for(i in seq_along(sq)){
		upper_lim <- sq[i]+ t0 -1
		if(upper_lim > mint){
			current_sequence <- seq(sq[i] - t0, sq[i] + t0 -1, 1)
			dt <- x[Frame %in% current_sequence, .(node = parse_nodes(pos)), 
				by = 'Frame'][, .(N = .N), by = c('Frame', 'node')]
			M <- dt[CJ(Frame = current_sequence, node = node, unique = T),
				on=. (Frame, node)][is.na(N) | N < 2, N:=-1][N > 0, N:=1]
			m <- dcast(M, Frame ~ node, value.var = 'N')
			M <- NULL
			m$Frame <- NULL
			result[[i]] <- mean(mutinformation(m))

		} else {
			result[[i]] <- 0
		}
		do.call('gc', args = list(verbose = FALSE))
	}
	result
}


madak <- CJ(colnames(m), colnames(m), unique = TRUE)


t0 <- Sys.time()
mean(unlist(parallel::mclapply(1:nrow(madak), 
		   function(i){
		   	mutinformation(m[[madak[i, V1]]], m[[madak[i, V2]]])
		   }, mc.cores = 8L)))
Sys.time()-t0

'Time difference of 8.959998 secs'


tMI <- function(x, tw = 2400, shift = 300){
	t0 <- tw %/% 2
	mint <- min(x[N > 0, Frame])
	sq <- seq(t0, max(x[['Frame']]) - tw, shift)
	result <- parallel::mclapply(seq_along(sq), function(i){
		upper_lim <- sq[i]+ t0 -1
		if(upper_lim > mint){
			current_sequence <- seq(sq[i] - t0, sq[i] + t0 -1, 1)
			dt <- x[Frame %in% current_sequence, .(node = parse_nodes(pos)), 
				by = 'Frame'][, .(N = .N), by = c('Frame', 'node')]
			M <- dt[CJ(Frame = current_sequence, node = node, unique = T),
				on=. (Frame, node)][is.na(N), N:=-1][N > 0, N:=1]
			m <- dcast(M, Frame ~ node, value.var = 'N')
			m$Frame <- NULL
			mean(mutinformation(m))
		} else {
			NA
		}
	}, mc.cores = 8L)
	names(result) <- sq
	result[lapply(result, is.na) == FALSE]
}