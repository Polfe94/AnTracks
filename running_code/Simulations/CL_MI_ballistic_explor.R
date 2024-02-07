library(data.table)
library(arrow)
library(infotheo)
library(parallel)

path <- '/home/usuaris/pol.fernandez/research/AutomatAnts/results/2024/updated_params/nest_avoid_explor/'
f <- list.files(path)
files <- f[grepl('data', f)]

nCores <- 20L # <--- Change suitably !

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
			result[[i]] <- mutinformation(m)
			
		} else {
			result[[i]] <- 0
		}
		do.call('gc', args = list(verbose = FALSE))
		print(paste0('Iter ', i, ' complete; MI = ', round(result[[i]], 4), ' at time = ', sq[i]))
	}
	result
}


sim_mi <- mclapply(X = seq_along(files),
		   FUN = function(i){
		   	x <- data.table(read_parquet(paste0(path, files[i])))
		   	MI(x, tw = 2400, shift = 300)
		   },
		   mc.cores = nCores)

save(sim_mi, file = '/home/usuaris/pol.fernandez/research/AnTracks/results/MI_ballistic_explor.RData')