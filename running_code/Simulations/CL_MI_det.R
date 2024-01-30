library(data.table)
library(infotheo)
library(parallel)

load('/home/usuaris/pol.fernandez/research/AnTracks/data/det.RData')
nCores <- 10L # <--- Change suitably !

MI <- function(x, tw = 2400, shift = 300){
	x <- x[, .(N = .N), by = c('Frame', 'node')]
	t0 <- tw %/% 2
	mint <- min(x[N > 0, Frame])
	sq <- seq(t0, max(x[['Frame']]) - tw, shift)
	result <- vector('list', length = length(sq))
	names(result) <- sq
	for(i in seq_along(sq)){
		upper_lim <- sq[i]+ t0 -1
		if(upper_lim > mint){
			current_sequence <- seq(sq[i] - t0, sq[i] + t0 -1, 1)
			dt <- x[Frame %in% current_sequence]
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
	}
	result
}

det_mi <- mclapply(X = seq_along(det),
			       FUN = function(i){
			       	x <- data.table(det[[i]]@data)[, .(N = .N, node = node), by = 'Frame']
			       	MI(x, tw = 2400, shift = 300)
			       },
			       mc.cores = nCores)

save(det_mi, file = '/home/usuaris/pol.fernandez/research/AnTracks/results/MI_det_full.RData')