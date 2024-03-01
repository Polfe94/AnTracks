t0 <- Sys.time()
files_ <- f[grepl('data', f)][1]
l <- length(files)
p <- lapply(seq_along(files_), function(x){
	do.call('gc', args = list(verbose = FALSE))
	data <- data.table(read_parquet(paste0(path_shutdown_DET, files_[x])))[, .(node = parse_nodes(pos)),
								       by = 'Frame']
	rho <- round(as.numeric(strsplit(files_[x], '_')[[1]][2]), 2)
	mdata <- merge(data, hex_sim[, c('x', 'y', 'node')])[order(Frame)]
	dists <- pdist(as.matrix(hex_sim[node == '(0, 22)', c('x', 'y')]), as.matrix(mdata[, c('x', 'y')]))
	indices <- vapply(R_sims, function(i){
		which.max(dists > i)
	}, numeric(1))
	times <- mdata[indices, Frame]
	cat(paste0('Finished iteration ', formatC(x, flag = '0', digits = 3), ' from ', l, '\r'))
	data.table(R = R_sims, t = times, rho = rho)
})
cat(paste('+++++ PROGRAM FINISHED +++++ Time ellapsed: ', round(Sys.time()-t0, 2)))

message_parallel <- function(...){
	system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
}
t0 <- Sys.time()
files_ <- f[grepl('data', f)]
l <- length(files)
p <- mclapply(seq_along(files_), function(x){
	do.call('gc', args = list(verbose = FALSE))
	data <- data.table(read_parquet(paste0(path_static_DET, files_[x])))[, .(node = parse_nodes(pos)),
									       by = 'Frame']
	rho <- round(as.numeric(strsplit(files_[x], '_')[[1]][2]), 2)
	mdata <- merge(data, hex_sim[, c('x', 'y', 'node')])[order(Frame)]
	dists <- pdist(as.matrix(hex_sim[node == '(0, 22)', c('x', 'y')]), as.matrix(mdata[, c('x', 'y')]))
	indices <- vapply(R_sims, function(i){
		which.max(dists > i)
	}, numeric(1))
	times <- mdata[indices, Frame]
	message_parallel('Iteration ', formatC(x, digits = 3, flag = '0'))
	data.table(R = R_sims, t = times, rho = rho)
}, mc.cores = nCores)
cat(paste('+++++ PROGRAM FINISHED +++++ Time ellapsed: ', round(Sys.time()-t0, 2)))

