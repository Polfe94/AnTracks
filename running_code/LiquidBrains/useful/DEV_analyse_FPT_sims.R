library(parallel)
library(data.table)
library(arrow)

N <- 10000

path <- '/home/usuaris/pol.fernandez/research/AutomatAnts/results/2024/movement_simulations/'
files <- list.files(path)
files <- files[grepl('data', files)]

mov_rho_R <- mclapply(files, function(i){
	data <- data.table(read_parquet(paste0(path, i))) 
	rho_R <- as.numeric(strsplit(i, '_')[[1]][c(2, 4)])
	LR <- round(N * rho_R[1])
	SR <- N - LR
	
	mov <- data[, .N, by = 'id']
	mov[['rho']] <- rho_R[1]
	mov[['R']] <- rho_R[2]
	mov[['scout']] <- c(rep('LR', LR), rep('SR', SR))
	
	mov
}, mc.cores = 16L)

save(mov_rho_R, file = '~/research/AnTracks/results/mov_rho_R.RData')


