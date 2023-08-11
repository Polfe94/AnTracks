library(data.table)

path <- '~/research/AutomatAnts/results/parameter_space/'
	
files <- list.files(path)
files <- files[grepl('.json', files)]

df <- data.frame(alpha = rep(0, length(files)), 
		 beta = 0, Si = 0, vSi = 0)

for(i in seq_along(files)){
	dir <- paste0(path, files[i])
	x <- data.table::setDT(RJSONIO::fromJSON(dir, nullValue = NA))[['SiOut']]
	df[i, c('Si', 'vSi')] <- c(mean(x, na.rm = T), var(x, na.rm = T))
	
	n <- strsplit(files[i], '_')[[1]][1:2]
	a <- strsplit(n, '=')
	
	for(j in seq_along(a)){
		df[i, a[[j]][1]] <- as.numeric(a[[j]][2])
	}
}
dt <- data.table::setDT(df)
dt_avg <- dt[, lapply(.SD, mean), .SDcols = c('Si', 'vSi'), by = c('alpha', 'beta')]

save(dt, file = paste0(path, 'phase_space_data.RData'))
save(dt_avg, file = paste0(path, 'phase_space_summary.RData'))
