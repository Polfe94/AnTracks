###################################################################
### SUPPORTING FUNCTIONS TO READ AND MANIPULATE SIMULATION DATA ###
###################################################################

## function to load output from (python) model simulations
read_ModelOutput <- function(path, verbose = TRUE){
	
	files <- list.files(path)
	files <- files[grepl('.json', files)]
	
	sims <- vector('list', length(files))
	
	for(i in seq_along(sims)){
		dir <- paste0(path, files[i])
		sims[[i]] <- data.table::setDT(RJSONIO::fromJSON(dir, nullValue = NA)[1:5])
		if(verbose){
			cat('\r Progress =', round(100*i/length(sims)), '%')
			flush.console()
		}
	}
	sims
}

## function to load output from (python) gain model simulations
read_gainsOutput <- function(path, verbose = TRUE){
	
	files <- list.files(path)
	files <- files[grepl('gain', files)]
	
	sims <- vector('list', length(files))
	
	for(i in seq_along(sims)){
		dir <- paste0(path, files[i])
		sims[[i]] <- data.table::setDT(RJSONIO::fromJSON(dir, nullValue = NA))
		set(sims[[i]], j = 'gain', value = strsplit(gsub('.json','', files[i]), '_')[[1]][1])
		if(verbose){
			cat('\r Progress =', round(100*i/length(sims)), '%')
			flush.console()
		}
	}
	rbindlist(sims, idcol = TRUE)
}

align_sequences <- function(df, sim_type = 'gain'){
	dt <- setDT(df)
	set(dt, j = 'Frame', value = cut(dt[['T']], seq(0, 10800, 0.5), labels = F))
	dt <- dt[!is.na(dt[['Frame']])]
	dt_clean <- dt[, .(N = mean(N), Si = mean(SiOut),
			   variable = unique(get(sim_type))), by = c('.id', 'Frame')]
	colnames(dt_clean)[colnames(dt_clean) == 'variable'] <- sim_type
	dt_clean
}