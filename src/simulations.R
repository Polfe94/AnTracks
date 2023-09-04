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

#' Takes two data tables and performes the DTW function
#' @param sim a data.table comprising N time series (concatenated with an .id label)
#' @param ref a data.table comprising M time series (as columns)
#' @param align if TRUE performs `align_sequences` before running the computations
#' @param only.average if TRUE only the averaged reference sequence is taken into account
#' @return a vector of length N x M of DTW distances between the simulated and reference time series
mass_dtw <- function(sim, ref, align = FALSE, only.average = TRUE){
	if(align){
		x <- colnames(sim)
		if('gain' %in% x){
			t <- 'gain'
		} else if(grepl('parameter', x)){
			t <- 'ps'
		} else {
			stop('Could not detect the type of simulation. Try setting align to FALSE.')
		}
		sim <- align_sequences(sim, t)
	}
	result <- data.frame(.id = seq_len(max(sim[['.id']])))
	if(only.average){
		tmp <- lapply(seq_along(result[['.id']]), function(i){
			k <- dtw::dtw(sim[sim[['.id']] == i, N], ref[['avg']])
			c(d = k[['distance']], nrm_d = k[["normalizedDistance"]])
		})
		result <- setDT(cbind(result, do.call('rbind', tmp)))
	} else {
		tmp <- lapply(seq_along(result[['.id']]), function(i){
			do.call('rbind', lapply(1:ncol(ref), function(j){
				k <- dtw::dtw(sim[sim[['.id']] == i, N], ref[,..j])
				c(d = k[['distance']], nrm_d = k[["normalizedDistance"]])
			}))
		})
		result <- setDT(cbind(result, ref = colnames(ref), do.call('rbind', tmp)))
	}

	result
}

#' Takes two data tables and performes the DTW function
#' @param sim a data.table comprising N time series (concatenated with an .id label)
#' @param ref a data.table comprising M time series (as columns)
#' @param align if TRUE performs `align_sequences` before running the computations
#' @param only.average if TRUE only the averaged reference sequence is taken into account
#' @return a vector of length N x M of DTW distances between the simulated and reference time series
pdtw <- function(sim, ref, align = FALSE, only.average = TRUE){
	if(align){
		x <- colnames(sim)
		if('gain' %in% x){
			t <- 'gain'
		} else if(grepl('parameter', x)){
			t <- 'ps'
		} else {
			stop('Could not detect the type of simulation. Try setting align to FALSE.')
		}
		sim <- align_sequences(sim, t)
	}
	result <- data.frame(.id = seq_len(max(sim[['.id']])))
	if(only.average){
		tmp <- parallel::mclapply(seq_along(result[['.id']]), function(i){
			k <- dtw::dtw(sim[sim[['.id']] == i, N], ref[['avg']])
			c(d = k[['distance']], nrm_d = k[["normalizedDistance"]])
		}, mc.cores = parallel::detectCores())
		result <- setDT(cbind(result, do.call('rbind', tmp)))
	} else {
		tmp <- parallel::mclapply(seq_along(result[['.id']]), function(i){
			do.call('rbind', parallel::mclapply(1:ncol(ref), function(j){
				k <- dtw::dtw(sim[sim[['.id']] == i, N], ref[,..j])
				c(d = k[['distance']], nrm_d = k[["normalizedDistance"]])
			}, mc.cores = parallel::detectCores()))
		}, mc.cores = parallel::detectCores())
		result <- setDT(cbind(result, ref = colnames(ref), do.call('rbind', tmp)))
	}
	
	result
}