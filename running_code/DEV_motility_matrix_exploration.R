library(data.table)

source('~/research/gits/AnTracks/src/Experiment.R')
load('~/research/gits/AnTracks/data/nf.RData')
load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/sto.RData')

#################################################################
#################################################################

tmax <- 20*120
min_pos <- 15
## ALL NF EXPERIMENTS
M_nf <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
for(e in seq_along(nf)){
	exp <- nf[[e]]
	dt <- setDT(exp@data)[Frame <= (min(Frame) + tmax), c('N_ind', 'node')]
	Ns <- dt[, .N, N_ind]
	inds <- Ns[N >= min_pos, N_ind]
	result_list <- vector('list', length(inds))
	M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
	for(j in seq_along(inds)){
		df_ind <- exp@data[N_ind == inds[j], node]
		nds_ind <- sapply(2:length(df_ind), function(i){
			if(df_ind[i] != df_ind[i-1]){
				return(df_ind[i])
			}
		})
		nds_ind <- append(df_ind[1], nds_ind)
		nds_ind <- nds_ind[lapply(nds_ind, length) > 0]
		if(length(nds_ind) > 2){
			dirs <- vapply(3:length(nds_ind), function(i){
				x <- list(hex[hex$node == nds_ind[[i-2]], c('x', 'y')],
					  hex[hex$node == nds_ind[[i-1]], c('x', 'y')],
					  hex[hex$node == nds_ind[[i]], c('x', 'y')])
				get_direction(x)
			}, numeric(1))
			result <- na.omit(dirs)
			if(length(result) > 1){
				M_ind <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
				for(ii in 2:length(result)){
					M[as.character(result[ii-1]), as.character(result[ii])] <-
						M[as.character(result[ii-1]), as.character(result[ii])]+1
					M_ind[as.character(result[ii-1]), as.character(result[ii])] <-
						M_ind[as.character(result[ii-1]), as.character(result[ii])]+1
					M_nf[as.character(result[ii-1]), as.character(result[ii])] <-
						M_nf[as.character(result[ii-1]), as.character(result[ii])]+1
				}
				
				
				M_ind / sum(M_ind)
				result_list[[j]] <- apply(M_ind / sum(M_ind), 2, function(i) i / sum(i))
			}
		}
	}
	M
	apply(M / sum(M), 2, function(i) i / sum(i))
}
apply(M_nf / sum(M_nf), 2, function(i) i / sum(i))


#################################################################
## ALL DETERMINIST EXPS


min_pos <- 15

M_det <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
for(e in seq_along(det)){
	exp <- det[[e]]
	tmax <- min(rbindlist(exp@food)[['t']])
	dt <- setDT(exp@data)[Frame <= tmax, c('N_ind', 'node')]
	Ns <- dt[, .N, N_ind]
	inds <- Ns[N >= min_pos, N_ind]
	result_list <- vector('list', length(inds))
	M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
	for(j in seq_along(inds)){
		df_ind <- exp@data[N_ind == inds[j], node]
		nds_ind <- sapply(2:length(df_ind), function(i){
			if(df_ind[i] != df_ind[i-1]){
				return(df_ind[i])
			}
		})
		nds_ind <- append(df_ind[1], nds_ind)
		nds_ind <- nds_ind[lapply(nds_ind, length) > 0]
		if(length(nds_ind) > 2){
			dirs <- vapply(3:length(nds_ind), function(i){
				x <- list(hex[hex$node == nds_ind[[i-2]], c('x', 'y')],
					  hex[hex$node == nds_ind[[i-1]], c('x', 'y')],
					  hex[hex$node == nds_ind[[i]], c('x', 'y')])
				get_direction(x)
			}, numeric(1))
			result <- na.omit(dirs)
			if(length(result) > 1){
				M_ind <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
				for(ii in 2:length(result)){
					M[as.character(result[ii-1]), as.character(result[ii])] <-
						M[as.character(result[ii-1]), as.character(result[ii])]+1
					M_ind[as.character(result[ii-1]), as.character(result[ii])] <-
						M_ind[as.character(result[ii-1]), as.character(result[ii])]+1
					M_det[as.character(result[ii-1]), as.character(result[ii])] <-
						M_det[as.character(result[ii-1]), as.character(result[ii])]+1
				}
				
				
				M_ind / sum(M_ind)
				result_list[[j]] <- apply(M_ind / sum(M_ind), 2, function(i) i / sum(i))
			}
		}
	}
	M
	apply(M / sum(M), 2, function(i) i / sum(i))
}
apply(M_det / sum(M_det), 2, function(i) i / sum(i))


#################################################################
## ALL STOCHASTIC EXPS
min_pos <- 15

M_sto <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
for(e in seq_along(sto)){
	exp <- sto[[e]]
	tmax <- min(rbindlist(exp@food)[['t']])
	dt <- setDT(exp@data)[Frame <= tmax, c('N_ind', 'node')]
	Ns <- dt[, .N, N_ind]
	inds <- Ns[N >= min_pos, N_ind]
	result_list <- vector('list', length(inds))
	M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
	for(j in seq_along(inds)){
		df_ind <- exp@data[N_ind == inds[j], node]
		nds_ind <- sapply(2:length(df_ind), function(i){
			if(df_ind[i] != df_ind[i-1]){
				return(df_ind[i])
			}
		})
		nds_ind <- append(df_ind[1], nds_ind)
		nds_ind <- nds_ind[lapply(nds_ind, length) > 0]
		if(length(nds_ind) > 2){
			dirs <- vapply(3:length(nds_ind), function(i){
				x <- list(hex[hex$node == nds_ind[[i-2]], c('x', 'y')],
					  hex[hex$node == nds_ind[[i-1]], c('x', 'y')],
					  hex[hex$node == nds_ind[[i]], c('x', 'y')])
				get_direction(x)
			}, numeric(1))
			result <- na.omit(dirs)
			if(length(result) > 1){
				M_ind <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
				for(ii in 2:length(result)){
					M[as.character(result[ii-1]), as.character(result[ii])] <-
						M[as.character(result[ii-1]), as.character(result[ii])]+1
					M_ind[as.character(result[ii-1]), as.character(result[ii])] <-
						M_ind[as.character(result[ii-1]), as.character(result[ii])]+1
					M_sto[as.character(result[ii-1]), as.character(result[ii])] <-
						M_sto[as.character(result[ii-1]), as.character(result[ii])]+1
				}
				
				
				M_ind / sum(M_ind)
				result_list[[j]] <- apply(M_ind / sum(M_ind), 2, function(i) i / sum(i))
			}
		}
	}
	M
	apply(M / sum(M), 2, function(i) i / sum(i))
}
apply(M_sto / sum(M_sto), 2, function(i) i / sum(i))
