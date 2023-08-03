source('~/research/2022/ANTS/AnTracks/src/Experiment.R')
load('~/research/2022/ANTS/AnTracks/data/nf.RData')
load('~/research/2022/ANTS/AnTracks/data/det.RData')
load('~/research/2022/ANTS/AnTracks/data/sto.RData')

exp <- nf[[1]]
df <- exp@data[, c('N_ind', 'node')]
dt <- data.table(df)
Ns <- dt[, .N, N_ind]
inds <- Ns$N_ind[which(Ns$N > 2000)]
nds <- sapply(2:length(exp@data$node), function(i){
	if(exp@data$node[i] != exp@data$node[i-1]){
		return(c(exp@data$node[i], exp@data$N_ind[i]))
	}
})


nds[[1]] <- c(exp@data$node[1], exp@data$N_ind[1])
nds <- nds[lapply(nds, length) > 0]

######## +++++++++++++++++++++++++++++++++++++++++++++++++++++ ########
## per individual
df_ind_122 <- exp@data[exp@data$N_ind == 122, 'node']
nds_ind_122 <- sapply(2:length(df_ind_122), function(i){
	if(df_ind_122[i] != df_ind_122[i-1]){
		return(df_ind_122[i])
	}
})
nds_ind_122[[2]] <- nds_ind_122[[1]]
nds_ind_122[[1]] <- 689
nds_ind_122 <- nds_ind_122[lapply(nds_ind_122, length) > 0]
dirs_122 <- vapply(3:length(nds_ind_122), function(i){
	x <- list(hex[hex$node == nds_ind_122[[i-2]], c('x', 'y')],
	     hex[hex$node == nds_ind_122[[i-1]], c('x', 'y')],
	     hex[hex$node == nds_ind_122[[i]], c('x', 'y')])
	get_direction(x)
}, numeric(1))
result <- dirs_122


M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
for(ii in 2:length(result)){
	M[as.character(result[ii-1]), as.character(result[ii])] <-
		M[as.character(result[ii-1]), as.character(result[ii])]+1
}


M / sum(M)
apply(M / sum(M), 2, function(i) i / sum(i))



#################################################################
## ONLY INDIVIDUALS WITH MORE THAN 2000 POSITIONS
## per individual
result_list <- vector('list', length(inds))
M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
for(j in seq_along(inds)){
	df_ind <- exp@data[exp@data$N_ind == inds[j], 'node']
	nds_ind <- sapply(2:length(df_ind), function(i){
		if(df_ind[i] != df_ind[i-1]){
			return(df_ind[i])
		}
	})
	nds_ind <- append(df_ind[1], nds_ind)
	nds_ind <- nds_ind[lapply(nds_ind, length) > 0]
	dirs <- vapply(3:length(nds_ind), function(i){
		x <- list(hex[hex$node == nds_ind[[i-2]], c('x', 'y')],
			  hex[hex$node == nds_ind[[i-1]], c('x', 'y')],
			  hex[hex$node == nds_ind[[i]], c('x', 'y')])
		get_direction(x)
	}, numeric(1))
	result <- na.omit(dirs)
	
	M_ind <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
	for(ii in 2:length(result)){
		M[as.character(result[ii-1]), as.character(result[ii])] <-
			M[as.character(result[ii-1]), as.character(result[ii])]+1
		M_ind[as.character(result[ii-1]), as.character(result[ii])] <-
			M_ind[as.character(result[ii-1]), as.character(result[ii])]+1
	}
	
	
	M_ind / sum(M_ind)
	result_list[[j]] <- apply(M_ind / sum(M_ind), 2, function(i) i / sum(i))
	
}
apply(M / sum(M), 2, function(i) i / sum(i))


#################################################################
## ALL INDIVIDUALS

inds <- Ns$N_ind[Ns$N > 20]
result_list <- vector('list', length(inds))
M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
for(j in seq_along(inds)){
	df_ind <- exp@data[exp@data$N_ind == inds[j], 'node']
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
			}
			
			
			M_ind / sum(M_ind)
			result_list[[j]] <- apply(M_ind / sum(M_ind), 2, function(i) i / sum(i))
		}
	}
}
M
apply(M / sum(M), 2, function(i) i / sum(i))



#################################################################
## ALL INDIVIDUALS, EXP AFTER STOCHASTICS
exp <- nf[[7]]
df <- exp@data[, c('N_ind', 'node')]
dt <- data.table(df)
Ns <- dt[, .N, N_ind]
inds <- Ns$N_ind[Ns$N > 20]
result_list <- vector('list', length(inds))
M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
for(j in seq_along(inds)){
	df_ind <- exp@data[exp@data$N_ind == inds[j], 'node']
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
			}
			
			
			M_ind / sum(M_ind)
			result_list[[j]] <- apply(M_ind / sum(M_ind), 2, function(i) i / sum(i))
		}
	}
}
M
apply(M / sum(M), 2, function(i) i / sum(i))


#################################################################
## ALL NF EXPERIMENTS
M_nf <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
for(e in seq_along(nf)){
	exp <- nf[[e]]
	df <- exp@data[, c('N_ind', 'node')]
	dt <- data.table(df)
	Ns <- dt[, .N, N_ind]
	inds <- Ns$N_ind[Ns$N > 20]
	result_list <- vector('list', length(inds))
	M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
	for(j in seq_along(inds)){
		df_ind <- exp@data[exp@data$N_ind == inds[j], 'node']
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

M_det <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
for(e in seq_along(det)){
	exp <- det[[e]]
	df <- exp@data[, c('N_ind', 'node')]
	dt <- data.table(df)
	Ns <- dt[, .N, N_ind]
	inds <- Ns$N_ind[Ns$N > 20]
	result_list <- vector('list', length(inds))
	M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
	for(j in seq_along(inds)){
		df_ind <- exp@data[exp@data$N_ind == inds[j], 'node']
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

M_sto <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
for(e in seq_along(sto)){
	exp <- sto[[e]]
	df <- exp@data[, c('N_ind', 'node')]
	dt <- data.table(df)
	Ns <- dt[, .N, N_ind]
	inds <- Ns$N_ind[Ns$N > 20]
	result_list <- vector('list', length(inds))
	M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
	for(j in seq_along(inds)){
		df_ind <- exp@data[exp@data$N_ind == inds[j], 'node']
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
