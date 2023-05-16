source('~/research/2022/ANTS/AnTracks/src/Experiment.R')
load('~/research/2022/ANTS/AnTracks/data/det.RData')
load('~/research/2022/ANTS/AnTracks/data/sto.RData')
load('~/research/2022/ANTS/AnTracks/data/nf.RData')


unpack <- function(x){
	# -1, 0, 1
	result <- numeric(length(x)-1)
	for(i in 2:length(x)){
		if(x[i-1] == 'right'){
			if(x[i] == 'left'){
				result[i-1] <- 0
			} else if(x[i] == 'right'){
				result[i-1] <- 1
			} else if(x[i] == 'up'){
				result[i-1] <- -1
			} else {
				result[i-1] <- 1
			}
		} else if(x[i-1] == 'left'){
			if(x[i] == 'left'){
				result[i-1] <- -1
			} else if(x[i] == 'right'){
				result[i] <- 0
			} else if(x[i] == 'up'){
				result[i-1] <- 1
			} else {
				result[i-1] <- -1
			}
		} else if(x[i-1] == 'up'){
			if(x[i] == 'left'){
				result[i-1] <- -1
			} else if(x[i] == 'right'){
				result[i] <- 1
			} else {
				result[i-1] <- 0
			}
		} else {
			if(x[i] == 'left'){
				result[i-1] <- 1
			} else if(x[i] == 'right'){
				result[i] <- -1
			} else {
				result[i-1] <- 0
			}
		}
		
	}
	result
}

NF_RESULT <- lapply(nf, function(x){
	df <- unique(x@data[, c('N_ind', 'node')])

	L <- vector('list', length(unique(df$N_ind)))
	for(i in seq_along(L)){
		idx <- df_unique$N_ind == i
		l <- sum(idx)
		if(l > 1){
			sbst <- df_unique[idx, ]
			x <- 2
			while(x <= l){
				L[[i]] <- c(L[[i]], get_direction(sbst$node[x-1], sbst$node[x]))
				x <- x + 1
			}
			
		}
	}
	L <- L[lapply(L, length) > 2]
	
	result <- lapply(L, unpack)
	M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
	for(i in seq_along(result)){
		for(ii in 2:length(result[[i]])){
			M[as.character(result[[i]][ii-1]), as.character(result[[i]][ii])] <-
				M[as.character(result[[i]][ii-1]), as.character(result[[i]][ii])]+1
		}
		
	}
	M / sum(M)
})

MNF <- matrix(lapply(1:9, function(i){
	v <- vapply(NF_RESULT, function(x){
		x[i]
	}, numeric(1))
	mean(v)
}), ncol = 3, nrow = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))






DET_RESULT <- lapply(det, function(x){
	df <- unique(x@data[, c('N_ind', 'node')])
	
	L <- vector('list', length(unique(df$N_ind)))
	for(i in seq_along(L)){
		idx <- df_unique$N_ind == i
		l <- sum(idx)
		if(l > 1){
			sbst <- df_unique[idx, ]
			x <- 2
			while(x <= l){
				L[[i]] <- c(L[[i]], get_direction(sbst$node[x-1], sbst$node[x]))
				x <- x + 1
			}
			
		}
	}
	L <- L[lapply(L, length) > 2]
	
	result <- lapply(L, unpack)
	M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
	for(i in seq_along(result)){
		for(ii in 2:length(result[[i]])){
			M[as.character(result[[i]][ii-1]), as.character(result[[i]][ii])] <-
				M[as.character(result[[i]][ii-1]), as.character(result[[i]][ii])]+1
		}
		
	}
	M
})

MDET <- matrix(unlist(lapply(1:9, function(i){
	v <- vapply(DET_RESULT, function(x){
		x[i]
	}, numeric(1))
	sum(v)
})), ncol = 3, nrow = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))

MDET / sum(MDET)


STO_RESULT <- lapply(sto, function(x){
	df <- unique(x@data[, c('N_ind', 'node')])
	
	L <- vector('list', length(unique(df$N_ind)))
	for(i in seq_along(L)){
		idx <- df_unique$N_ind == i
		l <- sum(idx)
		if(l > 1){
			sbst <- df_unique[idx, ]
			x <- 2
			while(x <= l){
				L[[i]] <- c(L[[i]], get_direction(sbst$node[x-1], sbst$node[x]))
				x <- x + 1
			}
			
		}
	}
	L <- L[lapply(L, length) > 2]
	
	result <- lapply(L, unpack)
	M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
	for(i in seq_along(result)){
		for(ii in 2:length(result[[i]])){
			M[as.character(result[[i]][ii-1]), as.character(result[[i]][ii])] <-
				M[as.character(result[[i]][ii-1]), as.character(result[[i]][ii])]+1
		}
		
	}
	M
})

MSTO <- matrix(unlist(lapply(1:9, function(i){
	v <- vapply(STO_RESULT, function(x){
		x[i]
	}, numeric(1))
	sum(v)
})), ncol = 3, nrow = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))

MSTO / sum(MSTO)

# df <- det[[1]]@data
# for(i in df$N_ind){
# 	st <- df$node[]
# }
# 
# df_risas <- df[, c('N_ind', 'node')]
# df_unique <- unique(df_risas)
# 
# L <- vector('list', length(unique(df_unique$N_ind)))
# for(i in seq_along(L)){
# 	idx <- df_unique$N_ind == i
# 	l <- sum(idx)
# 	if(l > 1){
# 		sbst <- df_unique[idx, ]
# 		x <- 2
# 		while(x <= l){
# 			L[[i]] <- c(L[[i]], get_direction(sbst$node[x-1], sbst$node[x]))
# 			x <- x + 1
# 		}
# 		
# 	}
# }
# L <- L[lapply(L, length) > 5]


# result <- lapply(L, unpack)
# M <- matrix(data = 0, nrow = 3, ncol = 3, dimnames = list(c(1, 0, -1), c(1, 0, -1)))
# for(i in seq_along(result)){
# 	for(ii in 2:length(result[[i]])){
# 		M[as.character(result[[i]][ii-1]), as.character(result[[i]][ii])] <-
# 			M[as.character(result[[i]][ii-1]), as.character(result[[i]][ii])]+1
# 	}
# 	
# }
# M / sum(M)

