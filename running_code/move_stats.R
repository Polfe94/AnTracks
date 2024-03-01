move_stats <- function(exp, scouts = 'LR', maxt = 2400, maxd = 500, min_pos = 10){
	data <- setDT(exp@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	inds <- unique(data[as.numeric(dmatrix) > maxd, N_ind])
	if(scouts == 'SR'){
		inds <- unique(data[!N_ind %in% inds, N_ind])
	}
	set(data, j = 'd', value = dmatrix)
	data <- data[N_ind %in% inds]
	idx <- data[, .(idx = which.min(abs((d - maxd)))), by = 'N_ind'][['idx']]
	result <- as.data.frame(t(vapply(seq_along(inds), function(ii){
		sbst <- data[N_ind == inds[ii]]
		xy <- TrajFromCoords(sbst[1:idx[[ii]], c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > min_pos){
			c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
			  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
			  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE),
			  nrow(xy) / 120, nrow(sbst)/120)
		} else {rep(0, 7)}
		
	}, numeric(7))))
	colnames(result) <- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time', 'Total_time')
	result[result[['Time']] > 0, ]
}

move_stats(det[[1]], scouts = 'SR')
