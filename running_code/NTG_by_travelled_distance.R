source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
source('~/research/gits/AnTracks/src/fit_functions.R')
load('~/research/gits/AnTracks/data/nf.RData')
load('~/research/gits/AnTracks/data/det.RData')

library(trajr)

# LRs <- lapply(seq_along(LR), function(i){
# 	sbst <- data[N_ind == LR[i]]
# 	xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
# 	Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
# 	maxL <- max(Ls)
# 	if(maxL >= 100){
# 		idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
# 	} else {
# 		idxs <- nrow(xy)
# 	}
# 	straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
# 	data.frame(l = Ls[idxs], s = straightness)
# })

maxt <- 3000
maxd <- 500
LRs <- lapply(det, function(p){
	data <- setDT(p@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	set(data, j = 'd', value = as.numeric(dmatrix))
	LR <- unique(data[d > maxd, N_ind])
	SR <- unique(data[!N_ind %in% LR, N_ind])
	data_LR <- data[N_ind %in% LR]
	data_SR <- data[N_ind %in% SR]
	rbindlist(lapply(seq_along(LR), function(i){
	sbst <- data[N_ind == LR[i]]
	xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
	Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
	maxL <- max(Ls)
	if(maxL >= 100){
		idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
	} else {
		idxs <- nrow(xy)
	}
	straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
	data.frame(l = Ls[idxs], s = straightness)

}), idcol = TRUE)})
names(LRs) <- paste0('Det_', seq_along(det))
LRs <- rbindlist(LRs, idcol = TRUE)
LRs[['tag']] <- apply(LRs[, 1:2], 1, paste, collapse = '_')

ggplot(data = LRs[,.(straightness = mean(s)), by = .(traj_length = round(l/100)*100)],
       aes(traj_length, straightness)) + geom_point()


# plot(LRs[[1]])

# SRs <- lapply(seq_along(SR), function(i){
# 	sbst <- data[N_ind == SR[i]]
# 	xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
# 	Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
# 	maxL <- max(Ls)
# 	if(maxL >= 100){
# 		idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
# 	} else {
# 		idxs <- nrow(xy)
# 	}
# 	straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
# 	
# 	data.frame(l = Ls[idxs], s = straightness)
# })
# plot(do.call('rbind', SRs))

SRs <- lapply(det, function(p){
	data <- setDT(p@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	set(data, j = 'd', value = as.numeric(dmatrix))
	SR <- unique(data[d > maxd, N_ind])
	SR <- unique(data[!N_ind %in% SR, N_ind])
	data_SR <- data[N_ind %in% SR]
	data_SR <- data[N_ind %in% SR]
	rbindlist(lapply(seq_along(SR), function(i){
		sbst <- data[N_ind == SR[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 3){
			Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
			maxL <- max(Ls)
			if(maxL >= 100){
				idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
			} else {
				idxs <- nrow(xy)
			}
			straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
			data.frame(l = Ls[idxs], s = straightness)
		}
		
	}), idcol = TRUE)})
names(SRs) <- paste0('Det_', seq_along(det))
SRs <- rbindlist(SRs, idcol = TRUE)
SRs[['tag']] <- apply(SRs[, 1:2], 1, paste, collapse = '_')

ggplot(data = SRs[,.(straightness = mean(s)), by = .(traj_length = round(l/100)*100)],
       aes(traj_length, straightness)) + geom_point()







ggplot(data = rbindlist(list(LR = LRs[,.(straightness = mean(s)), by = .(traj_length = round(l/100)*100)],
			     SR = SRs[,.(straightness = mean(s)), by = .(traj_length = round(l/100)*100)]),
			idcol = TRUE),
       aes(traj_length/ 10, straightness, color = factor(.id)), alpha = 0.5) + 
	geom_point(size = 4, shape = 21) +
	geom_smooth(formula = y ~ log(x), se = FALSE)+
	scale_y_continuous('<Net-to-Gross ratio>', limits = c(0, 1), breaks = seq(0, 1, 0.2))+
	scale_x_continuous('Travelled distance (cm)')+
	scale_color_manual('', values = c('mediumpurple', 'gold3'), 
			   labels = c('LR scouts', 'SR scouts'))+
	theme(legend.title = element_blank(), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.position = c(3/4, 3/4))







LRs_nf <- lapply(nf[-c(1,2)], function(p){
	data <- setDT(p@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	set(data, j = 'd', value = as.numeric(dmatrix))
	LR <- unique(data[d > maxd, N_ind])
	SR <- unique(data[!N_ind %in% LR, N_ind])
	data_LR <- data[N_ind %in% LR]
	data_SR <- data[N_ind %in% SR]
	rbindlist(lapply(seq_along(LR), function(i){
		sbst <- data[N_ind == LR[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
		maxL <- max(Ls)
		if(maxL >= 100){
			idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
		} else {
			idxs <- nrow(xy)
		}
		straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
		data.frame(l = Ls[idxs], s = straightness)
		
	}), idcol = TRUE)})
names(LRs_nf) <- paste0('nf_', seq_along(nf[-c(1,2)]))
LRs_nf <- rbindlist(LRs_nf, idcol = TRUE)
LRs_nf[['tag']] <- apply(LRs_nf[, 1:2], 1, paste, collapse = '_')


SRs_nf <- lapply(nf[-c(1, 2)], function(p){
	data <- setDT(p@data)[Frame <= maxt]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	set(data, j = 'd', value = as.numeric(dmatrix))
	SR <- unique(data[d > maxd, N_ind])
	SR <- unique(data[!N_ind %in% SR, N_ind])
	data_SR <- data[N_ind %in% SR]
	data_SR <- data[N_ind %in% SR]
	rbindlist(lapply(seq_along(SR), function(i){
		sbst <- data[N_ind == SR[i]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > 3){
			Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
			maxL <- max(Ls)
			if(maxL >= 100){
				idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
			} else {
				idxs <- nrow(xy)
			}
			straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
			data.frame(l = Ls[idxs], s = straightness)
		}
		
	}), idcol = TRUE)})
names(SRs_nf) <- paste0('nf_', seq_along(nf[-c(1, 2)]))
SRs_nf <- rbindlist(SRs_nf, idcol = TRUE)
SRs_nf[['tag']] <- apply(SRs_nf[, 1:2], 1, paste, collapse = '_')

data_list <- list(LR = LRs[,.(y = mean(s)), by = .(x = round(l/100)*100)],
		  SR = SRs[,.(y = mean(s)), by = .(x = round(l/100)*100)],
		  LR_nf = LRs_nf[,.(y = mean(s)), by = .(x = round(l/100)*100)],
		  SR_nf = SRs_nf[,.(y = mean(s)), by = .(x = round(l/100)*100)])

fits <- lapply(seq_along(data_list), function(i){
	fit_doub_nexp(data_list[[i]], a_0 = 0.5, k_0 = 10**-3, a2_0 = 0.1, k2_0 = 10**-4,
		      alims = c(-10, 10), klims = c(-10, 10), a2lims = c(-10, 10), k2lims = c(-10, 10))
})
names(fits) <- names(data_list)
pred_data <- lapply(seq_along(data_list), function(i){
	data.table(x = data_list[[i]][['x']], 
		   y = doub_nexp(data_list[[i]][['x']], a = fits[[i]][['a']], 
		   	      k = fits[[i]][['k']], a2 = fits[[i]][['a2']], k2 = fits[[i]][['k2']]))
})
names(pred_data) <- names(data_list)


ggplot(data = rbindlist(data_list, idcol = TRUE),
       aes(x/10, y, color = factor(.id)), alpha = 0.5) + 
	geom_point(size = 4, shape = 21) +
	geom_path(data = rbindlist(pred_data, idcol = TRUE), linewidth = 1.5)+
	scale_y_continuous('<Net-to-Gross ratio>', limits = c(0, 1), breaks = seq(0, 1, 0.2))+
	scale_x_continuous('Travelled distance (cm)')+
	scale_color_manual('', values = c('mediumpurple', 'gold3', 'brown4', 'grey'), 
			   labels = c('LR DET', 'LR NFD', 'SR DET', 'SR NFD'))+
	theme(legend.title = element_blank(), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.position = c(3/4, 3/4))

data <- LRs[,.(y = mean(s)), by = .(x = round(l/100)*100)]

fit <- fit_nexp(data, a_0 = 0.5, k_0 = 10**-3, alims = c(-10, 10), klims = c(-10, 10), plot = TRUE)
fit2 <- fit_nexpC(data, a_0 = 0.5, k_0 = 10**-3, C_0 = 0, 
		  Clims = c(-10, 10),alims = c(-10, 10), klims = c(-10, 10), plot = TRUE)
fit3 <- fit_doub_nexp(data, a_0 = 0.5, k_0 = 10**-3, a2_0 = 0.1, k2_0 = 10**-4,
		      alims = c(-10, 10), klims = c(-10, 10), a2lims = c(-10, 10), k2lims = c(-10, 10),
		      plot = TRUE)
fit4 <- fit_doub_nexpC(data, a_0 = 0.5, k_0 = 10**-3, a2_0 = 0.1, k2_0 = 10**-4,C_0 = 0, Clims = c(-10,10),ylim = c(0,1),
		      alims = c(-10, 10), klims = c(-10, 10), a2lims = c(-10, 10), k2lims = c(-10, 10),
		      plot = TRUE)
l <- nexpC(data[['x']], a = fit2[['a']], k = fit2[['k']], C = fit2[['C']])
plot(data)
lines(data[['x']],l)

l <- doub_nexp(data[['x']], a = fit3[['a']], k = fit3[['k']], a2 = fit3[['a2']], k2 = fit3[['k2']])
plot(data)
lines(data[['x']],l)

# l <- doub_nexpC(data[['x']], a = fit4[['a']], k = fit4[['k']], a2 = fit4[['a2']], k2 = fit4[['k2']], C = fit4[['C']])
# plot(data)
# lines(data[['x']],l)
# ggplot(data = rbindlist(list(LR = LRs[,.(straightness = .N), by = .(traj_length = round(l/100)*100)],
# 			     SR = SRs[,.(straightness = .N), by = .(traj_length = round(l/100)*100)],
# 			     LR_nf = LRs_nf[,.(straightness = .N), by = .(traj_length = round(l/100)*100)],
# 			     SR_nf = SRs_nf[,.(straightness = .N), by = .(traj_length = round(l/100)*100)]),
# 			idcol = TRUE),
#        aes(traj_length/ 10, straightness, color = factor(.id)), alpha = 0.5) + 
# 	geom_point(size = 4, shape = 21) +
# 	geom_smooth(formula = y ~ log(x), se = FALSE)+
# 	scale_y_continuous('<Net-to-Gross ratio>')+
# 	scale_x_continuous('Travelled distance (cm)')+
# 	scale_color_manual('', values = c('mediumpurple', 'gold3', 'brown4', 'grey'), 
# 			   labels = c('LR DET', 'LR NF', 'SR DET', 'SR NF'))+
# 	theme(legend.title = element_blank(), 
# 	      legend.background = element_rect(color = 'black', fill = NA), 
# 	      legend.position = c(3/4, 3/4))





sq <- seq(5*120, 30*120, 5*120)
LRs_nf_2d <- lapply(nf[-c(1, 2)], function(p){
	result <- lapply(sq, function(t){
		data <- setDT(p@data)[Frame <= t]
		dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
		set(data, j = 'd', value = as.numeric(dmatrix))
		LR <- unique(data[d > maxd, N_ind])
		SR <- unique(data[!N_ind %in% LR, N_ind])
		data_LR <- data[N_ind %in% LR]
		data_SR <- data[N_ind %in% SR]
		rbindlist(lapply(seq_along(LR), function(i){
			sbst <- data[N_ind == LR[i]]
			xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
			Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
			maxL <- max(Ls)
			if(maxL >= 100){
				idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
			} else {
				idxs <- nrow(xy)
			}
			straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
			diffusion <- vapply(idxs, function(l) TrajDistance(xy[seq_len(l)]), numeric(1))
			distance <- vapply(idxs, function(l) TrajEmax(xy[seq_len(l)]), numeric(1))
			data.frame(l = Ls[idxs], s = straightness, d = diffusion, emax = distance)
			
		}), idcol = TRUE)
	})
	names(result) <- sq/120
	rbindlist(result, idcol = TRUE)
	
})
names(LRs_nf_2d) <- paste0('nf_', seq_along(nf))
LRs_nf_2d <- rbindlist(LRs_nf_2d, idcol = TRUE)
LRs_nf_2d[['tag']] <- apply(LRs_nf_2d[, 1:3], 1, paste, collapse = '_')

# ggplot(data = LRs_nf_2d[,.(straightness = mean(s), t = as.numeric(.id.1),
# 			   diffusion = mean(d)), 
# 			by = .(traj_length = round(l/100)*100)],
#        aes(t, traj_length, fill = straightness, color = straightness)) + geom_tile()

ggplot(data = LRs_nf_2d[,.(straightness = mean(s), t = as.numeric(.id.1),
			   diffusion = mean(d)), 
			by = .(traj_length = round(l/100)*100)],
       aes(t, traj_length /10, fill = diffusion)) + geom_tile()+
	scale_fill_viridis_c('Diffusion',option = 'C', breaks = seq(0, 3000, 250))+
	scale_x_continuous('Time (min)', breaks = seq(5, 30, 5))+
	scale_y_continuous('Travelled distance (cm)', breaks = seq(0, 6000, 250))+
	theme(legend.position = c(0.2, 4/5), 
	      legend.background = element_rect(fill = NA, color = 'black'))

# ggplot(data = LRs_nf_2d[,.(straightness = mean(s), t = as.numeric(.id.1),
# 			   diffusion = mean(d), distance = mean(emax)), 
# 			by = .(traj_length = round(l/100)*100)],
#        aes(t, traj_length, fill = distance, color = distance)) + geom_tile()




ggplot(data = LRs_nf_2d[,.(straightness = mean(s),
			   diffusion = mean(d)), 
			by = .(traj_length = round(l/100)*100, t = .id.1)],
       aes(traj_length/ 10, diffusion, color = factor(t)), alpha = 0.5) + 
	# aes(log(traj_length/ 10), straightness, color = factor(.id)), alpha = 0.5) + 
	geom_point(size = 4, shape = 21) +
	geom_line()+
	scale_y_continuous('Diffusion', breaks = seq(0, 1500, 250))+
	scale_x_continuous('Travelled distance (cm)')+
	scale_color_manual('Time (min)', values = c('mediumpurple', 'gold3', 'brown4', 'grey', 'darkblue'),
			   labels = seq(10, 30, 5))+
	theme(legend.title = element_blank(), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.position = c(4/5, 1/2))






sq <- seq(5*120, 30*120, 5*120)
SRs_nf_2d <- lapply(nf, function(p){
	result <- lapply(sq, function(t){
		data <- setDT(p@data)[Frame <= t]
		dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
		set(data, j = 'd', value = as.numeric(dmatrix))
		SR <- unique(data[d > maxd, N_ind])
		SR <- unique(data[!N_ind %in% SR, N_ind])
		data_SR <- data[N_ind %in% SR]
		rbindlist(lapply(seq_along(SR), function(i){
			sbst <- data[N_ind == SR[i]]
			xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
			if(nrow(xy) > 3){
				Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
				maxL <- max(Ls)
				if(maxL >= 100){
					idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
				} else {
					idxs <- nrow(xy)
				}
				straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
				diffusion <- vapply(idxs, function(l) TrajDistance(xy[seq_len(l)]), numeric(1))
				distance <- vapply(idxs, function(l) TrajEmax(xy[seq_len(l)]), numeric(1))
				data.frame(l = Ls[idxs], s = straightness, d = diffusion, emax = distance)
			}

			
		}), idcol = TRUE)
	})
	names(result) <- sq/120
	rbindlist(result, idcol = TRUE)
	
})
names(SRs_nf_2d) <- paste0('nf_', seq_along(nf))
SRs_nf_2d <- rbindlist(SRs_nf_2d, idcol = TRUE)
SRs_nf_2d[['tag']] <- apply(SRs_nf_2d[, 1:3], 1, paste, collapse = '_')

# ggplot(data = SRs_nf_2d[,.(straightness = mean(s), t = as.numeric(.id.1),
# 			   diffusion = mean(d)), 
# 			by = .(traj_length = round(l/100)*100)],
#        aes(t, traj_length, fill = straightness, color = straightness)) + geom_tile()

ggplot(data = SRs_nf_2d[,.(straightness = mean(s), t = as.numeric(.id.1),
			   diffusion = mean(d)), 
			by = .(traj_length = round(l/100)*100)][t > 5],
       aes(t, traj_length /10, fill = diffusion)) + geom_tile()+
	scale_fill_viridis_c('Diffusion',option = 'C', breaks = seq(0, 300, 50))+
	scale_x_continuous('Time (min)', breaks = seq(5, 30, 5))+
	scale_y_continuous('Travelled distance (cm)', breaks = seq(0, 600, 150))+
	theme(legend.position = c(0.125, 4/5), 
	      legend.background = element_rect(fill = NA, color = 'black'))


# ggplot(data = SRs_nf_2d[,.(straightness = mean(s), t = as.numeric(.id.1),
# 			   diffusion = mean(d), distance = mean(emax)), 
# 			by = .(traj_length = round(l/100)*100)],
#        aes(t, traj_length, fill = distance, color = distance)) + geom_tile()




ggplot(data = SRs_nf_2d[,.(straightness = mean(s),
			   diffusion = mean(d)), 
			by = .(traj_length = round(l/100)*100, t = .id.1)][as.numeric(t) > 5],
       aes(traj_length/ 10, diffusion, color = factor(t)), alpha = 0.5) + 
	# aes(log(traj_length/ 10), straightness, color = factor(.id)), alpha = 0.5) + 
	geom_point(size = 4, shape = 21) +
	geom_line()+
	scale_y_continuous('Diffusion',  breaks = seq(0, 300, 75))+
	scale_x_continuous('Travelled distance (cm)')+
	scale_color_manual('Time (min)', values = c('mediumpurple', 'gold3', 'brown4', 'grey', 'darkblue', 'black'),
			   labels = seq(5, 30, 5))+
	theme(legend.title = element_blank(), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.position = c(1/2, 4/5))















sq <- seq(5*120, 30*120, 5*120)
SRs_det_2d <- lapply(det, function(p){
	result <- lapply(sq, function(t){
		data <- setDT(p@data)[Frame <= t]
		dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
		set(data, j = 'd', value = as.numeric(dmatrix))
		SR <- unique(data[d > maxd, N_ind])
		SR <- unique(data[!N_ind %in% SR, N_ind])
		data_SR <- data[N_ind %in% SR]
		rbindlist(lapply(seq_along(SR), function(i){
			sbst <- data[N_ind == SR[i]]
			xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
			if(nrow(xy) > 3){
				Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
				maxL <- max(Ls)
				if(maxL >= 100){
					idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
				} else {
					idxs <- nrow(xy)
				}
				straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
				diffusion <- vapply(idxs, function(l) TrajDistance(xy[seq_len(l)]), numeric(1))
				distance <- vapply(idxs, function(l) TrajEmax(xy[seq_len(l)]), numeric(1))
				data.frame(l = Ls[idxs], s = straightness, d = diffusion, emax = distance)
			}
			
			
		}), idcol = TRUE)
	})
	names(result) <- sq/120
	rbindlist(result, idcol = TRUE)
	
})
names(SRs_det_2d) <- paste0('det_', seq_along(det))
SRs_det_2d <- rbindlist(SRs_det_2d, idcol = TRUE)
SRs_det_2d[['tag']] <- apply(SRs_det_2d[, 1:3], 1, paste, collapse = '_')

# ggplot(data = SRs_det_2d[,.(straightness = mean(s), t = as.numeric(.id.1),
# 			   diffusion = mean(d)), 
# 			by = .(traj_length = round(l/100)*100)],
#        aes(t, traj_length, fill = straightness, color = straightness)) + geom_tile()

ggplot(data = SRs_det_2d[,.(straightness = mean(s), t = as.numeric(.id.1),
			   diffusion = mean(d)), 
			by = .(traj_length = round(l/100)*100)][t > 5],
       aes(t, traj_length /10, fill = diffusion)) + geom_tile()+
	scale_fill_viridis_c('Diffusion',option = 'C', breaks = seq(0, 300, 50))+
	scale_x_continuous('Time (min)', breaks = seq(5, 30, 5))+
	scale_y_continuous('Travelled distance (cm)', breaks = seq(0, 600, 50))+
	theme(legend.position = c(0.125, 4/5), 
	      legend.background = element_rect(fill = NA, color = 'black'))


# ggplot(data = SRs_det_2d[,.(straightness = mean(s), t = as.numeric(.id.1),
# 			   diffusion = mean(d), distance = mean(emax)), 
# 			by = .(traj_length = round(l/100)*100)],
#        aes(t, traj_length, fill = distance, color = distance)) + geom_tile()




ggplot(data = SRs_det_2d[,.(straightness = mean(s),
			   diffusion = mean(d)), 
			by = .(traj_length = round(l/100)*100, t = .id.1)][as.numeric(t) > 5],
       aes(traj_length/ 10, diffusion, color = factor(t)), alpha = 0.5) + 
	# aes(log(traj_length/ 10), straightness, color = factor(.id)), alpha = 0.5) + 
	geom_point(size = 4, shape = 21) +
	geom_line()+
	scale_y_continuous('Diffusion',  breaks = seq(0, 500, 75))+
	scale_x_continuous('Travelled distance (cm)')+
	scale_color_manual('Time (min)', values = c('mediumpurple', 'gold3', 'brown4', 'grey', 'darkblue', 'black'),
			   labels = seq(5, 30, 5))+
	theme(legend.title = element_blank(), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.position = c(1/2, 4/5))



sq <- seq(5*120, 30*120, 5*120)
LRs_det_2d <- lapply(det, function(p){
	result <- lapply(sq, function(t){
		data <- setDT(p@data)[Frame <= t]
		dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
		set(data, j = 'd', value = as.numeric(dmatrix))
		LR <- unique(data[d > maxd, N_ind])
		data_LR <- data[N_ind %in% LR]
		rbindlist(lapply(seq_along(LR), function(i){
			sbst <- data[N_ind == LR[i]]
			xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
			if(nrow(xy) > 3){
				Ls <- vapply(2:nrow(xy), function(x) TrajLength(xy, endIndex = x), numeric(1))
				maxL <- max(Ls)
				if(maxL >= 100){
					idxs <- vapply(seq(100, max(Ls), 100), function(x) which.min(abs(Ls - x)), numeric(1))
				} else {
					idxs <- nrow(xy)
				}
				straightness <- vapply(idxs, function(l) TrajStraightness(xy[seq_len(l)]), numeric(1))
				diffusion <- vapply(idxs, function(l) TrajDistance(xy[seq_len(l)]), numeric(1))
				distance <- vapply(idxs, function(l) TrajEmax(xy[seq_len(l)]), numeric(1))
				data.frame(l = Ls[idxs], s = straightness, d = diffusion, emax = distance)
			}
			
			
		}), idcol = TRUE)
	})
	names(result) <- sq/120
	rbindlist(result, idcol = TRUE)
	
})
names(LRs_det_2d) <- paste0('det_', seq_along(det))
LRs_det_2d <- rbindlist(LRs_det_2d, idcol = TRUE)
LRs_det_2d[['tag']] <- apply(LRs_det_2d[, 1:3], 1, paste, collapse = '_')

# ggplot(data = LRs_det_2d[,.(straightness = mean(s), t = as.numeric(.id.1),
# 			   diffusion = mean(d)), 
# 			by = .(traj_length = round(l/100)*100)],
#        aes(t, traj_length, fill = straightness, color = straightness)) + geom_tile()

ggplot(data = LRs_det_2d[,.(straightness = mean(s), t = as.numeric(.id.1),
			    diffusion = mean(d)), 
			 by = .(traj_length = round(l/100)*100)][t > 5],
       aes(t, traj_length /10, fill = diffusion)) + geom_tile()+
	scale_fill_viridis_c('Diffusion',option = 'C', breaks = seq(0, 3000, 150))+
	scale_x_continuous('Time (min)', breaks = seq(5, 30, 5))+
	scale_y_continuous('Travelled distance (cm)', breaks = seq(0, 6000, 300))+
	theme(legend.position = c(0.125, 4/5), 
	      legend.background = element_rect(fill = NA, color = 'black'))


# ggplot(data = LRs_det_2d[,.(straightness = mean(s), t = as.numeric(.id.1),
# 			   diffusion = mean(d), distance = mean(emax)), 
# 			by = .(traj_length = round(l/100)*100)],
#        aes(t, traj_length, fill = distance, color = distance)) + geom_tile()




ggplot(data = LRs_det_2d[,.(straightness = mean(s),
			    diffusion = mean(d)), 
			 by = .(traj_length = round(l/100)*100, t = .id.1)][as.numeric(t) > 5],
       aes(traj_length/ 10, diffusion, color = factor(t)), alpha = 0.5) + 
	# aes(log(traj_length/ 10), straightness, color = factor(.id)), alpha = 0.5) + 
	geom_point(size = 4, shape = 21) +
	geom_line()+
	scale_y_continuous('Diffusion',  breaks = seq(0, 3000, 150))+
	scale_x_continuous('Travelled distance (cm)')+
	scale_color_manual('Time (min)', values = c('mediumpurple', 'gold3', 'brown4', 'grey', 'darkblue', 'black'),
			   labels = seq(5, 30, 5))+
	theme(legend.title = element_blank(), 
	      legend.background = element_rect(color = 'black', fill = NA), 
	      legend.position = c(0.15, 4/5))





plot(do.call('rbind' ,LRs), pch = 19, ylim = c(0, 1))
lines(do.call('rbind' ,LRs))

points(do.call('rbind' ,SRs), pch = 19, col = 'red')
lines(do.call('rbind' ,SRs), col = 'red')

xy <- TrajFromCoords(data_LR[, c('Xmm', 'Ymm', 'Frame')], fps = 2)


plot(seq(100, max(Ls), 100), straightness)



	
	result <- as.data.frame(t(vapply(seq_along(inds), function(ii){
		sbst <- data[N_ind == inds[ii]]
		xy <- TrajFromCoords(sbst[, c('Xmm', 'Ymm', 'Frame')], fps = 2)
		if(nrow(xy) > min_pos){
			c(TrajStraightness(xy), TrajDistance(xy), TrajLength(xy), 
			  mean(as.numeric(Mod(TrajVelocity(xy))), na.rm = TRUE), 
			  mean(as.numeric(Mod(TrajAcceleration(xy))), na.rm = TRUE),
			  nrow(xy) / 120, nrow(sbst)/120)
		} else {rep(0, 7)}
		
	}, numeric(7))))
	colnames(result) <- c('Straightness', 'Diffusion', 'Distance', 'Mean_v', 'Mean_acc', 'Time', 'Total_time')
	result[result[['Time']] > 0, ]
})
