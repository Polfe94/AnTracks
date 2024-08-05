source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
load('~/research/gits/AnTracks/data/nf.RData')
load('~/research/gits/AnTracks/data/det.RData')
# load('~/research/gits/AnTracks/results/mov_rho_R.RData')


maxt <- max(sapply(det, function(i){
	min(rbindlist(i@food)[['t']]) - min(i@data[['Frame']])
})) # 3000
maxd <- 500 # 500

library(trajr)
library(arrow)
library(egg)
library(latex2exp)

R_seq <- seq(100, 600, 10)

DET_scouts <- rbindlist(lapply(det, function(p){
	data <- setDT(p@data)[Frame < maxt + min(Frame)]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	set(data, j = 'd', value = as.numeric(dmatrix))
	LR <- unique(data[d > maxd, N_ind])
	SR <- unique(data[!N_ind %in% LR, N_ind])
	data_LR <- data[N_ind %in% LR]
	data_SR <- data[N_ind %in% SR]
	result_LR <- rbindlist(lapply(seq_along(LR), function(i){
		sbst <- data[N_ind == LR[i]]
		p <- vapply(R_seq, function(x){
			nrow(sbst[d < x])
		}, numeric(1)) / nrow(sbst)
		data.frame(R = R_seq, p = p, 
			   tag = paste0('det_LR_',i), scout = 'LR', condition = 'det')
		
	}), idcol = TRUE)
	result_SR <- rbindlist(lapply(seq_along(SR), function(i){
		sbst <- data[N_ind == SR[i]]
		p <- vapply(R_seq, function(x){
			nrow(sbst[d < x])
		}, numeric(1)) / nrow(sbst)
		data.frame(R = R_seq, p = p, 
			   tag = paste0('det_SR_',i), scout = 'SR', condition = 'det')
		
	}), idcol = TRUE)
	rbind(result_LR, result_SR)
}), idcol = TRUE)
DET_scouts[['tag']] <- apply(DET_scouts[, c(1, 5)], 1, paste0, collapse = '_')


NFD_scouts <- rbindlist(lapply(nf[-c(1, 2)], function(p){
	data <- setDT(p@data)[Frame < maxt + min(Frame)]
	dmatrix <- pdist(as.matrix(data[, c('Xmm', 'Ymm')]), as.matrix(hex[hex$node == 634, c('x', 'y')]))
	set(data, j = 'd', value = as.numeric(dmatrix))
	LR <- unique(data[d > maxd, N_ind])
	SR <- unique(data[!N_ind %in% LR, N_ind])
	data_LR <- data[N_ind %in% LR]
	data_SR <- data[N_ind %in% SR]
	result_LR <- rbindlist(lapply(seq_along(LR), function(i){
		sbst <- data[N_ind == LR[i]]
		p <- vapply(R_seq, function(x){
			nrow(sbst[d < x])
		}, numeric(1)) / nrow(sbst)
		data.frame(R = R_seq, p = p, 
			   tag = paste0('nf_LR_',i), scout = 'LR', condition = 'nf')
		
	}), idcol = TRUE)
	result_SR <- rbindlist(lapply(seq_along(SR), function(i){
		sbst <- data[N_ind == SR[i]]
		p <- vapply(R_seq, function(x){
			nrow(sbst[d < x])
		}, numeric(1)) / nrow(sbst)
		data.frame(R = R_seq, p = p, 
			   tag = paste0('nf_SR_',i), scout = 'SR', condition = 'nf')
		
	}), idcol = TRUE)
	rbind(result_LR, result_SR)
}), idcol = TRUE)
NFD_scouts[['tag']] <- apply(NFD_scouts[, c(1, 5)], 1, paste0, collapse = '_')

# ggplot(data = DET_scouts[, c('R', 'p', 'tag')], aes(R, p, color = tag)) + 
# 	geom_point(show.legend = FALSE)

# ggplot(data = DET_scouts[order(R, p), c('R', 'p', 'tag', 'scout')], aes(R, p, color = scout))  + 
# 	geom_line(aes(group = tag), linewidth = 2)+
# 	geom_point(data = DET_scouts[order(R, p),.(p = mean(p)),by = c('R', 'scout')],
# 		   aes(shape = scout), size = 5, color = 'black')+
# 	geom_line(data = DET_scouts[order(R, p),.(p = mean(p)),by = c('R', 'scout')], 
# 		  color = 'black', aes(group = scout))+
# 	scale_color_manual('', values = c('gold3', 'brown4'))
# 
# ggplot(data = DET_scouts[order(R, p), c('R', 'p', 'tag', 'scout')], aes(R, p, group = tag, color = scout)) + 
# 	geom_point(size = 4) + geom_line(show.legend = FALSE)
# 
# ggplot(data = NFD_scouts[order(R, p), c('R', 'p', 'tag', 'scout')], aes(R, p, group = tag, color = scout)) + 
# 	geom_point(size = 4) + geom_line(show.legend = FALSE)


png('/home/polfer/research/gits/AnTracks/plots/figures_LiquidBrains/Fig_S1.png', 6000, 3000, res = 350)
ggplot(data = rbind(DET_scouts, NFD_scouts)[order(R, p), c('R', 'p', 'tag', 'scout', 'condition')],
       aes(R, 100*p, group = tag, color = scout, shape = scout)) + 
	geom_point(size = 4) + geom_line(show.legend = FALSE)+
	facet_wrap(~ factor(condition, labels = c('DET', 'NFD')))+
	scale_y_continuous('Site fidelity (%)')+
	scale_x_continuous('Radius (mm)', breaks = seq(0, 600, 150))+
	theme(legend.title = element_blank(), legend.background = element_rect(fill = NA, color = 'black'))
dev.off()

# ggplot(data = NFD_scouts[tag == 'nf_SR_1', c('R', 'p', 'tag')], aes(R, p, group = tag)) + 
# 	geom_point(show.legend = FALSE) + geom_line(show.legend = FALSE)
# 
# ggplot(data = NFD_scouts[, .(p = mean(p)), by = c('R', 'scout')], aes(R, p, color = scout)) + 
# 	geom_point() + geom_path()
# 
# ggplot(data = NFD_scouts[R < 600, c('R', 'p')], aes(R, rank(p))) + 
# 	geom_density_2d_filled()
# 
# 
# ggplot(data = NFD_scouts[R < 400, c('R', 'p')], aes(R, p)) +
# 	geom_tile(aes(fill = ..density..), stat = "density2d", contour = FALSE) +
# 	scale_fill_viridis_c() +  # You can change the color scale as needed
# 	labs(x = "R", y = "p", fill = "Density") +
# 	theme_minimal()

