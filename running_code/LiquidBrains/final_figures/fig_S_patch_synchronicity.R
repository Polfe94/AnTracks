source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/Simulation.R')
source('~/research/gits/AnTracks/src/fit_functions.R')
# load('~/research/gits/AnTracks/data/nf.RData')
# load('~/research/gits/AnTracks/data/det.RData')

library(parallel)

path <- '/home/polfer/research/gits/AutomatAnts/results/2024/both_no_listen/'
files <- list.files(path)
food_files <- files[grepl('food', files)]

parallel_foraging <- function(path){
	files <- list.files(path)
	food_files <- files[grepl('food', files)]
	
	l <- mclapply(seq_along(food_files), function(i){
		txt <- data.table(read_parquet(paste0(path, food_files[i])))
		splt <- strsplit(food_files[i], split = '_')
		rho <- round(as.numeric(splt[[1]][[2]]), 2)
		eps <- round(as.numeric(splt[[1]][[4]]), 2)
		tleft <- txt[1:6, min(t, na.rm = TRUE)]
		tright <- txt[7:12, min(t, na.rm = TRUE)]
		if(tleft < tright){
			p1 <- 1:6
			p2 <- 7:12
		} else {
			p1 <- 7:12
			p2 <- 1:6
		}
		overlap <- max(txt[p1, t])-min(txt[p2, t])
		p2_duration <- max(txt[p2, t]) + min(txt[p1, t])
		sync <- sign(overlap)
		data.frame(overlap = overlap / p2_duration, sync = sync, rho = rho, eps = eps)
	}, mc.cores = 4L)
	
	rbindlist(l)[order(rho, eps)]
}

result <- parallel_foraging(path)

overlap_hd <- enhance_land_res(land = result[sync > 0, .(overlap = median(overlap, na.rm = TRUE)), 
					     by = c('rho', 'eps')],
		 xvar = 'rho', yvar = 'eps', zvar = 'overlap', spar = 0.75)

sync_hd <- enhance_land_res(land = result[sync > 0, .(sync = sum(sync, na.rm = TRUE)), 
					  by = c('rho', 'eps')],
						 xvar = 'rho', yvar = 'eps',
			    zvar = 'sync', spar = 0.75)



ggplot(data = overlap_hd,
       aes(rho, eps, fill = overlap / 60)) +
	geom_raster() + geom_contour(linetype = 4, aes(z = overlap), color = 'grey5') +
	scale_fill_viridis(TeX('Synchronicity duration (min)'), 
			   option = 'C')+
	xlab('Proportion of scouts') + ylab ('Social copying')+
	theme(aspect.ratio = 1, legend.position = 'top', 
			    legend.key.width = unit(1.4, 'cm'))+
	guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))


ggplot(data = sync_hd,
       aes(rho, eps, fill = sync)) +
	geom_raster() + geom_contour(linetype = 4, aes(z = sync), color = 'grey5') +
	scale_fill_viridis(TeX('Patch synchronicity (%)'), 
			   option = 'C')+
	xlab('Proportion of scouts') + ylab ('Social copying')+
	theme(aspect.ratio = 1, legend.position = 'top', 
	      legend.key.width = unit(1.4, 'cm'))+
	guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))
