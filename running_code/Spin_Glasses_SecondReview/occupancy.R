source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/visualization.R')
load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/sto.RData')


## DETERMINIST - GREEN COLONY ##
for(det_exp in seq_along(det)){
	M <- get_matrix(det[[det_exp]], data = 0)
	p <- colMeans(M)
	df <- data.frame(node = hex[['node']][hex[['y']] > 1000], p = 0)
	n <- as.numeric(names(p))
	idx <- sapply(n, function(i) which(i == df[['node']]))
	df[['p']][idx] <- p
	if(det_exp == 7){
		write.table(x = df, 
			    file = paste0('~/research/papers/SPIN_GLASSES/review_2/occupancy/discarded/det/', 
						  det[[det_exp]]@date, '_g.csv'),
			    sep = ',', dec = '.', row.names = F)
	} else {
		write.table(x = df, file = paste0('~/research/papers/SPIN_GLASSES/review_2/occupancy/det/data/', 
						  det[[det_exp]]@date, '_g.csv'),
			    sep = ',', dec = '.', row.names = F)
	}
}

fls <- list.files('~/research/papers/SPIN_GLASSES/review_2/occupancy/det/data/')
fls <- fls[grepl('.csv', fls)]
for(x in fls){
	df <- read.csv(paste0('~/research/papers/SPIN_GLASSES/review_2/occupancy/det/data/', x))
	png(paste0('~/research/papers/SPIN_GLASSES/review_2/occupancy/det/plots/', 
		   gsub('.csv', '.png',x)),
	    1920, 1080, res = 150)
	print(draw_nodes(z = rank(df[['p']]), show.legend = F, size = 7.5, 
			 add = draw_hexagons(add = draw_FoodPatches(det[[1]], 
			 					   alpha = 0.5, fill = 'purple'))) + 
			 	theme(aspect.ratio = 0.5)+
			 	scico::scale_fill_scico(palette = 'bilbao', direction = 1)+
			 	ggtitle(gsub('_g.csv', '', x)))
	dev.off()
}

## STOCHASTIC - GREEN COLONY ##
for(sto_exp in seq_along(sto)){
	M <- get_matrix(sto[[sto_exp]], data = 0)
	p <- colMeans(M)
	df <- data.frame(node = hex[['node']][hex[['y']] > 1000], p = 0)
	n <- as.numeric(names(p))
	idx <- sapply(n, function(i) which(i == df[['node']]))
	df[['p']][idx] <- p
	if(sto_exp == 2){
		write.table(x = df, 
			    file = paste0('~/research/papers/SPIN_GLASSES/review_2/occupancy/discarded/sto/', 
						  sto[[sto_exp]]@date, '_g.csv'),
			    sep = ',', dec = '.', row.names = F)
	} else {
		write.table(x = df, file = paste0('~/research/papers/SPIN_GLASSES/review_2/occupancy/sto/data/', 
						  sto[[sto_exp]]@date, '_g.csv'),
			    sep = ',', dec = '.', row.names = F)
	}

}

fls <- list.files('~/research/papers/SPIN_GLASSES/review_2/occupancy/sto/data/')
fls <- fls[grepl('.csv', fls)]
l <- 1
for(x in fls){
	df <- read.csv(paste0('~/research/papers/SPIN_GLASSES/review_2/occupancy/sto/data/', x))
	png(paste0('~/research/papers/SPIN_GLASSES/review_2/occupancy/sto/plots/', gsub('.csv', '.png',x)),
	    1920, 1080, res = 150)
	print(draw_nodes(z = rank(df[['p']]), show.legend = F, size = 7.5, 
			 add = draw_hexagons(add = draw_FoodPatches(sto[[l]], 
			 					   alpha = 0.5, fill = 'purple'))) + 
			 	theme(aspect.ratio = 0.5)+
			 	scico::scale_fill_scico(palette = 'bilbao', direction = 1)+
			 	ggtitle(gsub('_g.csv', '', x)))
	dev.off()
	l <- l+1
}
