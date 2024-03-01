source('~/research/gits/AnTracks/src/Experiment.R')

### TABLE DIMENSIONS ###
setDT(hex)
subhex <- hex[y > 1000]
nudge <- -40

data <- data.frame(xmin = c(min(subhex[['x']])+nudge, min(subhex[['x']])),
		   xmax = c(min(subhex[['x']])+nudge, max(subhex[['x']])),
		   ymin = c(min(subhex[['y']]), min(subhex[['y']])+nudge),
		   ymax = c(max(subhex[['y']]), min(subhex[['y']])+nudge),
		   label = c(1, 1, 2, 2))
		   	
		   	
draw_hexagons() + 
	geom_segment(data = data, aes(x = xmin, xend = xmax,
				      y = ymin, yend = ymax, group = label),
		     arrow = arrow(length= unit(0.015, 'npc'), ends = 'both', type = 'closed'))+
	annotate('text', label = paste0(max(subhex[['y']])-min(subhex[['y']]), ' mm'), 
		 x = min(subhex[['x']])+nudge*2, 
		 y = (max(subhex[['y']])+min(subhex[['y']]))/2, angle = 90)+
	annotate('text', label = paste0(max(subhex[['x']])-min(subhex[['x']]), ' mm'), 
		 y = min(subhex[['y']])+nudge*2, 
		 x = (max(subhex[['x']])+min(subhex[['x']]))/2)


### EXP FOOD POSITIONS ###
load('~/research/gits/AnTracks/data/sto.RData')
set(subhex, j = 'label', value = seq_len(nrow(subhex)))
foodlist <- lapply(sto, function(i){
	x <- do.call('rbind', i@food)
	y <- apply(x[, c('x', 'y')], 1, get_node)
	vapply(y, function(ii){subhex[node == ii, label]}, integer(1))
})
foodlist <- foodlist[-10]

dt <- data.table(do.call('cbind', foodlist))
colnames(dt) <- paste0('Exp_', seq_len(ncol(dt)))

detect_times <- data.table(do.call('cbind', lapply(1:9, function(i){
	# relative to the first ant that left the nest
	do.call('rbind', sto[[i]]@food)[['t']]/2 - min(sto[[i]]@data[['Time_sec']])
	# absolute times
	# do.call('rbind', sto[[i]]@food)[['t']]/2 - min(sto[[i]]@data[['Time_sec']])
})))
colnames(detect_times) <- paste0('Exp_', seq_len(ncol(detect_times)))


### SAVES ###

write.table(dt, file = '/home/polfer/research/papers/MATHEMATICS_OF_MOVEMENT/food_nodes.csv',
		sep = ',', dec = '.', row.names = FALSE)
write.table(detect_times, 
	    file = '/home/polfer/research/papers/MATHEMATICS_OF_MOVEMENT/detection_times.csv',
	    sep = ',', dec = '.', row.names = FALSE)
write.table(subhex[, c('label', 'x', 'y')], 
	    file = '/home/polfer/research/papers/MATHEMATICS_OF_MOVEMENT/coords.csv',
	    sep = ',', dec = '.', row.names = FALSE)


png('/home/polfer/research/papers/MATHEMATICS_OF_MOVEMENT/lattice_representation.png',
    1920, 1080, res = 100)
draw_hexagons(alpha = 0.15) + 
	geom_segment(data = data, aes(x = xmin, xend = xmax,
				      y = ymin, yend = ymax, group = label),
		     arrow = arrow(length= unit(0.015, 'npc'), ends = 'both', type = 'closed'))+
	annotate('text', label = paste0(max(subhex[['y']])-min(subhex[['y']]), ' mm'), 
		 x = min(subhex[['x']])+nudge*2, 
		 y = (max(subhex[['y']])+min(subhex[['y']]))/2, angle = 90)+
	annotate('text', label = paste0(max(subhex[['x']])-min(subhex[['x']]), ' mm'), 
		 y = min(subhex[['y']])+nudge*2, 
		 x = (max(subhex[['x']])+min(subhex[['x']]))/2)+
	annotate('text', label = subhex[['label']], x= subhex[['x']], y=subhex[['y']])+
	geom_point(data = subhex[node == 634], aes(x, y-20), size = 4, shape = 24, fill = 'blue')
dev.off()


# ggarrange(
# 	plotlist = lapply(seq_along(foodlist), function(i){
# draw_hexagons() +
# 	geom_segment(data = data, aes(x = xmin, xend = xmax,
# 				      y = ymin, yend = ymax, group = label),
# 		     arrow = arrow(length= unit(0.015, 'npc'), ends = 'both', type = 'closed'))+
# 	annotate('text', label = paste0(max(subhex[['y']])-min(subhex[['y']]), ' mm'),
# 		 x = min(subhex[['x']])+nudge*2,
# 		 y = (max(subhex[['y']])+min(subhex[['y']]))/2, angle = 90)+
# 	annotate('text', label = paste0(max(subhex[['x']])-min(subhex[['x']]), ' mm'),
# 		 y = min(subhex[['y']])+nudge*2,
# 		 x = (max(subhex[['x']])+min(subhex[['x']]))/2)+
# 	annotate('text', label = subhex[label %in% foodlist[[i]], label],
# 		 x = subhex[label %in% foodlist[[i]], x],
# 		 y = subhex[label %in% foodlist[[i]], y])
# }),labels = seq_along(foodlist))

for(i in seq_along(foodlist)){
	png(paste0('/home/polfer/research/papers/MATHEMATICS_OF_MOVEMENT/food_distrib_',i,'.png'),
	    1920, 1080, res = 100)
	print(draw_hexagons(alpha = 0.15, add = geom_foodpatches(sto[[i]]@food, alpha = 0.3)) + 
		geom_segment(data = data, aes(x = xmin, xend = xmax,
					      y = ymin, yend = ymax, group = label),
			     arrow = arrow(length= unit(0.015, 'npc'), ends = 'both', type = 'closed'))+
		annotate('text', label = paste0(max(subhex[['y']])-min(subhex[['y']]), ' mm'), 
			 x = min(subhex[['x']])+nudge*2, 
			 y = (max(subhex[['y']])+min(subhex[['y']]))/2, angle = 90)+
		annotate('text', label = paste0(max(subhex[['x']])-min(subhex[['x']]), ' mm'), 
			 y = min(subhex[['y']])+nudge*2, 
			 x = (max(subhex[['x']])+min(subhex[['x']]))/2)+
		annotate('text', label = subhex[['label']], x= subhex[['x']], y=subhex[['y']])+
		geom_point(data = subhex[node == 634], aes(x, y-20), size = 4, shape = 24, fill = 'blue'))
	dev.off()
}
png('/home/polfer/research/papers/MATHEMATICS_OF_MOVEMENT/lattice_representation.png',
    1920, 1080, res = 100)
draw_hexagons(alpha = 0.15) + 
	geom_segment(data = data, aes(x = xmin, xend = xmax,
				      y = ymin, yend = ymax, group = label),
		     arrow = arrow(length= unit(0.015, 'npc'), ends = 'both', type = 'closed'))+
	annotate('text', label = paste0(max(subhex[['y']])-min(subhex[['y']]), ' mm'), 
		 x = min(subhex[['x']])+nudge*2, 
		 y = (max(subhex[['y']])+min(subhex[['y']]))/2, angle = 90)+
	annotate('text', label = paste0(max(subhex[['x']])-min(subhex[['x']]), ' mm'), 
		 y = min(subhex[['y']])+nudge*2, 
		 x = (max(subhex[['x']])+min(subhex[['x']]))/2)+
	annotate('text', label = subhex[['label']], x= subhex[['x']], y=subhex[['y']])+
	geom_point(data = subhex[node == 634], aes(x, y-20), size = 4, shape = 24, fill = 'blue')
dev.off()