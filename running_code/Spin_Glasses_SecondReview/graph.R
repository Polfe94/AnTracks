source('~/research/gits/AnTracks/src/Experiment.R')
source('~/research/gits/AnTracks/src/visualization.R')

png('~/research/papers/SPIN_GLASSES/review_2/neighbours/graph.png', 1920, 1080, res = 100)
draw_nodes(add = draw_hexagons(linewidth = 1.1), size = 12, fill = 'white') +
	annotate('text', x = hex[["x"]], y = hex[["y"]], label = hex[["node"]]) + 
	scale_y_continuous(limits = c(1000, 2000)) + theme_void()
dev.off()

edges <- compute_edges()
sbst <- hex[['node']][hex[['y']] > 1000]
neighbors <- vapply(sbst, function(i){
	paste0(get_neighbors(i), collapse = ',')
}, character(1))
network <- data.table(data.frame(node = sbst, neighbors = neighbors))

write.table(network, 
	    file = "/home/polfer/research/papers/SPIN_GLASSES/review_2/neighbours/network.csv",
	    sep = ',', dec = '.', row.names = F)
