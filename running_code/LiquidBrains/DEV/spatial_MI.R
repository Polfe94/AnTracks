#### LIBRARIES, DATA AND GENERIC FUNCTIONS ####
source('~/research/gits/AnTracks/src/Simulation.R')
source('~/research/gits/AnTracks/src/Experiment.R')
load('~/research/gits/AnTracks/data/det.RData')
load('~/research/gits/AnTracks/data/sto.RData')
load('~/research/gits/AnTracks/data/nf.RData')

t0 <- Sys.time()
x1 <- lapply(seq(20*120, (21600 - 20*120), 300), function(i){
	mutual_info(det[[1]], i:(i-1+20*120), edges)
})
Sys.time() - t0

path_uniform <- '/home/polfer/research/gits/AutomatAnts/results/without_recruitment/parameters/uniform/'
uniform <- data.table(read.csv(paste0(path_uniform, list.files(path_uniform)[3])))
food <- data.table(read.csv(paste0(path_uniform, list.files(path_uniform)[1])))
l <- new('Simulation', data = uniform)

y1 <- lapply(seq(1, (21600 - 20*120), 300), function(i){
	mutual_info(l, i:(i-1+20*120), sim_edges)
})

# y1 <- list()
for(i in seq(1, (21600 - 20*120), 300)){
	y1[[length(y1) + 1]] <- mutual_info(l, i:(i-1+20*120), sim_edges)
}
draw_hexagons(edges = sim_edges, z = rank(y1[[64]]), linewidth = 3) + 
	scale_color_viridis_c(option = 'magma')
draw_hexagons(edges = sim_edges, z = rank(apply(do.call('rbind', y1), 2, mean)),
	      linewidth = 3, lineend = 'round', 
	      add = draw_hexagons(sim_edges, linewidth = 3.5, lineend = 'round',
	      		    add = geom_foodpatches(list(GP1 = data.frame(hex_sim[hex_sim$node %in% food[1:6, node], ]),
	      		    				   GP2 = data.frame(hex_sim[hex_sim$node %in% food[7:12, node], ])),
	      		    				   fill = 'mediumpurple1'))) +
	scico::scale_color_scico(palette = 'bilbao')

