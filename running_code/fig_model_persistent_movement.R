source('~/research/2022/ANTS/AnTracks/src/Experiment.R')
library(magick)
library(png)
library(grid)

ant <- image_read('~/research/tesi/figures_ilustratives/ant_image.png')
rotated_ant <- image_rotate(ant, 217)
g_rotated <- rasterGrob(image_transparent(rotated_ant, 'white'), interpolate = T)

node <- c(560, 533, 532, 559, 558, 531)

png('/home/polfer/research/tesi/figures_ilustratives/homing_move.png', 4000, 4000, res = 400)
draw_hexagons(size = 1.2, color = 'grey50') + 
	scale_x_continuous(limits = c(660, 1000)) +
	scale_y_continuous(limits = c(1340, 1700)) +
	geom_polygon(data = data.frame(x = c(673, 673, 977, 977),
				       y = c(1364, 1691, 1691, 1364)),
		     fill = NA, aes(x, y), color = 'black', linetype = 5, size = 1.5)+
	
	annotation_custom(g_rotated, xmin = 830, xmax = 950, ymin = 1520, ymax = 1610)+
	theme(aspect.ratio = 1, plot.margin = margin(-1.5,-5,-3,-5, 'cm'))+

	geom_path(data = hex[hex$node %in% node, ][order(node), ], aes(x, y), color = 'purple',
		  size = 4, linetype = 2)+
	geom_line(data = hex[hex$node %in% node[c(1, length(node))], ],
		  aes(x, y),
		  arrow = arrow(length = unit(0.5, 'cm'),ends = 'first', 
		  	      type = 'closed'), size = 1.2)
	
dev.off()
