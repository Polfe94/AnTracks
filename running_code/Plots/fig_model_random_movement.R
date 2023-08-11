source('~/research/2022/ANTS/AnTracks/src/Experiment.R')
library(magick)
library(png)
library(grid)

ant <- image_read('~/research/tesi/figures_ilustratives/ant_image.png')
rotated_ant <- image_rotate(ant, 344)
g_rotated <- rasterGrob(image_transparent(rotated_ant, 'white'), interpolate = T)

node <- c(533, 479, 507)

png('/home/polfer/research/tesi/figures_ilustratives/rndm_move.png', 4000, 4000, res = 400)
draw_hexagons(size = 1.2, color = 'grey50') + 
	scale_x_continuous(limits = c(660, 1000)) +
	scale_y_continuous(limits = c(1340, 1700)) +
	geom_polygon(data = data.frame(x = c(673, 673, 977, 977),
				       y = c(1364, 1691, 1691, 1364)),
		     fill = NA, aes(x, y), color = 'black', linetype = 5, size = 1.5)+
	
	annotation_custom(g_rotated, xmin = 770, xmax = 860, ymin = 1480, ymax = 1585)+
	geom_line(data = data.frame(x = c(846.79, 830),
				    y = c(1515, 1524.302)),
		  aes(x, y),
		  arrow = arrow(length = unit(0.5, 'cm'),
		  	      type = 'closed'), size = 1.2)+
	geom_line(data = data.frame(x = c(776.98, 760.19),
				    y = c(1524.694, 1515)),
		  aes(x, y),
		  arrow = arrow(length = unit(0.5, 'cm'),ends = 'first',
		  	      type = 'closed'), size = 1.2)+
	geom_line(data = data.frame(x = c(803.79, 803.79),
				    y = c(1570 ,1590)),
		  aes(x, y),
		  arrow = arrow(length = unit(0.5, 'cm'),
		  	      type = 'closed'), size = 1.2)+
	theme(aspect.ratio = 1, plot.margin = margin(-1.5,-5,-3,-5, 'cm'))
dev.off()

