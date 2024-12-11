source('~/research/gits/AnTracks/src/Experiment.R')
library(magick)
library(png)
library(grid)

ant <- image_read('G:/research/gits/AnTracks/plots/figura_metodes/subplot_elements/ant_image.png')

rotated_ant_1 <- image_rotate(ant, 217)
rotated_ant_2 <- image_rotate(ant, 37)
g_ant1 <- rasterGrob(image_transparent(rotated_ant_1, 'white'), interpolate = T)
g_ant2 <- rasterGrob(image_transparent(rotated_ant_2, 'white'), interpolate = T)

node_1 <- c(560, 533, 532, 559, 558, 531)
node_2 <- c(505, 478, 479, 452, 453, 426)

png('C:/Users/POL/homing_move.png', 8000, 8000, res = 800)
draw_hexagons(size = 1.2, color = 'grey50') + 
	scale_x_continuous(limits = c(660, 1000)) +
	scale_y_continuous(limits = c(1340, 1700)) +
        
        ## arcs 1
        geom_arc(hex[hex$node == node_2[1], c('x', 'y')], r = 19, npoints = 300, size = 1.3, start = -1.17, end = -1.3,
                 fill = 'grey50', color = 'black', alpha = 0.3)+
        
        geom_arc(hex[hex$node == node_2[1], c('x', 'y')], r = 12, npoints = 300, size = 1.3, start = -0.5, end = -1.3,
                 fill = 'grey50', color = 'black', alpha = 0.9)+
        geom_arc(hex[hex$node == node_2[1], c('x', 'y')], r = 23, npoints = 300, size = 1.3, start =  -1.3, end = -1.83,
                 fill = 'grey50', color = 'black', alpha = 0.5)+
        
        geom_text(data = data.frame(x = c(785, 777, 808), y = c(1430, 1464, 1470),
                                    label = c('$\\alpha_{1}$', '$\\alpha_{2}$', '$\\alpha_{3}$')),
                  aes(x, y, label = lapply(label, TeX, output = 'character')), parse = TRUE,size = 13)+
        
        
        geom_arc(hex[hex$node == node_1[1], c('x', 'y')], r = 19, npoints = 300, size = 1.3, start = -0.83, end = -0.57,
                 fill = 'grey50', color = 'black', alpha = 0.3)+
        geom_arc(hex[hex$node == node_1[1], c('x', 'y')], r = 25, npoints = 300, size = 1.3, start = -0.57, end = -0.17,
                 fill = 'grey50', color = 'black', alpha = 0.5)+
        geom_arc(hex[hex$node == node_1[1], c('x', 'y')], r = 12, npoints = 300, size = 1.3, start = -1.5, end = -0.57,
                 fill = 'grey50', color = 'black', alpha = 0.9)+
        
        
        
        geom_text(data = data.frame(x = c(905, 870, 873), y = c(1509, 1550, 1518),
                                    label = c('$\\alpha_{1}$', '$\\alpha_{3}$', '$\\alpha_{2}$')),
                  aes(x, y, label = lapply(label, TeX, output = 'character')), parse = TRUE,size = 13)+
        
        # geom_arc(hex[hex$node == node_1[1], c('x', 'y')], r = 5, npoints = 300, size = 1.3, start = -0.5, end = -1.14,
        #          fill = NA, color = 'black')+
        # geom_arc(hex[hex$node == node_1[1], c('x', 'y')], r = 15, npoints = 300, size = 1.3, start =  -1.3, end = -1.83,
        #          fill = NA, color = 'black')+
        
        
	geom_polygon(data = data.frame(x = c(673, 673, 977, 977),
				       y = c(1364, 1691, 1691, 1364)),
		     fill = NA, aes(x, y), color = 'black', linetype = 5, size = 1.5)+
	
	annotation_custom(g_ant1, xmin = 830, xmax = 950, ymin = 1520, ymax = 1610)+
        annotation_custom(g_ant2, xmin = 830-86.5, xmax = 950-86.5, ymin = 1520-145, ymax = 1610-145)+
	theme(aspect.ratio = 1, plot.margin = margin(-1.5,-5,-3,-5, 'cm'))+

	geom_path(data = hex[hex$node %in% node_1, ][order(node_1), ], aes(x, y), color = 'grey30',
		  size = 4, linetype = 2)+
	geom_line(data = hex[hex$node %in% node_1[c(1, length(node_1))], ],
		  aes(x, y),
		  arrow = arrow(length = unit(0.5, 'cm'),ends = 'first', 
		  	      type = 'closed'), size = 1.2)+
        geom_path(data = hex[hex$node %in% node_2, ][rank(node_2), ], aes(x, y), color = 'grey30',
                  size = 4, linetype = 2)+
        geom_line(data = hex[hex$node %in% node_2[c(1, length(node_2))], ],
                  aes(x, y),
                  arrow = arrow(length = unit(0.5, 'cm'),ends = 'first', 
                                type = 'closed'), size = 1.2)

	
dev.off()
