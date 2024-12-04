source('~/research/2022/ANTS/AnTracks/src/Experiment.R')
library(magick)
library(png)
library(grid)

ant_scout <- image_read('C:/Users/POL/Downloads/ant_scout.png')
ant_recruit <- image_read('C:/Users/POL/Downloads/ant_recruit.png')

rotated_scout <- image_rotate(ant_scout, 344)
g_scout <- rasterGrob(image_transparent(rotated_scout, 'white'), interpolate = T)
rotated_recruit <- image_rotate(ant_recruit, 217)
g_recruit <- rasterGrob(image_transparent(rotated_recruit, 'white'), interpolate = T)

node <- c(533, 479, 507)

# png('/home/polfer/research/tesi/figures_ilustratives/rndm_move.png', 4000, 4000, res = 400)
png('~/research/gits/AnTracks/plots/rndm_move.png', 8000, 8000, res = 800)
draw_hexagons(size = 1.2, color = 'grey50', edges = compute_edges(hex[hex$y > 1300 & hex$y < 1700 & 
                                                                              hex$x > 600 & hex$x < 1050, ])) + 
        scale_x_continuous(limits = c(660, 1000)) +
        scale_y_continuous(limits = c(1340, 1700)) +
        geom_polygon(data = data.frame(x = c(673, 673, 977, 977),
                                       y = c(1364, 1691, 1691, 1364)),
                     fill = NA, aes(x, y), color = 'black', linetype = 5, size = 1.5)+
        
        ## scout
        
        annotation_custom(g_scout, xmin = 770, xmax = 860, ymin = 1480, ymax = 1585)+
        
        ## scout -- arrows
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
        
        ## scout -- probs
        geom_text(data = data.frame(x = c(792, 775, 835), 
                                    y =c(1573, 1532, 1531),
                                    label = c('$R_{s} = 0.3974$', '$L_{s} = 0.4185$', '$B_{s} = 0.1841$'),
                                    angle = c(90, 27, -28)),
                  aes(angle = angle, x = x, y = y,
                      label = lapply(label, TeX, output = 'character')), parse = TRUE,
                  size = 6) +
        
        
        ## recruit
        annotation_custom(g_recruit, xmin = 844, xmax = 936, ymin = 1350, ymax = 1455)+
        
        ## recruit -- arrows
        geom_line(data = data.frame(x = c(890.1, 890.1),
                                    y = c(1423, 1440)),
                  aes(x, y),
                  arrow = arrow(length = unit(0.5, 'cm'),
                                type = 'closed'), size = 1.2)+
        geom_line(data = data.frame(x = c(863.79, 846.79),
                                    y = c(1374.5, 1365)),
                  aes(x, y),
                  arrow = arrow(length = unit(0.5, 'cm'),ends = 'first',
                                type = 'closed'), size = 1.2)+
        geom_line(data = data.frame(x = c(916.4, 933.4),
                                    y = c(1374.5 ,1365)),
                  aes(x, y),
                  arrow = arrow(length = unit(0.5, 'cm'),
                                type = 'closed'), size = 1.2)+
        ## recruit -- probs
        
        geom_text(data = data.frame(x = c(862, 920, 875), 
                                    y =c(1383,1381 , 1426),
                                    label = c('$R_{r} = 0.3179$', '$L_{r} = 0.3783$', '$B_{r} = 0.3038$'),
                                    angle = c(28, -28, 90)),
                  aes(angle = angle, x = x, y = y,
                      label = lapply(label, TeX, output = 'character')), parse = TRUE,
                  size = 6) +
        theme(aspect.ratio = 1, plot.margin = margin(-1.5,-5,-3,-5, 'cm'))
dev.off()

