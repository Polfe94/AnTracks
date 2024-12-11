source('~/research/2022/ANTS/AnTracks/src/Experiment.R')
library(magick)
library(png)
library(grid)

# ant <- image_read('G:/research/gits/AnTracks/plots/figura_metodes/subplot_elements/ant_image.png')

ant_scout <- image_read('C:/Users/POL/Downloads/ant_scout.png')
ant_recruit <- image_read('C:/Users/POL/Downloads/ant_recruit.png')

rotated_scout <- image_rotate(ant_scout, 344)
g_scout <- rasterGrob(image_transparent(rotated_scout, 'white'), interpolate = T)
rotated_recruit <- image_rotate(ant_recruit, 217)
g_recruit <- rasterGrob(image_transparent(rotated_recruit, 'white'), interpolate = T)

l_scout <- rasterGrob(image_transparent(image_rotate(ant_scout, 127), 'white'), interpolate = T)
l_recruit <- rasterGrob(image_transparent(image_rotate(ant_recruit, 127), 'white'), interpolate = T)

l_ant <- rasterGrob(image_transparent(image_rotate(ant, 127), 'white'), interpolate = T)


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
        
        
        ## legend
        geom_rect(fill = 'white', color = 'black', 
                  aes(xmin = 870-18, xmax = 985-30, ymax = 1680, ymin = 1617))+
        annotation_custom(l_scout, xmin = 770 + 60 + 5, xmax = 860 + 60 - 5,
                          ymin = 1480 + 130 +5, ymax = 1585 + 130 -5)+
        annotation_custom(l_recruit, xmin = 770 + 60 + 5, xmax = 860 + 60 - 5,
                          ymin = 1480 + 100 +5, ymax = 1585 + 100 -5)+
        
        geom_text(data = data.frame(x = c(919.7, 925), 
                                    y =c(1663, 1634),
                                    label = c('Scout', 'Recruit')),
                  aes(x = x, y = y,
                      label = lapply(label, TeX, output = 'character')), parse = TRUE,
                  size = 11) +
        
        ## scout
        annotation_custom(g_scout, xmin = 765, xmax = 865, ymin = 1475, ymax = 1590)+
        
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
        geom_text(data = data.frame(x = c(794, 768, 839), 
                                    y =c(1580, 1534, 1533),
                                    label = c('R=0.3974', 'L=0.4185', 'B=0.1841'),
                                    angle = c(90, 27, -28)),
                  aes(angle = angle, x = x, y = y,
                      label = lapply(label, TeX, output = 'character')), parse = TRUE,
                  size = 9) +
        
        
        ## recruit
        annotation_custom(g_recruit, xmin = 839, xmax = 941, ymin = 1345, ymax = 1460)+
        
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
        
        geom_text(data = data.frame(x = c(853, 923, 880), 
                                    y =c(1385,1383 , 1437),
                                    label = c('R=0.3179', 'L=0.3783', 'B=0.3038'),
                                    angle = c(28, -28, 90)),
                  aes(angle = angle, x = x, y = y,
                      label = lapply(label, TeX, output = 'character')), parse = TRUE,
                  size = 9) +
        theme(aspect.ratio = 1, plot.margin = margin(-1.5,-5,-3,-5, 'cm'))
dev.off()

