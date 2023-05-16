source('~/research/2022/ANTS/AnTracks/src/Experiment.R')
load('~/research/2022/ANTS/AnTracks/data/det.RData')
load('~/research/2022/ANTS/AnTracks/data/sto.RData')


#### --- DETERMINIST EXPERIMENTS --- ####
det <- lapply(det, get_N)
det <- lapply(det, get_I)
mx <- max(sapply(det, function(i) length(i@N)))
det_N <- lapply(det, function(i){
        x <- rep(NA, mx)
        x[seq_along(i@N)] <- i@N
        x
})
det_N <- apply(do.call('rbind', det_N), 2, mean, na.rm = T)

det_I <- lapply(det, function(i){
        x <- rep(NA, mx)
        x[seq_along(i@I)] <- i@I
        x
})
det_I <- apply(do.call('rbind', det_I), 2, mean, na.rm = T)

# food_times <- apply(do.call('rbind', lapply(det, function(i){
#         x <- do.call('rbind', i@food)
#         range(x$t)
# })), 2, mean)

det_food_times <- apply(do.call('rbind', lapply(seq_along(det)[-7], function(i){
        x <- do.call('rbind', det[[i]]@food)
        range(x$t)
})), 2, mean)


det_df <- data.table(t = seq_along(det_N) / 120,
                     N = norm_range(det_N, a = 0),
                     I = norm_range(cumsum(det_I), a = 0))
det_df <- reshape2::melt(det_df, id.vars = 't')
det_ylim <- c(30 / max(det_N), 2000 / max(cumsum(det_I)))

# png('G:/research/2022/AnTracks/plots/spin_glasses_paper/a.png', width = 1920, height = 1080, res = 200)
a <- ggplot(data = det_df, aes(t, value, color = variable)) +
        geom_food(t = det_food_times / 120, ylim = c(0, 1),
                  rectangle_params = list(color = 'black', fill = 'lightgrey', alpha = 0.75,
                                          linetype = 2))+
        geom_line(size = 2)+
        annotate('text', x = 2,y = 1.05,label = 'Exploration', size = 4.75)+
        annotate('text', x = 28,y = 1.05,label = 'Recruitment', size = 4.75)+
        annotate('text', x = 110,y = 1.05,label = 'Post-recruitment', size = 4.75)+
        scale_color_viridis_d('', end = 0.6, labels = c('Activity', 'Interactions'))+
        scale_x_continuous('Time (min)', breaks = seq(0, 180, 20), limits = c(0, 180))+
        scale_y_continuous('Activity',
                           breaks = seq(0, det_ylim[1], length.out = 7),
                           labels = seq(0, 30, length.out = 7),
                           
                           sec.axis = sec_axis(trans = ~.*1,
                                               breaks = seq(0, det_ylim[2], length.out = 9),
                                               labels = seq(0, 2000, length.out = 9),
                                               name = 'Interactions'),
                           
                           limits = c(0, max(det_ylim)))+
        theme(legend.position = c(0.9, 0.6),
              legend.background = element_rect(fill = 'white', size = 1.2,
                                               linetype = 'solid', color = 'black'),
              legend.title = element_blank())
        
# dev.off()


det_phases <- list(p1 = vector('list', length(det)),
                   p2 = vector('list', length(det)),
                   p3 = vector('list', length(det)))

t0 <- Sys.time()
for(i in seq_along(det)){
        t <- range(do.call('rbind', det[[i]]@food)$t)
        det_phases$p1[[i]] <- mutual_info(det[[i]], t = 1:t[1])
        det_phases$p2[[i]] <- mutual_info(det[[i]], t = (t[1]+1):t[2])
        det_phases$p3[[i]] <- mutual_info(det[[i]], t = (t[2]+1):length(det[[i]]@N))
}
Sys.time() - t0
beepr::beep(3)


# p1 <- ggarrange(
# draw_FoodPatches(det[[1]], fill = 'black', alpha = 0.2, add = 
# draw_hexagons(z = colMeans(do.call('rbind', det_phases$p1)), size = 1.5, alpha = 0.9, add = draw_hexagons(size = 1)) +
#         theme(aspect.ratio = 0.5, legend.position = 'none', plot.margin = unit(c(-20,10,-260,-20), 'pt')) +
#         scico::scale_color_scico(palette = 'davos', direction = -1)
# ),
# draw_FoodPatches(det[[1]], fill = 'black', alpha = 0.2, add = 
#                          draw_hexagons(z = colMeans(do.call('rbind', det_phases$p2)), size = 1.5, alpha = 0.9, add = draw_hexagons(size = 1)) +
#                          theme(aspect.ratio = 0.5, legend.position = 'none',plot.margin = unit(c(-20,10,-260,-20), 'pt')) +
#                          scico::scale_color_scico(palette = 'davos', direction = -1)
# ),
# draw_FoodPatches(det[[1]], fill = 'black', alpha = 0.2, add = 
#                          draw_hexagons(z = colMeans(do.call('rbind', det_phases$p3)), size = 1.5, alpha = 0.97, add = draw_hexagons(size = 1)) +
#                          theme(aspect.ratio = 0.5, legend.position = 'none', plot.margin = unit(c(-20,10,-260,-20), 'pt')) +
#                          scico::scale_color_scico(palette = 'davos', direction = -1)
# ), ncol = 3, nrow = 1
# )


zdet <- list(z1 = apply(do.call('rbind', det_phases$p1), 2, mean),
             z2 = apply(do.call('rbind', det_phases$p2), 2, mean),
             z3 = apply(do.call('rbind', det_phases$p3), 2, mean))

zdetvec <- rank(c(unlist(zdet), recursive = T))
zdet <- list(z1 = zdetvec[1:(length(zdetvec)/3)],
             z2 = zdetvec[(length(zdetvec)/3 + 1):(2*length(zdetvec)/3)],
             z3 = zdetvec[(2*length(zdetvec)/3 + 1):length(zdetvec)])

# p1 <- ggarrange(
#         draw_hexagons(z = zdet$z1,
#                       size = 3, lineend = 'round',
#                       add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add = 
#                                                   geom_foodpatches(food = det[[1]]@food, fill = 'mediumpurple1', color = 'white'))) +
#                 theme(aspect.ratio = 0.5, legend.position = 'none', plot.margin = unit(c(-20,10,-260,-20), 'pt')) +
#                 scico::scale_color_scico(palette = 'bilbao', direction = 1),
#                 
#         draw_hexagons(z = zdet$z2,
#                       size = 3, lineend = 'round',
#                       add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add = 
#                                                   geom_foodpatches(food = det[[1]]@food, fill = 'mediumpurple1', color = 'white'))) +
#                 theme(aspect.ratio = 0.5, legend.position = 'none', plot.margin = unit(c(-20,10,-260,-20), 'pt')) +
#                 scico::scale_color_scico(palette = 'bilbao', direction = 1),
#         draw_hexagons(z = zdet$z3,
#                       size = 3, lineend = 'round',
#                       add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add = 
#                                                   geom_foodpatches(food = det[[1]]@food, fill = 'mediumpurple1',
#                                                                    color = 'white'))) +
#                 theme(aspect.ratio = 0.5, legend.position = 'none', plot.margin = unit(c(-20,10,-260,-20), 'pt')) +
#                 scico::scale_color_scico(palette = 'bilbao', direction = 1),
#         
# nrow = 1)

#### --- STOCHASTIC EXPERIMENTS --- ####
sto <- lapply(sto, get_N)
sto <- lapply(sto, get_I)
mx <- max(sapply(sto, function(i) length(i@N)))
sto_N <- lapply(sto, function(i){
        x <- rep(NA, mx)
        x[seq_along(i@N)] <- i@N
        x
})
sto_N <- apply(do.call('rbind', sto_N), 2, mean, na.rm = T)

sto_I <- lapply(sto, function(i){
        x <- rep(NA, mx)
        x[seq_along(i@I)] <- i@I
        x
})
sto_I <- apply(do.call('rbind', sto_I), 2, mean, na.rm = T)

# food_times <- apply(do.call('rbind', lapply(sto, function(i){
#         x <- do.call('rbind', i@food)
#         range(x$t)
# })), 2, mean)

sto_food_times <- apply(do.call('rbind', lapply(seq_along(sto)[-2], function(i){
        x <- do.call('rbind', sto[[i]]@food)
        range(x$t)
})), 2, mean)


sto_df <- data.table(t = seq_along(sto_N) / 120,
                     N = norm_range(sto_N, a = 0),
                     I = norm_range(cumsum(sto_I), a = 0))
sto_df <- reshape2::melt(sto_df, id.vars = 't')
sto_ylim <- c(35 / max(sto_N), 1000 / max(cumsum(sto_I)))

# png('G:/research/2022/AnTracks/plots/spin_glasses_paper/a.png', width = 1920, height = 1080, res = 200)
b <- ggplot(data = sto_df, aes(t, value, color = variable)) +
        geom_food(t = sto_food_times / 120, ylim = c(0, 1),
                  rectangle_params = list(color = 'black', fill = 'lightgrey', alpha = 0.75,
                                          linetype = 2))+
        geom_line(size = 2)+
        annotate('text', x = 5,y = 1.05,label = 'Exploration', size = 4.75)+
        annotate('text', x = 31,y = 1.05,label = 'Recruitment', size = 4.75)+
        annotate('text', x = 115,y = 1.05,label = 'Post-recruitment', size = 4.75)+
        scale_color_viridis_d('', end = 0.6, labels = c('Activity', 'Interactions'))+
        scale_x_continuous('Time (min)', breaks = seq(0, 180, 20), limits = c(0, 180))+
        scale_y_continuous('Activity',
                           breaks = seq(0, sto_ylim[1], length.out = 8),
                           labels = seq(0, 35, length.out = 8),
                           
                           sec.axis = sec_axis(trans = ~.*1,
                                               breaks = seq(0, sto_ylim[2], length.out = 9),
                                               labels = seq(0, 1000, length.out = 9),
                                               name = 'Interactions'),
                           
                           limits = c(0, max(sto_ylim)))+
        theme(legend.position = c(0.9, 0.6),
              legend.background = element_rect(fill = 'white', size = 1.2,
                                               linetype = 'solid', color = 'black'),
                legend.title = element_blank())
# dev.off()


sto_phases <- list(p1 = vector('list', length(sto)),
                   p2 = vector('list', length(sto)),
                   p3 = vector('list', length(sto)))

t0 <- Sys.time()
for(i in seq_along(sto)){
        t <- range(do.call('rbind', sto[[i]]@food)$t)
        sto_phases$p1[[i]] <- mutual_info(sto[[i]], t = 1:t[1])
        sto_phases$p2[[i]] <- mutual_info(sto[[i]], t = (t[1]+1):t[2])
        sto_phases$p3[[i]] <- mutual_info(sto[[i]], t = (t[2]+1):length(sto[[i]]@N))
}
Sys.time() - t0
beepr::beep(3)


# mod_edges[, c('x', 'xend')] <- mod_edges[, c('x', 'xend')] + 0.1
# draw_nodes(
#         add = draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p1)), size = 2, add = draw_hexagons(size = 2.5)) +
#         theme(aspect.ratio = 0.5) + scico::scale_color_scico(palette = 'davos', direction = -1),
#         fill = 'white', size = 1
# )

# draw_nodes(fill = 'white', color = 'black', size = 0.95, add = 
# draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p1)), size = 1.5, alpha = 0.95,
#              add = draw_hexagons(size = 1))) +
#         theme(aspect.ratio = 0.5) + scico::scale_color_scico(palette = 'davos', direction = -1)
# draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p1)), size = 1.5, alpha = 0.9, add = draw_hexagons(size = 1)) +
#         theme(aspect.ratio = 0.5) + scico::scale_color_scico(palette = 'davos', direction = -1)
# draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p2)), size = 1.5, alpha = 0.9, add = draw_hexagons(size = 1)) +
#         theme(aspect.ratio = 0.5) + scico::scale_color_scico(palette = 'davos', direction = -1)
# draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p3)), size = 1.5, alpha = 0.9, add = draw_hexagons(size = 1)) +
#         theme(aspect.ratio = 0.5) + scico::scale_color_scico(palette = 'davos', direction = -1)
# 
# sto_food <- lapply(sto, function(i) i@food)
# 
# 
# draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p2)), size = 3, lineend = 'round',
#               add = draw_hexagons(size = 3, color = 'black', lineend = 'round', add = 
#                                           geom_foodpatches(food = sto_food, fill = 'grey30', color = 'white'))) +
#         theme(aspect.ratio = 0.5, legend.position = 'none',plot.margin = unit(c(-260,10,0,-20), 'pt')) +
#         scico::scale_color_scico(palette = 'davos', direction = -1)
# geom_foodpatches(sto_food, fill = 'black', alpha = 0.2, add = 
#                          draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p2)), size = 3, lineend = 'round',
#                                        add = draw_hexagons(size = 3, color = 'black', lineend = 'round')) +
#                          theme(aspect.ratio = 0.5, legend.position = 'none',plot.margin = unit(c(-260,10,0,-20), 'pt')) +
#                          scico::scale_color_scico(palette = 'davos', direction = -1)
# )
# 
# zsto <- list(z1 = apply(do.call('rbind', sto_phases$p1[-2]), 2, mean),
#              z2 = apply(do.call('rbind', sto_phases$p2[-2]), 2, mean),
#              z3 = apply(do.call('rbind', sto_phases$p3[-2]), 2, mean))
# 
# zstovec <- rank(c(unlist(zsto), recursive = T))
# zsto <- list(z1 = zstovec[1:(length(zstovec)/3)],
#              z2 = zstovec[(length(zstovec)/3 + 1):(2*length(zstovec)/3)],
#              z3 = zstovec[(2*length(zstovec)/3 + 1):length(zstovec)])



# p2 <- ggarrange(
#         draw_hexagons(z = zsto$z1, size = 3, lineend = 'round',
#                       add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add = 
#                                                   geom_foodpatches(food = sto_food, fill = 'mediumpurple1', color = 'white'))) +
#                 theme(aspect.ratio = 0.5, legend.position = 'none', plot.margin = unit(c(-260,10,0,-20), 'pt')) +
#                 scico::scale_color_scico(palette = 'bilbao'),
#         draw_hexagons(z = zsto$z2, size = 3, lineend = 'round',
#                       add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add = 
#                                                   geom_foodpatches(food = sto_food, fill = 'mediumpurple1', color = 'white'))) +
#                 theme(aspect.ratio = 0.5, legend.position = 'none', plot.margin = unit(c(-260,10,0,-20), 'pt')) +
#                 scico::scale_color_scico(palette = 'bilbao'),
#         draw_hexagons(z = zsto$z3, size = 3, lineend = 'round',
#                       add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add = 
#                                                   geom_foodpatches(food = sto_food, fill = 'mediumpurple1', color = 'white'))) +
#                 theme(aspect.ratio = 0.5, legend.position = 'none', plot.margin = unit(c(-260,10,0,-20), 'pt')) +
#                 scico::scale_color_scico(palette = 'bilbao'),
# nrow = 1)

# p2 <- ggarrange(
#         geom_foodpatches(sto_food, fill = 'black', alpha = 0.2, add = 
#                                  draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p1)), size = 1.5, alpha = 0.9, add = draw_hexagons(size = 1)) +
#                                  theme(aspect.ratio = 0.5, legend.position = 'none', plot.margin = unit(c(-260,10,0,-20), 'pt')) +
#                                  scico::scale_color_scico(palette = 'davos', direction = -1)
#         ),
#         geom_foodpatches(sto_food, fill = 'black', alpha = 0.2, add = 
#                                  draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p2)), size = 1.5, alpha = 0.9, add = draw_hexagons(size = 1)) +
#                                  theme(aspect.ratio = 0.5, legend.position = 'none',plot.margin = unit(c(-260,10,0,-20), 'pt')) +
#                                  scico::scale_color_scico(palette = 'davos', direction = -1)
#         ),
#         geom_foodpatches(sto_food, fill = 'black', alpha = 0.2, add = 
#                                  draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p3)), size = 1.5, alpha = 0.97, add = draw_hexagons(size = 1)) +
#                                  theme(aspect.ratio = 0.5, legend.position = 'none', plot.margin = unit(c(-260,10,0,-20), 'pt')) +
#                                  scico::scale_color_scico(palette = 'davos', direction = -1)
#         ), ncol = 3, nrow = 1
# )

# png('G:/research/2022/AnTracks/plots/spin_glasses_paper/combined.png', width = 6000, height = 4000, res = 300)
# ggarrange(p1, p2, nrow = 2)
# dev.off()

# png('~/research/2022/ANTS/AnTracks/plots/combined.png', width = 6000, height = 4000, res = 300)
# ggarrange(p1, p2, nrow = 2)
# dev.off()


zdet <- list(z1 = apply(do.call('rbind', det_phases$p1[-7]), 2, mean),
             z2 = apply(do.call('rbind', det_phases$p2[-7]), 2, mean),
             z3 = apply(do.call('rbind', det_phases$p3[-7]), 2, mean))
zsto <- list(z1 = apply(do.call('rbind', sto_phases$p1[-2]), 2, mean),
             z2 = apply(do.call('rbind', sto_phases$p2[-2]), 2, mean),
             z3 = apply(do.call('rbind', sto_phases$p3[-2]), 2, mean))

supervec <- rank(c(unlist(zdet), unlist(zsto), recursive = T))
ztotal <- list(d1 = supervec[1:(length(supervec)/6)],
               d2 = supervec[(length(supervec)/6 + 1):(2*length(supervec)/6)],
               d3 = supervec[(2*length(supervec)/6 + 1):(3*length(supervec)/6)],
               s1 = supervec[(3*length(supervec)/6 + 1):(4*length(supervec)/6)],
               s2 = supervec[(4*length(supervec)/6 + 1):(5*length(supervec)/6)],
               s3 = supervec[(5*length(supervec)/6 + 1):length(supervec)])


# png('~/research/2022/ANTS/AnTracks/plots/cd.png', width = 6000, height = 4000, res = 300)
# ggarrange( draw_hexagons(z = ztotal$d1,
#                          size = 3, lineend = 'round',
#                          add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', 
#                                              add = geom_foodpatches(
#                                                      food = det[[1]]@food, fill = 'mediumpurple1', color = 'white')
#                                              )) +
#                    ylab('DET')+ ggtitle('Scouting')+
#                    theme(aspect.ratio = 0.5, legend.position = 'none', 
#                          plot.margin = unit(c(-20,2.5,-250,10), 'pt'),
#                          plot.title = element_text(hjust = 0.5))+
#                    scico::scale_color_scico(palette = 'bilbao', direction = 1),
#            
#            draw_hexagons(z = ztotal$d2,
#                          size = 3, lineend = 'round',
#                          add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add = 
#                                                      geom_foodpatches(food = det[[1]]@food, fill = 'mediumpurple1', color = 'white'))) +
#                    theme(aspect.ratio = 0.5, legend.position = 'none', 
#                          plot.margin = unit(c(-20,2.5,-250,-25), 'pt'),
#                          plot.title = element_text(hjust = 0.5)) +
#                    ggtitle('Recruitment')+
#                    scico::scale_color_scico(palette = 'bilbao', direction = 1),
#            draw_hexagons(z = ztotal$d3,
#                          size = 3, lineend = 'round',
#                          add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add = 
#                                                      geom_foodpatches(food = det[[1]]@food, fill = 'mediumpurple1',
#                                                                       color = 'white'))) +
#                    ggtitle('Post-recruitment')+
#                    theme(aspect.ratio = 0.5, legend.position = 'none', 
#                          plot.margin = unit(c(-20,2.5,-250,-25), 'pt'),
#                          plot.title = element_text(hjust = 0.5)) +
#                    scico::scale_color_scico(palette = 'bilbao', direction = 1),
#         draw_hexagons(z = ztotal$s1, size = 3, lineend = 'round',
#                       add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add = 
#                                                   geom_foodpatches(food = sto_food, fill = 'mediumpurple1', color = 'white'))) +
#                 theme(aspect.ratio = 0.5, legend.position = 'none', plot.margin = unit(c(-250,2.5,0,10), 'pt')) +
#                 scico::scale_color_scico(palette = 'bilbao')+
#                 ylab('STO'),
#         draw_hexagons(z = ztotal$s2, size = 3, lineend = 'round',
#                       add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add = 
#                                                   geom_foodpatches(food = sto_food, fill = 'mediumpurple1', color = 'white'))) +
#                 theme(aspect.ratio = 0.5, legend.position = 'none', plot.margin = unit(c(-250,2.5,0,-25), 'pt')) +
#                 scico::scale_color_scico(palette = 'bilbao'),
#         draw_hexagons(z = ztotal$s3, size = 3, lineend = 'round',
#                       add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add = 
#                                                   geom_foodpatches(food = sto_food, fill = 'mediumpurple1', color = 'white'))) +
#                 theme(aspect.ratio = 0.5, legend.position = 'none', plot.margin = unit(c(-250,2.5,0,-25), 'pt')) +
#                 scico::scale_color_scico(palette = 'bilbao'),
#         nrow = 2, ncol = 3)
# dev.off()
# 
# png('G:/research/2022/AnTracks/plots/spin_glasses_paper/ab.png', width = 1920*2, height = 1080*1.5, res = 200*1.5)
# ggarrange(a,b, nrow = 1, common.legend = T, legend = 'bottom')
# dev.off()
# png('G:/research/2022/AnTracks/plots/spin_glasses_paper/det_p1.png', width = 640, height = 480, res = 80)
# draw_hexagons(z = colMeans(do.call('rbind', det_phases$p1)), size = 2, 
#               add = draw_FoodPatches(det[[1]], alpha = 0.5, add = draw_hexagons(size = 2.5))) +
#         theme(aspect.ratio = 0.5) + scico::scale_color_scico(palette = 'davos', direction = -1)
# dev.off()
# 
# draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p2)), size = 2, add = draw_hexagons(size = 3))
# draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p3)), size = 2, add = draw_hexagons(size = 3))




L2 <- list(draw_hexagons(z = ztotal$d1,
                         size = 3, lineend = 'round',
                         add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', 
                                             add = geom_foodpatches(
                                                     food = det[[1]]@food, fill = 'mediumpurple1', color = 'white')
                         )) +
                   ylab('DET')+ ggtitle('Exploration')+
                   theme(aspect.ratio = 0.5, legend.position = 'none',
                         plot.title = element_text(hjust = 0.5, size = 16))+
                   scico::scale_color_scico(palette = 'bilbao', direction = 1)+
                   geom_point(data = data.frame(x = hex$x[hex$node == 634],
                                                y = hex$y[hex$node == 634] - 35),
                              aes(x, y), shape = 24, fill = 'mediumpurple1', size = 3),
           
           draw_hexagons(z = ztotal$d2,
                         size = 3, lineend = 'round',
                         add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add = 
                                                     geom_foodpatches(food = det[[1]]@food, fill = 'mediumpurple1', color = 'white'))) +
                   theme(aspect.ratio = 0.5, legend.position = 'none', 
                         plot.title = element_text(hjust = 0.5, size = 16)) +
                   ggtitle('Recruitment')+
                   scico::scale_color_scico(palette = 'bilbao', direction = 1)+
                   geom_point(data = data.frame(x = hex$x[hex$node == 634],
                                                y = hex$y[hex$node == 634] - 35),
                              aes(x, y), shape = 24, fill = 'mediumpurple1', size = 3),
           draw_hexagons(z = ztotal$d3,
                         size = 3, lineend = 'round',
                         add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add = 
                                                     geom_foodpatches(food = det[[1]]@food, fill = 'mediumpurple1',
                                                                      color = 'white'))) +
                   ggtitle('Post-recruitment')+
                   theme(aspect.ratio = 0.5, legend.position = 'none',
                         plot.title = element_text(hjust = 0.5, size = 16)) +
                   scico::scale_color_scico(palette = 'bilbao', direction = 1)+
                   geom_point(data = data.frame(x = hex$x[hex$node == 634],
                                                y = hex$y[hex$node == 634] - 35),
                              aes(x, y), shape = 24, fill = 'mediumpurple1', size = 3),
           draw_hexagons(z = ztotal$s1, size = 3, lineend = 'round',
                         add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add = 
                                                     geom_foodpatches(food = sto_food, fill = 'mediumpurple1', color = 'white'))) +
                   theme(aspect.ratio = 0.5, legend.position = 'none') +
                   scico::scale_color_scico(palette = 'bilbao')+
                   ylab('STO')+
                   geom_point(data = data.frame(x = hex$x[hex$node == 634],
                                                y = hex$y[hex$node == 634] - 35),
                              aes(x, y), shape = 24, fill = 'mediumpurple1', size = 3),
           draw_hexagons(z = ztotal$s2, size = 3, lineend = 'round',
                         add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add = 
                                                     geom_foodpatches(food = sto_food, fill = 'mediumpurple1', color = 'white'))) +
                   theme(aspect.ratio = 0.5, legend.position = 'none') +
                   scico::scale_color_scico(palette = 'bilbao')+
                   geom_point(data = data.frame(x = hex$x[hex$node == 634],
                                                y = hex$y[hex$node == 634] - 35),
                              aes(x, y), shape = 24, fill = 'mediumpurple1', size = 3),
           draw_hexagons(z = ztotal$s3, size = 3, lineend = 'round',
                         add = draw_hexagons(size = 3.1, color = 'black', lineend = 'round', add = 
                                                     geom_foodpatches(food = sto_food, fill = 'mediumpurple1', color = 'white'))) +
                   theme(aspect.ratio = 0.5, legend.position = 'none') +
                   scico::scale_color_scico(palette = 'bilbao')+
                   geom_point(data = data.frame(x = hex$x[hex$node == 634],
                                                y = hex$y[hex$node == 634] - 35),
                              aes(x, y), shape = 24, fill = 'mediumpurple1', size = 3)
)




png('~/research/2022/ANTS/AnTracks/plots/whole_plot.png', width = 8000, height = 6000, res = 400)
ggarrange(plotlist = list(ggarrange(a, b, nrow = 1, labels = c('DET', 'STO'),
                                    font.label = list(face = 'plain', size = 18),
                                    hjust = -5),
               ggarrange(plotlist = L2, nrow = 2, ncol = 3)), 
          nrow = 2)
dev.off()

