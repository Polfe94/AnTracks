source('G:/research/2022/AnTracks/src/Experiment.R')
load('G:/research/2022/AnTracks/data/det.RData')
load('G:/research/2022/AnTracks/data/sto.RData')


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

food_times <- apply(do.call('rbind', lapply(det, function(i){
        x <- do.call('rbind', i@food)
        range(x$t)
})), 2, mean)


det_df <- data.table(t = seq_along(det_N) / 120,
                     N = norm_range(det_N, a = 0),
                     I = norm_range(cumsum(det_I), a = 0))
det_df <- reshape2::melt(det_df, id.vars = 't')
det_ylim <- c(30 / max(det_N), 2000 / max(cumsum(det_I)))

# png('G:/research/2022/AnTracks/plots/spin_glasses_paper/a.png', width = 1920, height = 1080, res = 200)
a <- ggplot(data = det_df, aes(t, value, color = variable)) +
        geom_food(t = food_times / 120, ylim = c(0, 1),
                  rectangle_params = list(color = 'black', fill = 'lightgrey', alpha = 0.75,
                                          linetype = 2))+
        geom_line(size = 2)+
        scale_color_viridis_d('', end = 0.6, labels = c('N', 'CI'))+
        scale_x_continuous('Time (min)', breaks = seq(0, 180, 20), limits = c(0, 180))+
        scale_y_continuous('Activity (N)',
                           breaks = seq(0, det_ylim[1], length.out = 7),
                           labels = seq(0, 30, length.out = 7),
                           
                           sec.axis = sec_axis(trans = ~.*1,
                                               breaks = seq(0, det_ylim[2], length.out = 9),
                                               labels = seq(0, 2000, length.out = 9),
                                               name = 'Cumulated interactions (CI)'),
                           
                           limits = c(0, max(det_ylim)))+
        ggtitle('DET')
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


p1 <- ggarrange(
draw_FoodPatches(det[[1]], fill = 'black', alpha = 0.2, add = 
draw_hexagons(z = colMeans(do.call('rbind', det_phases$p1)), size = 1.5, alpha = 0.9, add = draw_hexagons(size = 1)) +
        theme(aspect.ratio = 0.5, legend.position = 'none', plot.margin = unit(c(-20,10,-260,-20), 'pt')) +
        scico::scale_color_scico(palette = 'davos', direction = -1)
),
draw_FoodPatches(det[[1]], fill = 'black', alpha = 0.2, add = 
                         draw_hexagons(z = colMeans(do.call('rbind', det_phases$p2)), size = 1.5, alpha = 0.9, add = draw_hexagons(size = 1)) +
                         theme(aspect.ratio = 0.5, legend.position = 'none',plot.margin = unit(c(-20,10,-260,-20), 'pt')) +
                         scico::scale_color_scico(palette = 'davos', direction = -1)
),
draw_FoodPatches(det[[1]], fill = 'black', alpha = 0.2, add = 
                         draw_hexagons(z = colMeans(do.call('rbind', det_phases$p3)), size = 1.5, alpha = 0.97, add = draw_hexagons(size = 1)) +
                         theme(aspect.ratio = 0.5, legend.position = 'none', plot.margin = unit(c(-20,10,-260,-20), 'pt')) +
                         scico::scale_color_scico(palette = 'davos', direction = -1)
), ncol = 3, nrow = 1
)

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

food_times <- apply(do.call('rbind', lapply(sto, function(i){
        x <- do.call('rbind', i@food)
        range(x$t)
})), 2, mean)


sto_df <- data.table(t = seq_along(sto_N) / 120,
                     N = norm_range(sto_N, a = 0),
                     I = norm_range(cumsum(sto_I), a = 0))
sto_df <- reshape2::melt(sto_df, id.vars = 't')
sto_ylim <- c(35 / max(sto_N), 1000 / max(cumsum(sto_I)))

# png('G:/research/2022/AnTracks/plots/spin_glasses_paper/a.png', width = 1920, height = 1080, res = 200)
b <- ggplot(data = sto_df, aes(t, value, color = variable)) +
        geom_food(t = food_times / 120, ylim = c(0, 1),
                  rectangle_params = list(color = 'black', fill = 'lightgrey', alpha = 0.75,
                                          linetype = 2))+
        geom_line(size = 2)+
        scale_color_viridis_d('', end = 0.6, labels = c('N', 'CI'))+
        scale_x_continuous('Time (min)', breaks = seq(0, 180, 20), limits = c(0, 180))+
        scale_y_continuous('Activity (N)',
                           breaks = seq(0, sto_ylim[1], length.out = 8),
                           labels = seq(0, 35, length.out = 8),
                           
                           sec.axis = sec_axis(trans = ~.*1,
                                               breaks = seq(0, sto_ylim[2], length.out = 9),
                                               labels = seq(0, 1000, length.out = 9),
                                               name = 'Cumulated interactions (CI)'),
                           
                           limits = c(0, max(sto_ylim)))+
        ggtitle('STO')
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
draw_nodes(
        add = draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p1)), size = 2, add = draw_hexagons(size = 2.5)) +
        theme(aspect.ratio = 0.5) + scico::scale_color_scico(palette = 'davos', direction = -1),
        fill = 'white', size = 1
)

# draw_nodes(fill = 'white', color = 'black', size = 0.95, add = 
# draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p1)), size = 1.5, alpha = 0.95,
#              add = draw_hexagons(size = 1))) +
#         theme(aspect.ratio = 0.5) + scico::scale_color_scico(palette = 'davos', direction = -1)
draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p1)), size = 1.5, alpha = 0.9, add = draw_hexagons(size = 1)) +
        theme(aspect.ratio = 0.5) + scico::scale_color_scico(palette = 'davos', direction = -1)
draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p2)), size = 1.5, alpha = 0.9, add = draw_hexagons(size = 1)) +
        theme(aspect.ratio = 0.5) + scico::scale_color_scico(palette = 'davos', direction = -1)
draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p3)), size = 1.5, alpha = 0.9, add = draw_hexagons(size = 1)) +
        theme(aspect.ratio = 0.5) + scico::scale_color_scico(palette = 'davos', direction = -1)

sto_food <- lapply(sto, function(i) i@food)

p2 <- ggarrange(
        geom_foodpatches(sto_food, fill = 'black', alpha = 0.2, add = 
                                 draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p1)), size = 1.5, alpha = 0.9, add = draw_hexagons(size = 1)) +
                                 theme(aspect.ratio = 0.5, legend.position = 'none', plot.margin = unit(c(-260,10,0,-20), 'pt')) +
                                 scico::scale_color_scico(palette = 'davos', direction = -1)
        ),
        geom_foodpatches(sto_food, fill = 'black', alpha = 0.2, add = 
                                 draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p2)), size = 1.5, alpha = 0.9, add = draw_hexagons(size = 1)) +
                                 theme(aspect.ratio = 0.5, legend.position = 'none',plot.margin = unit(c(-260,10,0,-20), 'pt')) +
                                 scico::scale_color_scico(palette = 'davos', direction = -1)
        ),
        geom_foodpatches(sto_food, fill = 'black', alpha = 0.2, add = 
                                 draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p3)), size = 1.5, alpha = 0.97, add = draw_hexagons(size = 1)) +
                                 theme(aspect.ratio = 0.5, legend.position = 'none', plot.margin = unit(c(-260,10,0,-20), 'pt')) +
                                 scico::scale_color_scico(palette = 'davos', direction = -1)
        ), ncol = 3, nrow = 1
)

png('G:/research/2022/AnTracks/plots/spin_glasses_paper/combined.png', width = 6000, height = 4000, res = 300)
ggarrange(p1, p2, nrow = 2)
dev.off()

png('G:/research/2022/AnTracks/plots/spin_glasses_paper/ab.png', width = 1920*2, height = 1080*1.5, res = 200*1.5)
ggarrange(a,b, nrow = 1, common.legend = T, legend = 'bottom')
dev.off()
# png('G:/research/2022/AnTracks/plots/spin_glasses_paper/det_p1.png', width = 640, height = 480, res = 80)
# draw_hexagons(z = colMeans(do.call('rbind', det_phases$p1)), size = 2, 
#               add = draw_FoodPatches(det[[1]], alpha = 0.5, add = draw_hexagons(size = 2.5))) +
#         theme(aspect.ratio = 0.5) + scico::scale_color_scico(palette = 'davos', direction = -1)
# dev.off()
# 
# draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p2)), size = 2, add = draw_hexagons(size = 3))
# draw_hexagons(z = colMeans(do.call('rbind', sto_phases$p3)), size = 2, add = draw_hexagons(size = 3))
